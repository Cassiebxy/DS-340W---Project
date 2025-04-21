# ==== Load required libraries ====
library(data.table)
library(lubridate)
library(ggplot2)
library(dplyr)
library(reshape2)  

# ==== Load and preprocess mortality data ====
# import the csv document, and change the document path as needed
#dat <- read.csv("local_mortality.csv") # for window system
dat <- read.csv("~/Desktop/local_mortality.csv") # for Mac system
setDT(dat)
dat[, log_deaths := log(deaths + 1)] # Log-transform to stabilize variance
dat[, date := as.Date(paste(year, time, "01", sep = "-"))]
dat[, region_id := paste(country_name, local_unit_name, sep = "_")]

# ==== Reshape to matrix form ====
mortality_surface <- dcast(dat, date ~ region_id, value.var = "log_deaths", fun.aggregate = mean)
setDT(mortality_surface)
mort_mat <- as.matrix(mortality_surface[, -1, with = FALSE])
dim(mort_mat)
rownames(mort_mat) <- as.character(mortality_surface$date)

# ==== Filter valid date range ====
mortality_surface[, date := as.Date(date)]

mort_mat <- as.matrix(mortality_surface[, -1, with = FALSE])
rownames(mort_mat) <- as.character(mortality_surface$date)  
all_dates <- as.Date(rownames(mort_mat))                    

valid_indices <- which(!is.na(all_dates) & all_dates > as.Date("2000-01-01")) 
all_dates <- all_dates[valid_indices]
mort_mat <- mort_mat[valid_indices, , drop = FALSE]

# ==== Split into train/test ====
cutoff <- floor(0.8 * length(all_dates))
train_dates <- all_dates[1:cutoff]
test_dates <- all_dates[(cutoff + 1):length(all_dates)]

dim(mort_mat)      
length(train_dates)    
length(test_dates)     
range(train_dates)      
range(test_dates)     

mort_mat_train <- mort_mat[rownames(mort_mat) %in% as.character(train_dates), , drop = FALSE]
mort_mat_test  <- mort_mat[rownames(mort_mat) %in% as.character(test_dates),  , drop = FALSE]

# ==== Select valid regions ====
region_stats <- data.table(
  region = colnames(mort_mat_train),
  sd = apply(mort_mat_train, 2, sd, na.rm = TRUE),
  n_obs = apply(mort_mat_train, 2, function(x) sum(!is.na(x)))
)
valid_regions <- region_stats[sd > 0.05 & n_obs > 30, region]
length(valid_regions) 
mort_mat_train <- mort_mat_train[, valid_regions, drop = FALSE]
mort_mat_test <- mort_mat_test[, valid_regions, drop = FALSE]

# ==== Remove incomplete rows ====
valid_rows <- complete.cases(mort_mat_train)
mort_mat_train <- mort_mat_train[valid_rows, , drop = FALSE]

# ==== PCA ====
pca_fit <- prcomp(mort_mat_train, center = TRUE, scale. = FALSE)

# ==== OU Process simulation ====
simulate_ou <- function(beta_hist, theta, mu, sigma, h = 36, dt = 1) {
  beta_sim <- numeric(h)
  beta_sim[1] <- tail(beta_hist, 1)
  for (i in 2:h) {
    dW <- rnorm(1, mean = 0, sd = sqrt(dt))
    beta_sim[i] <- beta_sim[i - 1] + theta * (mu - beta_sim[i - 1]) * dt + sigma * dW
  }
  return(beta_sim)
}

# OU params
theta <- 0.3

# ==== 1PC MODEL ====
# Extract first PC from PCA
beta_train_1pc <- pca_fit$x[, 1]
phi_1pc <- pca_fit$rotation[, 1]
region_names <- rownames(pca_fit$rotation)

# Simulate future β(t) for 1PC
horizon_1pc <- nrow(mort_mat_test)
mu_1pc <- mean(beta_train_1pc)
sigma_1pc <- sd(beta_train_1pc)
beta_future_1pc <- simulate_ou(beta_train_1pc, theta, mu_1pc, sigma_1pc, h = horizon_1pc)

# === Verification ===
cat("Length of beta_future_1pc (should match test months):", length(beta_future_1pc), "\n")
cat("Mean used for simulation (mu_1pc):", round(mu_1pc, 4), "\n")
cat("SD used for simulation (sigma_1pc):", round(sigma_1pc, 4), "\n")

# Reconstruct using outer product
forecast_matrix_1pc <- outer(beta_future_1pc, phi_1pc)
colnames(forecast_matrix_1pc) <- region_names
rownames(forecast_matrix_1pc) <- rownames(mort_mat_test)

# Evaluate 1PC model
mse_1pc <- mean((forecast_matrix_1pc - mort_mat_test)^2, na.rm = TRUE)
rmse_1pc <- sqrt(mse_1pc)

cat("---- 1PC Model ----\n")
cat("MSE:", round(mse_1pc, 4), "\n")
cat("RMSE:", round(rmse_1pc, 4), "\n\n")


# ==== 2PC MODEL ====
# Extract first 2 PCs
beta_train_2pc <- pca_fit$x[, 1:2]
phi_2pc <- pca_fit$rotation[, 1:2]

# Simulate both PCs
horizon_2pc <- nrow(mort_mat_test)
mu1 <- mean(beta_train_2pc[, 1])
mu2 <- mean(beta_train_2pc[, 2])
sd1 <- sd(beta_train_2pc[, 1])
sd2 <- sd(beta_train_2pc[, 2])
beta1_future <- simulate_ou(beta_train_2pc[, 1], theta, mu1, sd1, h = horizon_2pc)
beta2_future <- simulate_ou(beta_train_2pc[, 2], theta, mu2, sd2, h = horizon_2pc)

# === Verification ===
cat("Length of beta1_future:", length(beta1_future), "\n")
cat("Length of beta2_future:", length(beta2_future), "\n")
cat("Mu1, SD1:", round(mu1, 4), round(sd1, 4), "\n")
cat("Mu2, SD2:", round(mu2, 4), round(sd2, 4), "\n")

# Combine both components
forecast_matrix_2pc <- outer(beta1_future, phi_2pc[, 1]) + outer(beta2_future, phi_2pc[, 2])
colnames(forecast_matrix_2pc) <- region_names
rownames(forecast_matrix_2pc) <- rownames(mort_mat_test)

# Evaluate 2PC model
mse_2pc <- mean((forecast_matrix_2pc - mort_mat_test)^2, na.rm = TRUE)
rmse_2pc <- sqrt(mse_2pc)

cat("---- 2PC Model ----\n")
cat("MSE:", round(mse_2pc, 4), "\n")
cat("RMSE:", round(rmse_2pc, 4), "\n\n")



# ==== Visualization ====
# Actual log mortality
actual_dt <- as.data.table(mort_mat_test)
actual_dt[, date := as.Date(rownames(mort_mat_test))]
actual_dt <- melt(actual_dt, id.vars = "date", variable.name = "region", value.name = "actual_log")

# Predictions (1PC)
pred_1pc_dt <- as.data.table(forecast_matrix_1pc)
pred_1pc_dt[, date := as.Date(rownames(forecast_matrix_1pc))]
pred_1pc_dt <- melt(pred_1pc_dt, id.vars = "date", variable.name = "region", value.name = "pred_log")
setDT(pred_1pc_dt)
pred_1pc_dt[, model := "1PC"]

# Predictions (2PC)
pred_2pc_dt <- as.data.table(forecast_matrix_2pc)
pred_2pc_dt[, date := as.Date(rownames(forecast_matrix_2pc))]
pred_2pc_dt <- melt(pred_2pc_dt, id.vars = "date", variable.name = "region", value.name = "pred_log")
setDT(pred_2pc_dt)
pred_2pc_dt[, model := "2PC"]

# Combine
pred_all_dt <- rbind(pred_1pc_dt, pred_2pc_dt)
merged_all <- merge(actual_dt, pred_all_dt, by = c("date", "region"))
merged_all

# === SELECT TOP REGIONS ===
top_regions <- names(sort(apply(mort_mat_test, 2, sd, na.rm = TRUE), decreasing = TRUE))[1:3]
plot_dt <- merged_all %>% filter(region %in% top_regions)


# === PLOT INDIVIDUALS ===

# 1PC
ggplot(plot_dt %>% filter(model == "1PC"), aes(x = date)) +
  geom_line(aes(y = actual_log, color = region), linewidth = 1.2) +
  geom_line(aes(y = pred_log, color = region), linetype = "dashed", linewidth = 1.2) +
  labs(
    title = "1PC: Predicted vs Actual log(Mortality)",
    subtitle = "Dashed = Prediction",
    y = "log(mortality)", x = "Date"
  ) +
  theme_minimal()


# 2PC
ggplot(plot_dt %>% filter(model == "2PC"), aes(x = date)) +
  geom_line(aes(y = actual_log, color = region), linewidth = 1.2) +
  geom_line(aes(y = pred_log, color = region), linetype = "dashed", linewidth = 1.2) +
  labs(
    title = "2PC: Predicted vs Actual log(Mortality)",
    subtitle = "Dashed = Prediction",
    y = "log(mortality)",
    x = "Date"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(face = "italic", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

# === COMBINED FACET PLOT ===
ggplot(plot_dt, aes(x = date, group = region, color = region)) +
  # Actual log mortality: solid line
  geom_line(aes(y = actual_log), linewidth = 1) +
  
  # Predicted log mortality: dashed line
  geom_line(aes(y = pred_log), linetype = "dashed", linewidth = 1) +
  
  # Facet by model
  facet_wrap(~model, ncol = 2) +
  
  # Labels and theme
  labs(
    title = "Predicted vs Actual Log-Mortality by Model",
    subtitle = "Solid = Actual, Dashed = Forecasted (Top 3 Regions)",
    x = "Date",
    y = "log(Mortality)",
    color = "Region"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(face = "italic", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 13)
  )

