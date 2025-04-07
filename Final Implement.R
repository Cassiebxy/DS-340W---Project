# Load necessary packages
library(data.table)
library(lubridate)
library(ggplot2)
library(dplyr)

# Read and process the mortality dataset
dat <- read.csv("~/Desktop/local_mortality.csv")
setDT(dat)
dat[, log_deaths := log(deaths + 1)]
dat[, date := as.Date(paste(year, time, "01", sep = "-"))]
dat[, region_id := paste(country_name, local_unit_name, sep = "_")]

# Reshape data into a mortality surface (wide format)
mortality_surface <- dcast(dat, date ~ region_id, value.var = "log_deaths", fun.aggregate = mean)
setDT(mortality_surface)
mort_mat <- as.matrix(mortality_surface[, -1, with = FALSE])
rownames(mort_mat) <- as.character(mortality_surface$date)

# Split into training and testing sets based on time
all_dates <- sort(unique(as.Date(rownames(mort_mat))))
cutoff <- floor(0.8 * length(all_dates))
train_dates <- all_dates[1:cutoff]
test_dates <- all_dates[(cutoff + 1):length(all_dates)]

mort_mat_train <- mort_mat[rownames(mort_mat) %in% as.character(train_dates), , drop = FALSE]
mort_mat_test  <- mort_mat[rownames(mort_mat) %in% as.character(test_dates),  , drop = FALSE]

# Select valid regions based on standard deviation and data completeness
region_stats <- data.table(
  region = colnames(mort_mat_train),
  sd = apply(mort_mat_train, 2, sd, na.rm = TRUE),
  n_obs = apply(mort_mat_train, 2, function(x) sum(!is.na(x)))
)
valid_regions <- region_stats[sd > 0.05 & n_obs > 30, region]
mort_mat_train <- mort_mat_train[, valid_regions, drop = FALSE]
mort_mat_test <- mort_mat_test[, valid_regions, drop = FALSE]

# Remove incomplete rows in training data
valid_rows <- complete.cases(mort_mat_train)
mort_mat_train <- mort_mat_train[valid_rows, , drop = FALSE]

# Perform PCA on training data
pca_fit <- prcomp(mort_mat_train, center = TRUE, scale. = FALSE)
beta_train <- pca_fit$x[, 1]
phi_vec <- pca_fit$rotation[, 1]
region_names <- names(phi_vec)

# Define the Ornstein-Uhlenbeck simulation function
simulate_ou <- function(beta_hist, theta, mu, sigma, h = 36, dt = 1) {
  beta_sim <- numeric(h)
  beta_sim[1] <- tail(beta_hist, 1)
  for (i in 2:h) {
    dW <- rnorm(1, mean = 0, sd = sqrt(dt))
    beta_sim[i] <- beta_sim[i - 1] + theta * (mu - beta_sim[i - 1]) * dt + sigma * dW
  }
  return(beta_sim)
}

# Set OU parameters
theta <- 0.3
mu <- mean(beta_train)
sigma <- sd(beta_train)
h <- nrow(mort_mat_test)

# Simulate beta_t for the testing period
beta_future <- simulate_ou(beta_train, theta, mu, sigma, h = h)

# Combine training and testing dates and beta values
dates_train <- as.Date(rownames(mort_mat_train))
dates_test <- as.Date(rownames(mort_mat_test))
dates_all <- c(dates_train, dates_test)
beta_all <- c(beta_train, beta_future)

# Reconstruct log mortality for training and testing sets
forecast_matrix_train <- outer(beta_train, phi_vec)
forecast_matrix_test <- outer(beta_future, phi_vec)
colnames(forecast_matrix_test) <- region_names
rownames(forecast_matrix_test) <- as.character(dates_test)

# Compute MSE and RMSE on the test set
log_actual_test <- mort_mat_test
log_predicted_test <- forecast_matrix_test
mse_mat <- (log_predicted_test - log_actual_test)^2
mse_value <- mean(mse_mat, na.rm = TRUE)
rmse_value <- sqrt(mse_value)

cat("Test MSE (log mortality):", round(mse_value, 4), "\n")
cat("Test RMSE (log mortality):", round(rmse_value, 4), "\n")

# Visualization: compare predicted and actual log mortality for selected regions
top_regions <- names(sort(apply(mort_mat_test, 2, sd, na.rm = TRUE), decreasing = TRUE))[1:3]

actual_dt <- as.data.table(log_actual_test)
actual_dt[, date := dates_test]
actual_dt <- melt(actual_dt, id.vars = "date", variable.name = "region", value.name = "actual_log")

pred_dt <- as.data.table(log_predicted_test)
pred_dt[, date := dates_test]
pred_dt <- melt(pred_dt, id.vars = "date", variable.name = "region", value.name = "pred_log")

# Merge predicted and actual data
merged_dt <- merge(actual_dt, pred_dt, by = c("date", "region"))

merged_dt %>%
  filter(region %in% top_regions) %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = actual_log, color = region), linewidth = 1.2) +
  geom_line(aes(y = pred_log, color = region), linetype = "dashed", linewidth = 1.2) +
  labs(
    title = "Predicted vs Actual log(Mortality)",
    subtitle = "Solid = Actual, Dashed = Predicted",
    x = "Date", y = "log(Mortality)",
    color = "Region"
  ) +
  scale_color_brewer(palette = "Set1") +  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(face = "italic", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  scale_y_continuous(limits = c(0, NA))  

