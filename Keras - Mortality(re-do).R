#### Purpose: Fit autoencoders to HMD data
#### Author: Ronald Richman
#### License: MIT
#### Data: The data were sourced from the HMD by downloading the relevant text files

require(data.table)
require(dplyr)
require(ggplot2)
require(data.table)
require(reshape2)
require(HMDHFDplus)
require(gnm)


### Get England/Wales data from StMoMo
#dat = HMDHFDplus::readHMD("c:/r/Mx_1x1.txt") %>% data.table

## What I change:
library(data.table)
file_path <- "/Users/cathybao/Desktop/AI_in_Actuarial_Science-master/GBRCENW/STATS/Mx_1x1.txt"
dat <- fread(file_path, skip = 2, header = TRUE, sep = " ", strip.white = TRUE)

# Converts a character column (from the second to the fifth column) to a numeric type
dat[, `:=`(
  Age = as.numeric(Age),
  Female = as.numeric(Female),
  Male = as.numeric(Male),
  Total = as.numeric(Total))]


### transform data
dat[,logmx:=log(Male)]

scale_min_max = function(dat,dat_test)  {
  min_dat = min(dat)
  max_dat = max(dat)
  dat_scaled=(dat-min_dat)/(max_dat-min_dat)
  dat_scaled_test = (dat_test-min_dat)/(max_dat-min_dat)
  return(list(train = dat_scaled, test = dat_scaled_test, min = min_dat, max=max_dat))
}

scale_z = function(dat,dat_test)  {
  mean_dat = mean(dat)
  sd_dat = sd(dat)
  dat_scaled=(dat-mean_dat)/(sd_dat)
  dat_scaled_test = (dat_test-mean_dat)/(sd_dat)
  return(list(train = dat_scaled, test = dat_scaled_test, mean_dat = mean_dat, sd_dat=sd_dat))
}

dat = dat[Year>1949 & Age<100]
train = dat[Year < 2000]
test = dat[Year >= 2000]

scaled = scale_min_max(train$logmx, test$logmx)

train$mx_scale = scaled$train
test$mx_scale = scaled$test

train_rates = train %>% dcast.data.table(Year~Age, value.var = "mx_scale")

### Lee-Carter Baseline - Fit on Raw Rates
fit = gnm(Male~ as.factor(Age) + Mult(as.factor(Age), as.factor(Year), inst = 1) -1,family = poisson(link = "log"),data=train)

train[,pred_LC:=predict(fit, type="response")]

coefs = data.table(names = fit %>% coef %>% names, coef=fit %>% coef)
coefs[,row:=.I]
coefs[row %in% c(1:100),var:="ax"]
coefs[row %in% c(101:200),var:="bx"]
coefs[row %in% c(201:250),var:="k"]

ax =coefs[var == "ax"]$coef
bx =coefs[var == "bx"]$coef
k =coefs[var == "k"]$coef

c1 = mean(k)
c2 = sum(bx)
ax = ax+c1*bx
bx = bx/c2
k = (k-c1)*c2

forecast_k <- k %>% forecast::forecast(17)
k_forecast = forecast_k[[2]]

fitted = (ax+(bx)%*%t(k))
fitted_test = (ax+(bx)%*%t(k_forecast)) %>% melt
##
test$pred_LC <- rep(fitted_test$value %>% exp, length.out = nrow(test))

results = data.table(model = "Lee-Carter", MSE_OutOfSample = 
                       test[,sum((Male - pred_LC)^2)])

test[Year==2016]%>% ggplot(aes(x=Age, y = log(pred_LC)))+ geom_point(size = 0.5, alpha=0.5)+facet_wrap(~Year)+geom_line(aes(x=Age, y=log(Male)))
