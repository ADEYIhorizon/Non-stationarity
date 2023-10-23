library(extRemes)
library(dplyr)
library(tibble)
library(tidyselect)
library(readr)

AMS_data_red <- read.csv("AMS_data_red.csv")
plot(AMS_data_red$H_1, type = "l",  lwd = 1.5, cex.lab = 1.25,
     xlab = "Year", ylab = "Annual Maximum Precipitation")
#model I
fit0 <- fevd(H_1, AMS_data_red, units="mm")

fit0
plot(fit0)
plot(fit0, "trace")
return.level(fit0)
return.level(fit0, do.ci=TRUE)
as <- ci(fit0, return.period = c(2, 5, 10, 25, 50, 100))


#model 2
# Extract intercept and slope
fit1 <- fevd(H_1, AMS_data_red, location.fun=~time, units="mm", rperiods=c(2, 5, 10, 25, 50, 100))
fit1
#erlevd(fit1, period = 20)
plot(fit1)
plot(fit1, "trace")

v <- make.qcov(fit1, vals=list(mu1=c(1, 2)))
return.level(fit1, return.period=c(2, 5, 10, 25, 50, 100))


#model 3
# Extract intercept and slope
time <- AMS_data_red$time
fit2 <- fevd(H_1, AMS_data_red, location.fun=~time, scale.fun=~time, units="mm",  use.phi=TRUE)
fit2
plot(fit2)
plot(fit2, "trace")
v <- make.qcov(fit1, vals=list(mu1=c(1, 2)))
return.level(fit2, return.period=c(2, 5, 10, 25, 50, 100))


#model 4
# Extract intercept and slope
time <- AMS_data_red$time
fit3 <- fevd(H_1, AMS_data_red, location.fun=~time, scale.fun=~exp(time), units="mm", use.phi=TRUE)
fit3
plot(fit3)
return.level(fit3)
plot(fit3, "trace")
return.level(fit3, return.period=c(2, 5, 10, 25, 50, 100))

