#import the necessary libraries
library(ismev)
library(DEoptim)
library(extRemes)
library(writexl)
library(dplyr)

##RCP45

H_RCP <- read.csv("./Data/Subdaily data/RCP45_sub_daily_update.csv")
plot(H_RCP$H2, type = "l",  lwd = 1.5, cex.lab = 1.25,
     xlab = "Year", ylab = "Annual Maximum Precipitation")
cov <- matrix(ncol=2,nrow=nrow(H_RCP))
cov[,1] <- H_RCP$Years
cov[,2] <- cov[,1]^2

AICc <- function(model=GEV_H_I_1H){
  n <- length(H_RCP[[2]])
  k <- length(model$results$par)
  AIC <- 2*model$results$value+2*length(model$results$par)
  lhs <- (2*k*(k+1))/(n-k-1)
  return(AIC + lhs)
}


delta_fun <- function(col, H=H){
  col <- col
  cov[,1] <- H$Years
  #Model 1
  GEV_H_I_1H <- fevd(col, H, units="mm")
  
  #Model II
  GEV_H_II_1H <- fevd(col, H,location.fun = ~ cov[,1])
  
  #Model III
  GEV_H_III_1H <- fevd(col,H,location.fun = ~ cov[,1], scale.fun = ~cov[,1])
  
  #Model IV
  
  GEV_H_IV_1H <- fevd(col,H,location.fun = ~ cov[,1],scale.fun = ~ cov[,1], 
                      use.phi = TRUE)
  
  #Model V
  GEV_H_V_1H <- fevd(col,H,scale.fun = ~ cov[,1])
  
  #Model VI
  GEV_H_VI_1H <- fevd(col,H, scale.fun = ~ cov[,1],use.phi = TRUE)
  
  AICc_H_1H <-
    c(AICc(GEV_H_I_1H), AICc(GEV_H_II_1H), AICc(GEV_H_III_1H), AICc(GEV_H_IV_1H),
      AICc(GEV_H_V_1H), AICc(GEV_H_VI_1H))
  min_1 <- min(AICc_H_1H)
  
  #Delta
  Delta_H_1H <- AICc_H_1H-min_1
  
  Delta_H_1H
  
  return(Delta_H_1H)
}

lr_fun <- function(col, H=H){
  col <- col
  cov[,1] <- H$Years
  #Model 1
  GEV_H_I_1H <- fevd(col, H, units="mm")
  
  #Model II
  GEV_H_II_1H <- fevd(col, H,location.fun = ~ cov[,1])
  
  #Model III
  GEV_H_III_1H <- fevd(col,H,location.fun = ~ cov[,1], scale.fun = ~cov[,1])
  
  #Model IV
  
  GEV_H_IV_1H <- fevd(col,H,location.fun = ~ cov[,1],scale.fun = ~ cov[,1], 
                      use.phi = TRUE)
  
  #Model V
  GEV_H_V_1H <- fevd(col,H,scale.fun = ~ cov[,1])
  
  #Model VI
  GEV_H_VI_1H <- fevd(col,H, scale.fun = ~ cov[,1],use.phi = TRUE)
  
  LR_H_1H <- c(lr.test(GEV_H_I_1H,GEV_H_I_1H)$p.value,lr.test(GEV_H_I_1H,GEV_H_II_1H)$p.value,lr.test(GEV_H_I_1H,GEV_H_III_1H)$p.value, 
               lr.test(GEV_H_I_1H,GEV_H_IV_1H)$p.value,lr.test(GEV_H_I_1H,GEV_H_V_1H)$p.value,lr.test(GEV_H_I_1H,GEV_H_VI_1H)$p.value)
  return(LR_H_1H)
}

#initialize dataframe
# Create a matrix filled with 1s
matrix_with_ones <- matrix(1, nrow = 1, ncol = 6)
delta_df <- data.frame(matrix_with_ones)
lr_df <- data.frame(matrix_with_ones)


for (i in 2:49) {
  print(i)
  Delta_H_1H <- delta_fun(col= H_RCP[[i]], H_RCP)
  delta_df <- rbind(delta_df, Delta_H_1H)
  print("done")
}

for (i in 2:49) {
  LR_H_1H <- lr_fun(col= H_RCP[[i]], H_RCP)
  lr_df <- rbind(lr_df, LR_H_1H)
  print("done")
}

colnames(delta_df) <- c("Model1", "Model2", "Model3", "Model4", "Model5", "Model6")
colnames(lr_df) <- c("Model1", "Model2", "Model3", "Model4", "Model5", "Model6")

write_xlsx(delta_df, path = "Delta_SSP585.xlsx")
write_xlsx(lr_df, path = "LR Test SSP585.xlsx")

