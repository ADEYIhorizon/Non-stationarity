#import the necessary libraries
library(ismev)
library(DEoptim)
library(extRemes)
library(writexl)
library(dplyr)


H <- read.csv("./Data/Subdaily data/SSP245_sub_daily_update.csv")
#H <- H[, -which(names(H) == "Years")]

model_spec <- read.csv("model SSP585.csv")

level_prep <- function(col = H[[i+1]], H=H, model_spec = model_spec[i,2]){
  col <- col
  cov <- matrix(ncol=2,nrow=nrow(H))
  cov[,1] <- H$Years
  if (model_spec == 1){
    #Model 1
    GEV <- fevd(col, H, units="mm") 
  }
  else if (model_spec == 2){
    #Model II
    GEV <- fevd(col, H,location.fun = ~ cov[,1])
  }
  else if (model_spec == 3){
    #Model III
    GEV <- fevd(col,H,location.fun = ~ cov[,1], scale.fun = ~cov[,1])
  }
  else if (model_spec == 4){
    #Model IV
    GEV <- fevd(col,H,location.fun = ~ cov[,1],scale.fun = ~ cov[,1], 
                        use.phi = TRUE)
  }
  else if (model_spec == 5){
    #Model V
    GEV <- fevd(col,H,scale.fun = ~ cov[,1])
  }
  else {
    #Model VI
    GEV <- fevd(col,H, scale.fun = ~ cov[,1],use.phi = TRUE)
  }
  
  RP <- c(2,5,10,20,25,50,100)
  PTOT_MODEL_RCP8.5_5 <- erlevd(GEV, period=RP)
  quant <- c(quantile(PTOT_MODEL_RCP8.5_5[1,],0.95),quantile(PTOT_MODEL_RCP8.5_5[2,],0.95),quantile(PTOT_MODEL_RCP8.5_5[3,],0.95),quantile(PTOT_MODEL_RCP8.5_5[4,],0.95),quantile(PTOT_MODEL_RCP8.5_5[5,],0.95),quantile(PTOT_MODEL_RCP8.5_5[6,],0.95),quantile(PTOT_MODEL_RCP8.5_5[7,],0.95))

  return(quant)
}

#initialize dataframe
# Create a matrix filled with 1s
matrix_with_ones <- matrix(1, nrow = 1, ncol = 7)
delta_df <- data.frame(matrix_with_ones)

for (i in 1:48) {
  print(i)
  ret_lev <- level_prep(col = H[[i+1]], H=H, model_spec = model_spec[i,2])
  delta_df <- rbind(delta_df, ret_lev)
  print("done")
}

colnames(delta_df) <- c("2-year", "5-year", "10-year", "20-year", "25-year", "50-year", "100-year")

write_xlsx(delta_df, path = "SSP245_IDF_new2.xlsx")



cov[,1] <- H$Years







