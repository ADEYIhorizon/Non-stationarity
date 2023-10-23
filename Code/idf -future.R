#import the necessary libraries
library(ismev)
library(DEoptim)
library(extRemes)
library(writexl)
library(dplyr)


H <- read.csv("./AMS_data.csv")
#H <- H[, -which(names(H) == "Years")]

model_spec <- read.csv("model SSP585.csv")

level_prep_hist <- function(col = H[[i]], H=H){
  col <- col
    #Model 1
  GEV <- fevd(col, H, units="mm") 
  RP <- c(2,5,10,20,25,50,100)
  PTOT_MODEL_RCP8.5_5 <- erlevd(GEV, period=RP)
  return(PTOT_MODEL_RCP8.5_5)
}

#initialize dataframe
# Create a matrix filled with 1s
matrix_with_ones <- matrix(1, nrow = 1, ncol = 7)
delta_df <- data.frame(matrix_with_ones)

for (i in 3:50) {
  print(i)
  ret_lev <- level_prep_hist(col = H[[i]], H=H)
  delta_df <- rbind(delta_df, ret_lev[,1])
  print("done")
}

colnames(delta_df) <- c("2", "5", "10", "20", "25", "50", "100")

write_xlsx(delta_df, path = "IDF_Historical1.xlsx")



cov[,1] <- H$Years







