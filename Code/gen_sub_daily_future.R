#import the necessary libraries
library(ismev)
library(DEoptim)
library(extRemes)
library(writexl)
library(dplyr)
library(readxl)


#Read GCM data#########################################
H <- read.csv("AMS_data.csv")
#Historical
Hist_5 <- data.frame(read.csv("Data/CMIP5_Hist_annual_max.csv"))
colnames(Hist_5) <- c("Year", "Rainfall")

Hist_6 <-  data.frame(read.csv("Data/CMIP6_Hist_annual_max.csv"))
colnames(Hist_6) <- c("Year", "Rainfall")

#future
Fut_5_RCP45 <- data.frame(read.csv("Data/CMIP5_RCP45_annual_max.csv"))
colnames(Fut_5_RCP45) <- c("Year", "Rainfall")

Fut_5_RCP85 <- data.frame(read.csv("Data/CMIP5_RCP85_annual_max.csv"))
colnames(Fut_5_RCP85) <- c("Year", "Rainfall")

######################################CMIP5  ###################################################################
#Historical
GEV_Hist_MODEL_CMIP5 <- fevd(Hist_5$Rainfall,Hist_5)
CDF_Hist_MODEL_CMIP5 <- exp(-(1+GEV_Hist_MODEL_CMIP5$results$par[3]*
                              (Hist_5$Rainfall-GEV_Hist_MODEL_CMIP5$ results$par[1])/
                              GEV_Hist_MODEL_CMIP5$ results$par[2])^(-1/GEV_Hist_MODEL_CMIP5$ results$par[3]))

#future
GEV_Fut_RCP45_CMIP5 <- fevd(Fut_5_RCP45$Rainfall,Fut_5_RCP45)
CDF_Fut_RCP45_CMIP5 <- exp(-(1+GEV_Fut_RCP45_CMIP5$results$par[3]*
                                (Fut_5_RCP45$Rainfall-GEV_Fut_RCP45_CMIP5$ results$par[1])/
                                GEV_Fut_RCP45_CMIP5$ results$par[2])^(-1/GEV_Fut_RCP45_CMIP5$ results$par[3]))

GEV_Fut_RCP85_CMIP5 <- fevd(Fut_5_RCP85$Rainfall,Fut_5_RCP85)
CDF_Fut_RCP85_CMIP5 <- exp(-(1+GEV_Fut_RCP85_CMIP5$results$par[3]*
                               (Fut_5_RCP85$Rainfall-GEV_Fut_RCP85_CMIP5$ results$par[1])/
                               GEV_Fut_RCP85_CMIP5$ results$par[2])^(-1/GEV_Fut_RCP85_CMIP5$ results$par[3]))


##########################################CMIP6 ###############################################################

#future
Fut_6_SSP245 <- data.frame(read.csv("Data/CMIP6_SSP245_annual_max.csv"))
colnames(Fut_6_SSP245) <- c("Year", "Rainfall")

Fut_6_SSP585 <- data.frame(read.csv("Data/CMIP6_SSP585_annual_max.csv"))
colnames(Fut_6_SSP585) <- c("Year", "Rainfall")


#Historical
GEV_Hist_MODEL_CMIP6 <- fevd(Hist_6$Rainfall,Hist_6)
CDF_Hist_MODEL_CMIP6 <- exp(-(1+GEV_Hist_MODEL_CMIP6$results$par[3]*
                                (Hist_5$Rainfall-GEV_Hist_MODEL_CMIP6$ results$par[1])/
                                GEV_Hist_MODEL_CMIP6$ results$par[2])^(-1/GEV_Hist_MODEL_CMIP6$ results$par[3]))

#future
GEV_Fut_SSP245_CMIP6 <- fevd(Fut_6_SSP245$Rainfall,Fut_6_SSP245)
CDF_Fut_SSP245_CMIP6 <- exp(-(1+GEV_Fut_SSP245_CMIP6$results$par[3]*
                               (Fut_6_SSP245$Rainfall-GEV_Fut_SSP245_CMIP6$ results$par[1])/
                                GEV_Fut_SSP245_CMIP6$ results$par[2])^(-1/GEV_Fut_SSP245_CMIP6$ results$par[3]))

GEV_Fut_SSP585_CMIP6 <- fevd(Fut_6_SSP585$Rainfall,Fut_6_SSP585)
CDF_Fut_SSP585_CMIP6 <- exp(-(1+GEV_Fut_SSP585_CMIP6$results$par[3]*
                               (Fut_6_SSP585$Rainfall-GEV_Fut_SSP585_CMIP6$ results$par[1])/
                                GEV_Fut_SSP585_CMIP6$ results$par[2])^(-1/GEV_Fut_SSP585_CMIP6$ results$par[3]))






gen_sub_daily_rcp45 <- function(col=col){
  #spatial downscaling ##############################
  GEV_H_I_1H <- fevd(col, H, units="mm")
  
  SD_MODEL_Hist_CMIP5_1H<-
    qevd(CDF_Hist_MODEL_CMIP5,GEV_H_I_1H$results$par[[1]],
         GEV_H_I_1H$results$par[[2]],GEV_H_I_1H$results$par[[3]])
  #find the adjusted parameters
  fr_MODEL_Hist_CMIP5_1H <- function(x){x1<-x[1];x2<-x[2];x3<-x[3];x4<-x[4]
  sum((SD_MODEL_Hist_CMIP5_1H-(((x1+Hist_5$Rainfall)/(x2+x3*Hist_5$Rainfall))+ (x4/Hist_5$Rainfall)))^2)}
  
  Rel1_MODEL_Hist_CMIP5_1H <- DEoptim(fr_MODEL_Hist_CMIP5_1H,lower = c(-100,-1000,-1,-100),upper = 
                                        c(1000,1000,1,1000),DEoptim.control(NP=100,itermax = 1000,trace=FALSE))
  #temporal downscaling ##############################
  TD_MODEL_RCP45 <- qevd(CDF_Fut_RCP45_CMIP5, GEV_Hist_MODEL_CMIP5$results$par[[1]],GEV_Hist_MODEL_CMIP5$results$par[[2]],GEV_Hist_MODEL_CMIP5
                         $results$par[[3]])
  
  TD_MODEL_RCP85 <- qevd(CDF_Fut_RCP85_CMIP5, GEV_Hist_MODEL_CMIP5$results$par[[1]],GEV_Hist_MODEL_CMIP5$results$par[[2]],GEV_Hist_MODEL_CMIP5
                         $results$par[[3]])
  #Scaling Factor
  Rel2_MODEL_RCP45_CMIP5_1H <- Fut_5_RCP45$Rainfall/TD_MODEL_RCP45
  
  Rel2_MODEL_RCP85_CMIP5_1H <- Fut_5_RCP85$Rainfall/TD_MODEL_RCP85
  
  #Future sub-daily data
  NEW_MODEL_RCP45_1H <-
    (((Rel1_MODEL_Hist_CMIP5_1H$optim$bestmem[[1]]+TD_MODEL_RCP45)/(Rel1_MODEL_Hist_CMIP5_1H$optim$bestmem[[2]]+
                                                                      (Rel1_MODEL_Hist_CMIP5_1H$optim$bestmem[[3]]*TD_MODEL_RCP45)))+
       ((Rel1_MODEL_Hist_CMIP5_1H$optim$bestmem[[4]]/TD_MODEL_RCP45)))*Rel2_MODEL_RCP45_CMIP5_1H
  
  return(NEW_MODEL_RCP45_1H)
}

gen_sub_daily_rcp85 <- function(col=col){
  #spatial downscaling ##############################
  GEV_H_I_1H <- fevd(col, H, units="mm")
  
  SD_MODEL_Hist_CMIP5_1H<-
    qevd(CDF_Hist_MODEL_CMIP5,GEV_H_I_1H$results$par[[1]],
         GEV_H_I_1H$results$par[[2]],GEV_H_I_1H$results$par[[3]])
  
  #find the adjusted parameters
  fr_MODEL_Hist_CMIP5_1H <- function(x){x1<-x[1];x2<-x[2];x3<-x[3];x4<-x[4]
  sum((SD_MODEL_Hist_CMIP5_1H-(((x1+Hist_5$Rainfall)/(x2+x3*Hist_5$Rainfall))+ (x4/Hist_5$Rainfall)))^2)}
  
  Rel1_MODEL_Hist_CMIP5_1H <- DEoptim(fr_MODEL_Hist_CMIP5_1H,lower = c(-100,-1000,-1,-100),upper = 
                                        c(1000,1000,1,1000),DEoptim.control(NP=100,itermax = 1000,trace=FALSE))
  #temporal downscaling ##############################
  TD_MODEL_RCP85 <- qevd(CDF_Fut_RCP85_CMIP5, GEV_Hist_MODEL_CMIP5$results$par[[1]],GEV_Hist_MODEL_CMIP5$results$par[[2]],GEV_Hist_MODEL_CMIP5
                         $results$par[[3]])
  #Scaling Factor
  Rel2_MODEL_RCP85_CMIP5_1H <- Fut_5_RCP85$Rainfall/TD_MODEL_RCP85
  
  #Future sub-daily data
  NEW_MODEL_RCP85_1H <-
    (((Rel1_MODEL_Hist_CMIP5_1H$optim$bestmem[[1]]+TD_MODEL_RCP85)/(Rel1_MODEL_Hist_CMIP5_1H$optim$bestmem[[2]]+
                                                                      (Rel1_MODEL_Hist_CMIP5_1H$optim$bestmem[[3]]*TD_MODEL_RCP85)))+
       ((Rel1_MODEL_Hist_CMIP5_1H$optim$bestmem[[4]]/TD_MODEL_RCP85)))*Rel2_MODEL_RCP85_CMIP5_1H
  
  return(NEW_MODEL_RCP85_1H)
}

gen_sub_daily_ssp245 <- function(col=col){
  #spatial downscaling ##############################
  GEV_H_I_1H <- fevd(col, H, units="mm")
  
  SD_MODEL_Hist_CMIP6_1H<-
    qevd(CDF_Hist_MODEL_CMIP6,GEV_H_I_1H$results$par[[1]],
         GEV_H_I_1H$results$par[[2]],GEV_H_I_1H$results$par[[3]])
  #find the adjusted parameters
  fr_MODEL_Hist_CMIP6_1H <- function(x){x1<-x[1];x2<-x[2];x3<-x[3];x4<-x[4]
  sum((SD_MODEL_Hist_CMIP6_1H-(((x1+Hist_6$Rainfall)/(x2+x3*Hist_6$Rainfall))+ (x4/Hist_6$Rainfall)))^2)}
  
  Rel1_MODEL_Hist_CMIP6_1H <- DEoptim(fr_MODEL_Hist_CMIP6_1H,lower = c(-100,-1000,-1,-100),upper = 
                                        c(1000,1000,1,1000),DEoptim.control(NP=100,itermax = 1000,trace=FALSE))
  #temporal downscaling ##############################
  TD_MODEL_SSP245 <- qevd(CDF_Fut_SSP245_CMIP6, GEV_Hist_MODEL_CMIP6$results$par[[1]],GEV_Hist_MODEL_CMIP6$results$par[[2]],GEV_Hist_MODEL_CMIP6
                         $results$par[[3]])
  
  TD_MODEL_SSP585 <- qevd(CDF_Fut_SSP585_CMIP6, GEV_Hist_MODEL_CMIP6$results$par[[1]],GEV_Hist_MODEL_CMIP6$results$par[[2]],GEV_Hist_MODEL_CMIP6
                         $results$par[[3]])
  #Scaling Factor
  Rel2_MODEL_SSP245_CMIP6_1H <- Fut_6_SSP245$Rainfall/TD_MODEL_SSP245
  
  Rel2_MODEL_SSP585_CMIP6_1H <- Fut_6_SSP585$Rainfall/TD_MODEL_SSP585
  
  #Future sub-daily data
  NEW_MODEL_SSP245_1H <-
    (((Rel1_MODEL_Hist_CMIP6_1H$optim$bestmem[[1]]+TD_MODEL_SSP245)/(Rel1_MODEL_Hist_CMIP6_1H$optim$bestmem[[2]]+
                                                                      (Rel1_MODEL_Hist_CMIP6_1H$optim$bestmem[[3]]*TD_MODEL_SSP245)))+
       ((Rel1_MODEL_Hist_CMIP6_1H$optim$bestmem[[4]]/TD_MODEL_SSP245)))*Rel2_MODEL_SSP245_CMIP6_1H
  
  return(NEW_MODEL_SSP245_1H)
}

gen_sub_daily_ssp585 <- function(col=col){
  #spatial downscaling ##############################
  GEV_H_I_1H <- fevd(col, H, units="mm")
  
  SD_MODEL_Hist_CMIP6_1H<-
    qevd(CDF_Hist_MODEL_CMIP6,GEV_H_I_1H$results$par[[1]],
         GEV_H_I_1H$results$par[[2]],GEV_H_I_1H$results$par[[3]])
  
  #find the adjusted parameters
  fr_MODEL_Hist_CMIP6_1H <- function(x){x1<-x[1];x2<-x[2];x3<-x[3];x4<-x[4]
  sum((SD_MODEL_Hist_CMIP6_1H-(((x1+Hist_6$Rainfall)/(x2+x3*Hist_6$Rainfall))+ (x4/Hist_6$Rainfall)))^2)}
  
  Rel1_MODEL_Hist_CMIP6_1H <- DEoptim(fr_MODEL_Hist_CMIP6_1H,lower = c(-100,-1000,-1,-100),upper = 
                                        c(1000,1000,1,1000),DEoptim.control(NP=100,itermax = 1000,trace=FALSE))
  
  #temporal downscaling ##############################
  TD_MODEL_SSP585 <- qevd(CDF_Fut_SSP585_CMIP6, GEV_Hist_MODEL_CMIP6$results$par[[1]],GEV_Hist_MODEL_CMIP6$results$par[[2]],GEV_Hist_MODEL_CMIP6
                          $results$par[[3]])
  #Scaling Factor
  Rel2_MODEL_SSP585_CMIP6_1H <- Fut_6_SSP585$Rainfall/TD_MODEL_SSP585
  
  #Future sub-daily data
  NEW_MODEL_SSP585_1H <-
    (((Rel1_MODEL_Hist_CMIP6_1H$optim$bestmem[[1]]+TD_MODEL_SSP585)/(Rel1_MODEL_Hist_CMIP6_1H$optim$bestmem[[2]]+
                                                                       (Rel1_MODEL_Hist_CMIP6_1H$optim$bestmem[[3]]*TD_MODEL_SSP585)))+
       ((Rel1_MODEL_Hist_CMIP6_1H$optim$bestmem[[4]]/TD_MODEL_SSP585)))*Rel2_MODEL_SSP585_CMIP6_1H
  
  return(NEW_MODEL_SSP585_1H)
}
H <- read.csv("AMS_data.csv")


#initialize dataframe
# Generate a sequence of numbers from 2006 to 2100
values <- seq(2006, 2100)

# Create a 95-row matrix from the sequence
years <- matrix(values, nrow = 95)
NEW_MODEL_RCP45_1H_df <- data.frame(years)
NEW_MODEL_RCP85_1H_df <- data.frame(years)

for (i in 3:50) {
  NEW_MODEL_RCP45_1H_v <- gen_sub_daily_rcp45(col= H[[i]])
  NEW_MODEL_RCP45_1H_df <- cbind(NEW_MODEL_RCP45_1H_df, NEW_MODEL_RCP45_1H_v)
}

for (i in 3:50) {
  NEW_MODEL_RCP85_1H_v <- gen_sub_daily_rcp85(col= H[[i]])
  NEW_MODEL_RCP85_1H_df <- cbind(NEW_MODEL_RCP85_1H_df, NEW_MODEL_RCP85_1H_v)
}


colnames(NEW_MODEL_RCP45_1H_df) <- c("Years", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10",
                     "H11", "H12", "H13", "H14", "H15", "H16", "H17", "H18", "H19", "H20",
                     "H21", "H22", "H23", "H24", "H25", "H26", "H27", "H28", "H29", "H30",
                     "H31", "H32", "H33", "H34", "H35", "H36", "H37", "H38", "H39", "H40",
                     "H41", "H42", "H43", "H44", "H45", "H46", "H47", "H48")
colnames(NEW_MODEL_RCP85_1H_df) <- c("Years", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10",
                                     "H11", "H12", "H13", "H14", "H15", "H16", "H17", "H18", "H19", "H20",
                                     "H21", "H22", "H23", "H24", "H25", "H26", "H27", "H28", "H29", "H30",
                                     "H31", "H32", "H33", "H34", "H35", "H36", "H37", "H38", "H39", "H40",
                                     "H41", "H42", "H43", "H44", "H45", "H46", "H47", "H48")

write.csv(NEW_MODEL_RCP45_1H_df, file = "RCP45_sub_daily_update.csv", row.names = FALSE)
write.csv(NEW_MODEL_RCP85_1H_df, file = "RCP85_sub_daily2_update.csv", row.names = FALSE)



#initialize dataframe
# Generate a sequence of numbers from 2006 to 2100
values1 <- seq(2015, 2100)

# Create a 95-row matrix from the sequence
years <- matrix(values1, nrow = 86)
NEW_MODEL_SSP245_1H_df <- data.frame(years)
NEW_MODEL_SSP585_1H_df <- data.frame(years)

for (i in 3:50) {
  NEW_MODEL_SSP245_1H_v <- gen_sub_daily_ssp245(col= H[[i]])
  NEW_MODEL_SSP245_1H_df <- cbind(NEW_MODEL_SSP245_1H_df, NEW_MODEL_SSP245_1H_v)
}

for (i in 3:50) {
  NEW_MODEL_SSP585_1H_v <- gen_sub_daily_ssp585(col= H[[i]])
  NEW_MODEL_SSP585_1H_df <- cbind(NEW_MODEL_SSP585_1H_df, NEW_MODEL_SSP585_1H_v)
}

colnames(NEW_MODEL_SSP245_1H_df) <- c("Years", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10",
                                     "H11", "H12", "H13", "H14", "H15", "H16", "H17", "H18", "H19", "H20",
                                     "H21", "H22", "H23", "H24", "H25", "H26", "H27", "H28", "H29", "H30",
                                     "H31", "H32", "H33", "H34", "H35", "H36", "H37", "H38", "H39", "H40",
                                     "H41", "H42", "H43", "H44", "H45", "H46", "H47", "H48")

colnames(NEW_MODEL_SSP585_1H_df) <- c("Years", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10",
                                     "H11", "H12", "H13", "H14", "H15", "H16", "H17", "H18", "H19", "H20",
                                     "H21", "H22", "H23", "H24", "H25", "H26", "H27", "H28", "H29", "H30",
                                     "H31", "H32", "H33", "H34", "H35", "H36", "H37", "H38", "H39", "H40",
                                     "H41", "H42", "H43", "H44", "H45", "H46", "H47", "H48")

write.csv(NEW_MODEL_SSP245_1H_df, file = "SSP245_sub_daily_update.csv", row.names = FALSE)
write.csv(NEW_MODEL_SSP585_1H_df, file = "SSP585_sub_daily_update.csv", row.names = FALSE)








