library(HyetosMinute)
### ** Examples

# Example 1: Disaggregate daily rainfall into hourly values
# with the random parameter Bartlett-Lewis gamma model
# Import daily time series from "HistDailyData.txt" file

# To load the data set "HistDailyData" use
data(HistDailyData)

# To export the daily rainfall depths of "HistDailyData"
# data set in the chosen working directory use
write.table(HistDailyData,file="HistDailyData.txt",sep="\t",
            quote=FALSE,row.names=FALSE,col.names=FALSE)

# To disaggregate the daily rainfall depths of first 5 clusters of wet days use
ex21 <- DisagSimul(TimeScale=1,BLpar=list(lambda=0.569748,phi=0.048387,kappa=0.5996395,
                                          alpha=7.2933199,v=0.052517913,mx=30.4825,sxmx=32.391/30.4825),CellIntensityProp=list(Weibull=FALSE,
                                                                                                                               iota=NA),RepetOpt=list(DistAllowed=0.1,FacLevel1Rep=20,MinLevel1Rep=50,
                                                                                                                                                      TotalRepAllowed=5000),NumOfSequences=5,Statistics=list(print=TRUE,plot=FALSE),
                   ExportSynthData=list(exp=TRUE,FileContent=c("AllDays"),file="HistDailyData_Disag.txt"),
                   ImportHistData=list(file="HistDailyData.txt",na.values="NA",FileContent=c("AllDays"),
                                       DaysPerSeason=30),PlotHyetographs=FALSE,RandSeed=5)
 
# Example 2: Disaggregate daily rainfall into 10-min values
# with the random parameter Bartlett-Lewis model
# Import daily time series from "HistDailyData2.txt" file

# To load the data set "HistDailyData2" use
data(HistDailyData2)

# To export the daily rainfall depths of "HistDailyData2" data set
# in the chosen working directory use
write.table(HistDailyData2,file="HistDailyData2.txt",sep="\t",
            quote=FALSE,row.names=FALSE,col.names=FALSE)

# To disaggregate the daily rainfall depths of first 5 clusters of wet days use
x22 <- DisagSimul(TimeScale=1/6,BLpar=list(lambda=0.9396,phi=0.0568,kappa=1.05819,
                                           alpha=2.69519,v=0.006282916666,mx=24.33408,sxmx=1),CellIntensityProp=list(Weibull=FALSE,
                                                                                                                     iota=NA),RepetOpt=list(DistAllowed=0.1,FacLevel1Rep=20,MinLevel1Rep=50,
                                                                                                                                            TotalRepAllowed=5000),NumOfSequences=5,Statistics=list(print=TRUE,plot=FALSE),
                  ExportSynthData=list(exp=TRUE,FileContent=c("WetDays"),file="HistDailyData2_Disag.txt"),
                  ImportHistData=list(file="HistDailyData2.txt",na.values="NA",
                                      FileContent=c("WetDays"),DaysPerSeason=31),
                  PlotHyetographs=FALSE,RandSeed=5)
