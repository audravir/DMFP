rm(list=ls(all=TRUE))

load('data/10data.Rdata')

library(HARModel)
library(xts)

nn = dim(RVs)[1]
dm = dim(RVs)[2]

##-------------
## Compare models: 
## M1: HAR(1,5); M2: HAR(1,5,22)
##-------------

oos=252



# M1
RVsforc1_m1 = matrix(NA,ncol=ncol(RVs),nrow=oos)
RVsforc5_m1 = matrix(NA,ncol=ncol(RVs),nrow=oos)
for(i in 1:dm){
  x=as.xts(RVs[,i],order.by=date)
  ForecastHAR = HARForecast(x, periods = c(1,5), nRoll =oos,
                            nAhead = 5, type = "HAR")
  RVsforc1_m1[,i] = ForecastHAR@forecast[1,]
  RVsforc5_m1[,i] = ForecastHAR@forecast[5,]
}

# M2
RVsforc1_m2 = matrix(NA,ncol=ncol(RVs),nrow=oos)
RVsforc5_m2 = matrix(NA,ncol=ncol(RVs),nrow=oos)
for(i in 1:dm){
  x=as.xts(RVs[,i],order.by=date)
  ForecastHAR = HARForecast(x, periods = c(1,5,22), nRoll =oos,
                            nAhead = 5, type = "HAR")
  RVsforc1_m2[,i] = ForecastHAR@forecast[1,]
  RVsforc5_m2[,i] = ForecastHAR@forecast[5,]
}

RMSE = matrix(NA,ncol=dm,nrow=4)
MAE  = matrix(NA,ncol=dm,nrow=4)

for(i in 1:dm){
  RMSE[1,i] = sqrt(mean((RVsforc1_m1[,i] - tail(RVs[,i],oos))^2))
  RMSE[2,i] = sqrt(mean((RVsforc1_m2[,i] - tail(RVs[,i],oos))^2))
  RMSE[3,i] = sqrt(mean((RVsforc5_m1[,i] - tail(RVs[,i],oos))^2))
  RMSE[4,i] = sqrt(mean((RVsforc5_m2[,i] - tail(RVs[,i],oos))^2))

  MAE[1,i] = mean(abs(RVsforc1_m1[,i] - tail(RVs[,i],oos)))
  MAE[2,i] = mean(abs(RVsforc1_m2[,i] - tail(RVs[,i],oos)))
  MAE[3,i] = mean(abs(RVsforc5_m1[,i] - tail(RVs[,i],oos)))
  MAE[4,i] = mean(abs(RVsforc5_m2[,i] - tail(RVs[,i],oos)))
}


