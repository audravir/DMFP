rm(list=ls(all=TRUE))

load('data/10data.Rdata')

library(HARModel)
library(xts)
library(xtable)

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

ONEsa = matrix(NA,ncol=dm,nrow=4)
FIVEsa  = matrix(NA,ncol=dm,nrow=4)

for(i in 1:dm){
  ONEsa[1,i] = sqrt(mean((RVsforc1_m1[,i] - tail(RVs[,i],oos))^2))
  ONEsa[2,i] = sqrt(mean((RVsforc1_m2[,i] - tail(RVs[,i],oos))^2))
  ONEsa[3,i] = mean(abs(RVsforc1_m1[,i] - tail(RVs[,i],oos)))
  ONEsa[4,i] = mean(abs(RVsforc1_m2[,i] - tail(RVs[,i],oos)))

  FIVEsa[1,i] = sqrt(mean((RVsforc5_m1[,i] - tail(RVs[,i],oos))^2))
  FIVEsa[2,i] = sqrt(mean((RVsforc5_m2[,i] - tail(RVs[,i],oos))^2))
  FIVEsa[3,i] = mean(abs(RVsforc5_m1[,i] - tail(RVs[,i],oos)))
  FIVEsa[4,i] = mean(abs(RVsforc5_m2[,i] - tail(RVs[,i],oos)))
  }

colnames(ONEsa) = names
rownames(ONEsa) = c('HAR(1,5) RMSE','HAR(1,5) MAE',
                    'HAR(1,5,22) RMSE', 'HAR(1,5,22) MAE')


print(xtable(ONEsa,
             caption = "1-step-ahead realized volatility prediction results 
             for  2009/01/02-2009/12/31 out-of-sample period ($K=252$ observations).
             The table reports the root mean squared error (RMSE) and 
             mean absolute error (MAE) for HAR(1,5) and HAR (1,5,22) models of Corsi (2009).",
             label = 'table:1sa_RV', digits = 3),
      file='tables_and_figures/1sa_RV.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,scalebox=0.8)


colnames(FIVEsa) = names
rownames(FIVEsa) = c('HAR(1,5) RMSE','HAR(1,5) MAE',
                    'HAR(1,5,22) RMSE', 'HAR(1,5,22) MAE')

print(xtable(FIVEsa,
             caption = "5-step-ahead realized volatility prediction results 
             for  2009/01/06-2009/12/31 out-of-sample period ($K=248$ observations).
             The table reports the root mean squared error (RMSE) and 
             mean absolute error (MAE) for HAR(1,5) and HAR (1,5,22) models of Corsi (2009).",
             label = 'table:5sa_RV', digits = 3),
      file='tables_and_figures/5sa_RV.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,scalebox=0.8)



max(abs((ONEsa[1,]/ONEsa[2,]-1)*100))
max(abs((ONEsa[3,]/ONEsa[4,]-1)*100))


max(abs((FIVEsa[1,]/FIVEsa[2,]-1)*100))
max(abs((FIVEsa[3,]/FIVEsa[4,]-1)*100))


