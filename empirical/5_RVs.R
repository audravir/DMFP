rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')

library(HARModel)
library(xts)
library(xtable)
library(car)
library(MCS)
library(rugarch)

nn = dim(RVs)[1]
dm = dim(RVs)[2]

##-------------
## Compare models: 
##-------------

oos=797

# M1-HAR(1,5)
RVsforc1_m1 = matrix(NA,ncol=ncol(RVs),nrow=oos)
for(i in 1:dm){
  x=as.xts(RVs[,i]^2,order.by=date)
  ForecastHAR = HARForecast(x, periods = c(1,5), nRoll =oos,
                            nAhead = 1, type = "HAR")
  RVsforc1_m1[,i] = sqrt(ForecastHAR@forecast[1,])
}

# M2-full HAR(1,5,22)
RVsforc1_m2 = matrix(NA,ncol=ncol(RVs),nrow=oos)
for(i in 1:dm){
  x=as.xts(RVs[,i]^2,order.by=date)
  ForecastHAR = HARForecast(x, periods = c(1,5,22), nRoll =oos,
                            nAhead = 1, type = "HAR")
  RVsforc1_m2[,i] = sqrt(ForecastHAR@forecast[1,])
}


# M3-log-HAR(1,5,22)
RVsforc1_m3 = matrix(NA,ncol=dm,nrow=oos)
for(i in 1:dm){
  rvt   = log(RVs[,i])
  rvt1  = c(NA,rvt[-length(rvt)])
  rvt5  = rollmean(rvt1,5,fill = c(NA,NA,NA),align = 'right')
  rvt22 = rollmean(rvt1,22,fill = c(NA,NA,NA),align = 'right')
  dtf   = data.frame(rvt,rvt1,rvt5,rvt22)
  lrvm4 = lm(rvt~rvt1+rvt5+rvt22,data = dtf[1:(nn-oos),])
  RVsforc1_m3[,i]  = exp(predict(lrvm4, newdata = dtf[(nn-oos+1):nn,]))
}

par(mfrow=c(1,1))
plot(tail(RVs[,i],oos),type='l',lwd=2)
lines(RVsforc1_m1[,i],col=2)
lines(RVsforc1_m2[,i],col=3)
lines(RVsforc1_m3[,i],col=4)

##-------
## Calculate all kinds of loss functions 
##-------


# # https://www.rdocumentation.org/packages/MCS/versions/0.1.3/topics/LossVol
# 

lossf= c('SE1','SE2','QLIKE','R2LOG','AE1','AE2')

losses=matrix(NA,ncol=dm,nrow=3*length(lossf))
rel.losses=matrix(NA,ncol=dm,nrow=3*length(lossf))

for(i in 1:dm){
  for(j in 1:length(lossf)){
    losses[(j*3-2),i]=mean(LossVol(tail(RVs[,i],oos), RVsforc1_m1[,i], which = lossf[j]))
    losses[(j*3-1),i]= mean(LossVol(tail(RVs[,i],oos), RVsforc1_m2[,i], which = lossf[j]))
    losses[(j*3),i]=mean(LossVol(tail(RVs[,i],oos), RVsforc1_m3[,i], which = lossf[j]))
    m = min(losses[(j*3-2):(j*3),i])
    rel.losses[(j*3-2),i]=losses[(j*3-2),i]/m
    rel.losses[(j*3-1),i]=losses[(j*3-1),i]/m
    rel.losses[(j*3),i]=losses[(j*3),i]/m
  }
}

round(losses,2)
round(rel.losses,4)


df=data.frame(c(NA,'MSE1',NA,NA,'MSE2',NA,NA,'QLIKE',NA,NA,'R2LOG',NA,NA,'MAE1',NA,NA,'MAE2',NA)
              ,rep(c('HAR(1,5)','HAR(1,5,22)','log-HAR(1,5,22)'),6),
              rel.losses)


names(df) = c('','',assets)

df



print(xtable(df,
             caption = "1-step-ahead realized volatility prediction results
             for  2020/01/02-2023/0/31 out-of-sample period ($K=797$ observations).
             The table reports the relative values (with respect to the lowest value)
             of the six loss functions for  HAR(1,5), HAR (1,5,22) 
             and log-HAR(1,5,22) models.",
             label = 'table:1sa_RV', digits = 3),
      file='tables_and_figures/1sa_RV.tex',
      include.rownames = FALSE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,scalebox=0.5,
      hline.after = getOption("xtable.hline.after", c(-1,0,3,6,9,12,15,nrow(df))))


##-----------------
## Mincer Zarnowitz regressions to evaluate prediction accuracy
##-----------------

MZR2 =MZpv =matrix(NA,nrow=3,ncol = dm)

for(i in 1:dm){
  
  Y = tail(RVs[,i],oos)
  exog = (RVsforc1_m1[,i])
  reg.MZ = lm(Y~exog)
  MZR2[1,i] = summary(reg.MZ)[[8]]
  htest   = linearHypothesis(reg.MZ, c("(Intercept) = 0", "exog = 1"))   
  MZpv[1,i] = htest[[6]][2]
  
  Y = tail(RVs[,i],oos)
  exog = (RVsforc1_m2[,i])
  reg.MZ = lm(Y~exog)
  MZR2[2,i] = summary(reg.MZ)[[8]]
  htest   = linearHypothesis(reg.MZ, c("(Intercept) = 0", "exog = 1"))   
  MZpv[2,i] = htest[[6]][2]
  
  Y = log(tail(RVs[,i],oos))
  exog = log(RVsforc1_m3[,i])
  reg.MZ = lm(Y~exog)
  MZR2[3,i] = summary(reg.MZ)[[8]]
  htest   = linearHypothesis(reg.MZ, c("(Intercept) = 0", "exog = 1"))
  MZpv[3,i] = htest[[6]][2]
}

round(MZR2,2)
round(MZpv,3)


plot(log(tail(RVs[,i],oos)),type='l',lwd=2)
lines(log(RVsforc1_m1[,i]),col=2,lwd=2)
lines(log(RVsforc1_m2[,i]),col=3,lwd=2)
lines(log(RVsforc1_m3[,i]),col=4,lwd=2)



colnames(MZR2) = assets
colnames(MZpv) = assets
rownames(MZR2) = c('HAR(1,5)','HAR(1,5,22)','log-HAR(1,5,22)')
rownames(MZpv) = c('HAR(1,5)','HAR(1,5,22)','log-HAR(1,5,22)')

print(xtable(MZR2*100,
             caption = "The  Mincer Zarnowitz regression results.
             The table reports the $R^2$ for a regression of
              the true observed realized volatilities on 1-step-ahead forecasts from 
             HAR(1,5), HAR(1,5,22) and log-HAR(1,5,22) models.",
             label = 'table:MZR2', digits = 1),
      file='tables_and_figures/MZR2.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,scalebox=0.6)


print(xtable(MZpv,
             caption = "The  Mincer Zarnowitz regression results.
             The table reports the $p-$values for the joint  hypothesis
             that the intercept is zero and the slope is one
             in a regression of
              the true observed realized volatilities on 1-step-ahead 
             forecasts from HAR(1,5), HAR(1,5,22) and log-HAR(1,5,22) models.",
             label = 'table:MZpv', digits = 3),
      file='tables_and_figures/MZpv.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,scalebox=0.6)




##-----------------
## test for forecast equality via Diebold-Mariano type test
##-----------------

library(forecast)

DMtests = matrix(NA,ncol=dm,nrow=2)

em1ONEsa=em2ONEsa=matrix(NA,ncol=dm,nrow=oos)

for(i in 1:dm){
  em1ONEsa[,i] = tail(RVs[,i],oos)-RVsforc1_m2[,i]
  em2ONEsa[,i] = tail(RVs[,i],oos)-RVsforc1_m3[,i]
  
  DMtests[1,i] = dm.test(em1ONEsa[,i],em2ONEsa[,i],h = 1,power = 1,
                         alternative = "greater")$p.value
  DMtests[2,i] = dm.test(em1ONEsa[,i],em2ONEsa[,i],h = 1,power = 2,
                         alternative = "greater")$p.value
}

colnames(DMtests) = assets
rownames(DMtests) = c('1-s-a power=1','1-s-a power=2')
round(DMtests,3)
# 

# 
# print(xtable(DMtests,
#              caption = "The equality of prediction accuracy test results.
#              The table reports the $p$-values for the Diebold-Mariano 
#              test for 1-step-ahead  realized volatility predictions from HAR(1,5)
#              and HAR (1,5,22) models. The null hypothesis is equal forecast 
#              accuracy for two forecasts,
#              the alternative is two-sided. The  power controls the
#              power used in the loss function (either 1 or 2).",
#              label = 'table:DM_RV', digits = 3),
#       file='tables_and_figures/DM_RV.tex',
#       include.rownames = TRUE,latex.environments = "center" ,
#       caption.placement = "top",
#       include.colnames= TRUE,
#       rotate.colnames = FALSE,scalebox=0.8)


##-------------
## select predictions produced by the best model
##-------------


# RV_forc = list(tail(date,oos),
#                RVsforc1_m1)

RV_forc = list(tail(date,oos),RVsforc1_m3)
names(RV_forc) = c('Date','1sa')
save(RV_forc,file='empirical/temp/RV_forc.Rdata')

date_ticks <- seq(from=min(tail(date,oos)), to=max(tail(date,oos)), by="month")

pdf(file='tables_and_figures/har_rvs.pdf',width=20,height=8)
par(mfrow=c(3,5), mar=c(5, 3, 3, 3))
for(i in 1:dm) {
  plot(tail(date,oos),tail(RVs[,i],oos),type='l', xaxt='n',col='gray80',
       main=assets[i],lwd=4,ylab='',xlab='')
  axis.Date(1, at=date_ticks,  las=2,cex.axis=1, format="%Y/%m") # Rotate labels if needed
  lines(tail(date,oos),RV_forc$`1sa`[,i],lwd=2)
}
dev.off()





