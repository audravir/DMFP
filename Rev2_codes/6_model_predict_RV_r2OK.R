rm(list=ls(all=TRUE))

load('data/10data.Rdata')

library(HARModel)
library(xts)
library(xtable)
library(car)

nn = dim(RVs)[1]
dm = dim(RVs)[2]

##-------------
## Compare models: 
## M1: HAR(1,5); M2: HAR(1,5,22)
##-------------

oos=252

# M1
RVsforc1_m1 = matrix(NA,ncol=ncol(RVs),nrow=oos)
#RVsforc5_m1 = matrix(NA,ncol=ncol(RVs),nrow=oos)
for(i in 1:dm){
  x=as.xts(RVs[,i],order.by=date)
  ForecastHAR = HARForecast(x, periods = c(1,5), nRoll =oos,
                            nAhead = 5, type = "HAR")
  RVsforc1_m1[,i] = ForecastHAR@forecast[1,]
  #RVsforc5_m1[,i] = ForecastHAR@forecast[5,]
}

# M2
RVsforc1_m2 = matrix(NA,ncol=ncol(RVs),nrow=oos)
#RVsforc5_m2 = matrix(NA,ncol=ncol(RVs),nrow=oos)
for(i in 1:dm){
  x=as.xts(RVs[,i],order.by=date)
  ForecastHAR = HARForecast(x, periods = c(1,5,22), nRoll =oos,
                            nAhead = 5, type = "HAR")
  RVsforc1_m2[,i] = ForecastHAR@forecast[1,]
  #RVsforc5_m2[,i] = ForecastHAR@forecast[5,]
}

# ##-------------
# ## log(RV) models M3 and M4; same or worse
# ##-------------
# 
# RVsforc1_m3 = RVsforc1_m4 = matrix(NA,ncol=dm,nrow=oos)
# 
# for(i in 1:dm){
#   rvt   = log(RVs[,i]) 
#   rvt1  = c(NA,rvt[-length(rvt)])  
#   rvt5  = rollmean(rvt1,5,fill = c(NA,NA,NA),align = 'right')
#   rvt22 = rollmean(rvt1,22,fill = c(NA,NA,NA),align = 'right')
#   dtf   = data.frame(rvt,rvt1,rvt5,rvt22)
#   lrvm3 = lm(rvt~rvt1,data = dtf[1:(nn-oos),])
#   RVsforc1_m3[,i]  = exp(predict(lrvm3, newdata = dtf[(nn-oos+1):nn,]))
#   lrvm4 = lm(rvt~rvt1+rvt5+rvt22,data = dtf[1:(nn-oos),])
#   RVsforc1_m4[,i]  = exp(predict(lrvm4, newdata = dtf[(nn-oos+1):nn,]))
# }

##

ONEsa_mae  = matrix(NA,ncol=dm,nrow=2)
ONEsa_rmse = matrix(NA,ncol=dm,nrow=2)
#FIVEsa  = matrix(NA,ncol=dm,nrow=4)

em1ONEsa = em2ONEsa = em3ONEsa = em4ONEsa = matrix(NA,nrow=oos,ncol=dm)
#em1FIVEsa = em2FIVEsa = matrix(NA,nrow=oos,ncol=dm)


for(i in 1:dm){
  em1ONEsa[,i] = RVsforc1_m1[,i] - tail(RVs[,i],oos)
  em2ONEsa[,i] = RVsforc1_m2[,i] - tail(RVs[,i],oos)
  # em3ONEsa[,i] = RVsforc1_m3[,i] - tail(RVs[,i],oos)
  # em4ONEsa[,i] = RVsforc1_m4[,i] - tail(RVs[,i],oos)
  
  
  ONEsa_rmse[1,i] = sqrt(mean((em1ONEsa[,i])^2))
  ONEsa_rmse[2,i] = sqrt(mean((em2ONEsa[,i])^2))
  # ONEsa_rmse[3,i] = sqrt(mean((em3ONEsa[,i])^2))
  # ONEsa_rmse[4,i] = sqrt(mean((em4ONEsa[,i])^2))
  
  ONEsa_mae[1,i] = mean(abs(em1ONEsa[,i]))
  ONEsa_mae[2,i] = mean(abs(em2ONEsa[,i]))
  # ONEsa_mae[3,i] = mean(abs(em3ONEsa[,i]))
  # ONEsa_mae[4,i] = mean(abs(em4ONEsa[,i]))

  # em1FIVEsa[,i] = RVsforc5_m1[,i] - tail(RVs[,i],oos)
  # em2FIVEsa[,i] = RVsforc5_m2[,i] - tail(RVs[,i],oos)
  # 
  # FIVEsa[1,i] = sqrt(mean((em1FIVEsa[,i])^2))
  # FIVEsa[2,i] = sqrt(mean((em2FIVEsa[,i])^2))
  # FIVEsa[3,i] = mean(abs(em1FIVEsa[,i]))
  # FIVEsa[4,i] = mean(abs(em2FIVEsa[,i]))
  }

ONEsa = rbind(ONEsa_rmse,ONEsa_mae)

colnames(ONEsa) = names
rownames(ONEsa) = c('HAR(1,5) RMSE','HAR(1,5,22) RMSE',
                    'HAR(1,5) MAE','HAR(1,5,22) MAE')


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


# colnames(FIVEsa) = names
# rownames(FIVEsa) = c('HAR(1,5) RMSE','HAR(1,5,22) RMSE',
#                     'HAR(1,5) MAE','HAR(1,5,22) MAE')
# 
# print(xtable(FIVEsa,
#              caption = "5-step-ahead realized volatility prediction results 
#              for  2009/01/06-2009/12/31 out-of-sample period ($K=248$ observations).
#              The table reports the root mean squared error (RMSE) and 
#              mean absolute error (MAE) for HAR(1,5) and HAR (1,5,22) models of Corsi (2009).",
#              label = 'table:5sa_RV', digits = 3),
#       file='tables_and_figures/5sa_RV.tex',
#       include.rownames = TRUE,latex.environments = "center" ,
#       caption.placement = "top",
#       include.colnames= TRUE,
#       rotate.colnames = FALSE,scalebox=0.8)



max(abs((ONEsa[1,]/ONEsa[2,]-1)*100))
max(abs((ONEsa[3,]/ONEsa[4,]-1)*100))


# max(abs((FIVEsa[1,]/FIVEsa[2,]-1)*100))
# max(abs((FIVEsa[3,]/FIVEsa[4,]-1)*100))

##-----------------
## test for forecast equality via Diebold-Mariano type test
##-----------------

library(forecast)

# DMtests = matrix(NA,ncol=dm,nrow=4)
DMtests = matrix(NA,ncol=dm,nrow=2)

for(i in 1:dm){
  DMtests[1,i] = dm.test(em1ONEsa[,i],em2ONEsa[,i],h = 1,power = 1,
                         alternative = "two.sided")$p.value
  DMtests[2,i] = dm.test(em1ONEsa[,i],em2ONEsa[,i],h = 1,power = 2,
                         alternative = "two.sided")$p.value
  # DMtests[3,i] = dm.test(em1FIVEsa[,i],em2FIVEsa[,i],h = 5,power = 1,
  #                        alternative = "two.sided")$p.value
  # DMtests[4,i] = dm.test(em1FIVEsa[,i],em2FIVEsa[,i],h = 5,power = 2,
  #                        alternative = "two.sided")$p.value
}

colnames(DMtests) = names
rownames(DMtests) = c('1-s-a power=1','1-s-a power=2')


round(DMtests,3)

print(xtable(DMtests,
             caption = "The equality of prediction accuraccy test results.
             The table reports the $p$-values for the Diebold-Mariano 
             test for 1-step-ahead (1-s-a) realized volatility predictions from HAR(1,5)
             and HAR (1,5,22) models. The null hypothesis is equal forecast 
             accuracy for two forecasts,
             the alternative is two-sided. The  power controls the
             power used in the loss function (either 1 or 2).",
             label = 'table:DM_RV', digits = 3),
      file='tables_and_figures/DM_RV.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,scalebox=0.8)

##-----------------
## Mincer Zarnowitz regressions to evaluate prediction accuracy
##-----------------

# MZR2 =MZpv =matrix(NA,nrow=4,ncol = dm)
MZR2 =MZpv =matrix(NA,nrow=2,ncol = dm)

for(i in 1:dm){
  exog = RVsforc1_m1[,i]
  reg.MZ = lm(tail(RVs[,i],oos)~exog)
  MZR2[1,i] = summary(reg.MZ)[[8]]
  htest   = linearHypothesis(reg.MZ, c("(Intercept) = 0", "exog = 1"))   
  MZpv[1,i] = htest[[6]][2]
  
  exog = RVsforc1_m2[,i]
  reg.MZ = lm(tail(RVs[,i],oos)~exog)
  MZR2[2,i] = summary(reg.MZ)[[8]]
  htest   = linearHypothesis(reg.MZ, c("(Intercept) = 0", "exog = 1"))   
  MZpv[2,i] = htest[[6]][2]
  
  # exog = RVsforc1_m3[,i]
  # reg.MZ = lm(tail(RVs[,i],oos)~exog)
  # MZR2[3,i] = summary(reg.MZ)[[8]]
  # htest   = linearHypothesis(reg.MZ, c("(Intercept) = 0", "exog = 1"))   
  # MZpv[3,i] = htest[[6]][2]
  # 
  # 
  # exog = RVsforc1_m4[,i]
  # reg.MZ = lm(tail(RVs[,i],oos)~exog)
  # MZR2[4,i] = summary(reg.MZ)[[8]]
  # htest   = linearHypothesis(reg.MZ, c("(Intercept) = 0", "exog = 1"))   
  # MZpv[4,i] = htest[[6]][2]
  
  # exog = RVsforc5_m1[,i]
  # reg.MZ = lm(tail(RVs[,i],oos)~exog)
  # MZR2[3,i] = summary(reg.MZ)[[8]]
  # htest   = linearHypothesis(reg.MZ, c("(Intercept) = 0", "exog = 1"))   
  # MZpv[3,i] = htest[[6]][2]
  # 
  # exog = RVsforc5_m2[,i]
  # reg.MZ = lm(tail(RVs[,i],oos)~exog)
  # MZR2[4,i] = summary(reg.MZ)[[8]]
  # htest   = linearHypothesis(reg.MZ, c("(Intercept) = 0", "exog = 1"))   
  # MZpv[4,i] = htest[[6]][2]
}

colnames(MZR2) = names
colnames(MZpv) = names
rownames(MZR2) = c('1-s-a HAR(1,5)','1-s-a HAR(1,5,22)')
rownames(MZpv) = c('1-s-a HAR(1,5)','1-s-a HAR(1,5,22)')



print(xtable(MZR2*100,
             caption = "The  Mincer Zarnowitz regression results.
             The table reports the $R^2$ for a regression of 
              the true observed realized volatilities on 1-step-ahead (1-s-a)  forecasts from HAR(1,5) and HAR(1,5,22) models.",
             label = 'table:MZR2', digits = 1),
      file='tables_and_figures/MZR2.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,scalebox=0.8)


print(xtable(MZpv,
             caption = "The  Mincer Zarnowitz regression results.
             The table reports the $p-$values for the joint  hypothesis 
             that the intercept is zero and the slope is one
             in a regression of 
              the true observed realized volatilities on 1-step-ahead (1-s-a) 
             forecasts from HAR(1,5) and HAR(1,5,22) models.",
             label = 'table:MZpv', digits = 3),
      file='tables_and_figures/MZpv.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,scalebox=0.8)



##-------------
## select predictions produced by HAR(1,5)
##-------------


RV_forc = list(tail(date,oos),
               RVsforc1_m1)

# RV_forc = list(tail(date,oos),
#                RVsforc1_m1,
#                RVsforc5_m1)

# names(RV_forc) = c('Date','1sa','5sa')
names(RV_forc) = c('Date','1sa')


save(RV_forc,file='temp/RV_forc.Rdata')



