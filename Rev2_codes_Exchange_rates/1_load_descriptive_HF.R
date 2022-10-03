rm(list=ls(all=TRUE))
library(readr)
library(OpenMx)
library(sfsmisc)
library(moments)
library(xtable)
library(tseries)
library(xts)

##------------------------
## load log returns
##------------------------

load('Rev2_codes_Exchange_rates/EX.Rdata')
T     = dim(rets)[1]
dm    = dim(rets)[2]
date  = as.Date(index(rets),format="%Y-%m-%d")
rets  = data.matrix((rets))
names = c("EUR/USD","USD/JPY","EUR/GBP")


## descriptive stats for crude returns

desc = matrix(NA,ncol=10,nrow=dm)

for(i in 1:dm){
  desc[i,1] = mean(rets[,i])
  desc[i,2] = median(rets[,i])
  desc[i,3] = sd(rets[,i])
  desc[i,4] = skewness(rets[,i])
  desc[i,5] = kurtosis(rets[,i])
  t = shapiro.test(rets[,i])
  desc[i,6] = t$p.value
  t = ks.test(rets[,i],y='pnorm')
  desc[i,7] = t$p.value
  t = jarque.bera.test(rets[,i])
  desc[i,8] = t$p.value
  t = Box.test (rets[,i], lag = 5, type = "Ljung")
  desc[i,9] = t$p.value
  t = Box.test (rets[,i], lag = 10, type = "Ljung")
  desc[i,10] = t$p.value
}



rownames(desc) = names
colnames(desc) = c('mean','median','sd','skew.','kurt.',
                   'SW','KS','JB','LB(5)','LB(10)')
desc

print(xtable(desc, type = "latex",digits = c(1,2,2,2,2,2,4,4,4,4,4),
             label = 'table:rets_desc_EX',
             caption='Descriptive statistics for the return data and
              $p$-values for  Shapiro-Wilk (SW), Kolmogorov-Smirnov (KS) and 
             Jarque-Bera (JB) tests for Normality as well as Ljung-Box Q-test
             for autocorrelation for lags 5 and lags 10.'), 
      file = 'tables_and_figures/rets_desc_EX.tex')



##------------------------
## load Realized Measures
##------------------------

stand =  matrix(NA,ncol=dm,nrow=T)

for(t in 1:T){
  stand[t,] = rets[t,]/sqrt(diag(RCov[,,t]))
}

apply(stand,2,mean)
apply(stand,2,sd)

pdf('tables_and_figures/rets_desc_EX.pdf',width = 10,height = 3)
par(mfrow=c(1,3))
for(i in 1:dm) {plot(date,rets[,i],type='l',
    main=names[i],ylab = '',xlab = '',col='gray60')
  lines(date,sqrt(RCov[i,i,]))}
dev.off()


par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {plot(date,stand[,i],type='l',main=paste(colnames(rets)[i]))}

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {
  acf(stand[,i],main=paste(colnames(rets)[i]),lag.max = 50)
}

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {
  pacf(stand[,i],main=paste(colnames(rets)[i]),lag.max = 50)
}

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {
  acf(stand[,i]^2,main=paste(colnames(rets)[i]))
}

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {
  pacf(stand[,i]^2,main=paste(colnames(rets)[i]))
}

stand_means = apply(stand, 2, mean)
stand_sds   = apply(stand, 2, sd)

marginals = list(stand_means,stand_sds)
names(marginals) = c('mean','sd')
save(marginals,file='temp/marginals_EX.Rdata')

stand = stand-matrix(rep(apply(stand, 2, mean),T),nrow=T,byrow = TRUE)
stand = stand/matrix(rep(apply(stand, 2, sd),T),nrow=T,byrow = TRUE)

apply(stand,2,mean)
apply(stand,2,sd)


library(QTLRel)

pdf('tables_and_figures/standrets_qq_EX.pdf',width = 10,height = 3)
par(mfrow=c(1,3))
for(i in 1:dm) {qqPlot(stand[,i],x="norm",main=names[i],
                       ylab='',xlab='',pch=20,cex=1.2,
                       plot.it=TRUE,confidence=.95,
                       ylim=c(-4,4),xlim=c(-4,4))}
dev.off()


pdf('tables_and_figures/standrets_hist_EX.pdf',width = 10,height = 3)
par(mfrow=c(1,3))
for(i in 1:dm) {hist(stand[,i],nclass = 30,main=names[i],
  freq=FALSE,ylab='',xlab='',xlim=c(-4,4),ylim=c(0,0.5))
lines(seq(-4,4,length=500),dnorm(seq(-4,4,length=500)),
      lwd=2) }
dev.off()


pdf('tables_and_figures/standrets_pnorm_EX.pdf',width = 10,height = 3)
par(mfrow=c(1,3))
for(i in 1:dm) {hist(pnorm(stand[,i]),main=names[i],
                     ylab = '',xlab = '')}
dev.off()


desc = matrix(NA,ncol=10,nrow=dm)

for(i in 1:dm){
  desc[i,1] = mean(stand[,i])
  desc[i,2] = median(stand[,i])
  desc[i,3] = sd(stand[,i])
  desc[i,4] = skewness(stand[,i])
  desc[i,5] = kurtosis(stand[,i])
  t = shapiro.test(stand[,i])
  desc[i,6] = t$p.value
  t = ks.test(stand[,i],y='pnorm')
  desc[i,7] = t$p.value
  t = jarque.bera.test(stand[,i])
  desc[i,8] = t$p.value
  t = Box.test (stand[,i], lag = 5, type = "Ljung")
  desc[i,9] = t$p.value
  t = Box.test (stand[,i], lag = 10, type = "Ljung")
  desc[i,10] = t$p.value
}

rownames(desc) = names
colnames(desc) = c('mean','median','sd','skew.','kurt.',
                   'SW','KS','JB','LB(5)','LB(10)')

print(xtable(desc, type = "latex",digits = c(1,2,2,2,2,2,4,4,4,4,4),
             label = 'table:standrets_desc_EX',
             caption='Descriptive statistics for the standardized return data and
              $p$-values for  Shapiro-Wilk (SW), Kolmogorov-Smirnov (KS) and 
             Jarque-Bera (JB) tests for Normality as well as Ljung-Box Q-test
             for autocorrelation for lags 5 and lags 10..'), 
      file = 'tables_and_figures/standrets_desc_EX.tex')

round(desc,4)

udata  = matrix(NA,ncol=dm,nrow=dim(stand)[1])
for(i in 1:dm){
  udata[,i] = pnorm(stand[,i])
}

library(plyr)
Sigma = alply(RCor,3) #

save(Sigma,rets,RVs,stand,udata,RCov,RCor,names,date,file='data/3data_EX.Rdata')




