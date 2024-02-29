rm(list=ls(all=TRUE))

library(moments)
library(xtable)
load('data/FXdata.Rdata')

nn = dim(rets)[1]
dm = dim(rets)[2]

assets = assets[-c(10,15,16)]

date_ticks <- seq(from=min(date), to=max(date), by="year")
year_ticks <- format(date_ticks, "%Y") # Extracts the full year
last_two_digits <- substr(year_ticks, 3, 4)

##----------------------------------------
## descriptive stats for crude returns
##----------------------------------------


lbM     = function(x,M){t= Box.test(x,lag=M,type = "Ljung-Box");pv=t$p.value;return(pv)}
ks.pv   = function(x) {t = ks.test(x,y='pnorm');pv=t$p.value;return(pv)}
jb.pv   = function(x) {t = tseries::jarque.bera.test(x);pv=t$p.value;return(pv)}
ddst.pv = function(x) {t = ddst::ddst.uniform.test(x,compute.p = TRUE);pv=t$p.value;return(pv)}

desc = cbind(apply(rets,2,mean),
             apply(rets,2,quantile,0.5),
             apply(rets,2,sd),
             apply(rets,2,skewness),
             apply(rets,2,kurtosis),
             apply(rets,2,ks.pv),
             apply(rets,2,jb.pv),
             apply(rets,2,lbM,10),
             apply(rets^2,2,lbM,10))

ncol(desc)
nrow(desc)

rownames(desc) = assets
colnames(desc) = c('mean','median','sd','skew.','kurt.',
                   'KS','JB','LB(10)','ARCH(10)')
desc

print(xtable(desc, type = "latex",digits = c(1,2,2,2,2,2,4,4,4,4),
             label = 'table:rets_desc_FX',
             caption='Descriptive statistics for the return data and
              $p$-values for Kolmogorov-Smirnov (KS), 
             and Jarque-Bera (JB) tests for Normality, 
             Ljung-Box Q-test for autocorrelation and ARCH test for lag 10.'), 
      file = 'tables_and_figures/rets_desc_FX.tex')


pdf('tables_and_figures/rets_desc_FX.pdf',width = 11,height = 7)
par(mfrow=c(3,ceiling(dm/3)))
for(i in 1:dm) {
  plot(date,rets[,i],type='l',main=assets[i],ylab = '',xlab = '',col='gray60', xaxt='n')
  lines(date,sqrt(RCov[i,i,]))
  axis(1, at=date_ticks, labels=paste0("'", last_two_digits), las=2,cex.axis=1) # Rotate labels if needed
}
dev.off()

##----------------------------------------
## descriptive stats for standardized returns
##----------------------------------------


library(QTLRel)

pdf('tables_and_figures/standrets_qq_FX.pdf',width = 11,height = 7)
par(mfrow=c(3,ceiling(dm/3)))
for(i in 1:dm) {qqPlot(stand[,i],x="norm",main=assets[i],
                       ylab='',xlab='',pch=20,cex=1.2,
                       plot.it=TRUE,confidence=.95)}
dev.off()


pdf('tables_and_figures/standrets_hist_FX.pdf',width = 11,height = 7)
par(mfrow=c(3,ceiling(dm/3)))
for(i in 1:dm) {hist(stand[,i],main=assets[i],
                     freq=FALSE,ylab='',xlab='',xlim=c(-4,4),ylim=c(0,0.5))
  lines(seq(-4,4,length=500),dnorm(seq(-4,4,length=500)),
        lwd=2) }
dev.off()


pdf('tables_and_figures/standrets_pnorm_FX.pdf',width = 11,height = 7)
par(mfrow=c(3,ceiling(dm/3)))
for(i in 1:dm) {hist(pnorm(stand[,i]),main=assets[i],
                     ylab = '',xlab = '',freq=FALSE)
  abline(h=1,lwd=2)}
dev.off()


desc = cbind(apply(stand,2,mean),
             apply(stand,2,quantile,0.5),
             apply(stand,2,sd),
             apply(stand,2,skewness),
             apply(stand,2,kurtosis),
             apply(stand,2,ks.pv),
             apply(stand,2,jb.pv),
             apply(pnorm(stand),2,ddst.pv),
             apply(stand,2,lbM,10),
             apply(stand^2,2,lbM,10))

ncol(desc)

rownames(desc) = assets
colnames(desc) = c('mean','median','sd','skew.','kurt.',
                   'KS','JB','DD','LB(10)','ARCH(10)')
desc
round(desc[,10],3)

print(xtable(desc, type = "latex",digits = c(1,2,2,2,2,2,4,4,4,4,4),
             label = 'table:standrets_desc_FX',
             caption='Descriptive statistics for the standardized return data and
              $p$-values for Kolmogorov-Smirnov (KS) and 
             Jarque-Bera (JB) tests for Normality,
             data-driven smooth (DD) test for Uniformity, 
             Ljung-Box Q-test for autocorrelation and ARCH test for lag 10.'), 
      file = 'tables_and_figures/standrets_desc_FX.tex')

round(desc,4)

##----------------------------------------
## descriptive plots for EUR/USD and EUR/GBP
##----------------------------------------

pdf('tables_and_figures/desc_EURUSDGBP.pdf',width = 14,height = 6)
par(mfrow=c(2,4))

plot(date,rets[,1],type='l',main='EUR/USD: rt and RVt',ylab = '',xlab = '',col='gray60', xaxt='n')
lines(date,sqrt(RCov[1,1,]))
axis(1, at=date_ticks, labels=paste0("'", last_two_digits), las=2,cex.axis=1) # Rotate labels if needed
qqPlot(stand[,1],x="norm",main='EUR/USD: QQ-plot',ylab='',xlab='',pch=20,cex=1.2,plot.it=TRUE,confidence=.95)
hist(stand[,1],main='EUR/USD: rt and N(0,1)',freq=FALSE,ylab='',xlab='',xlim=c(-4,4),ylim=c(0,0.5))
lines(seq(-4,4,length=500),dnorm(seq(-4,4,length=500)), lwd=2)
hist(pnorm(stand[,1]),main='EUR/USD: ut', ylab = '',xlab = '',freq=FALSE)
abline(h=1,lwd=2)

plot(date,rets[,2],type='l',main='EUR/GBP: rt and RVt',ylab = '',xlab = '',col='gray60', xaxt='n')
lines(date,sqrt(RCov[2,2,]))
axis(1, at=date_ticks, labels=paste0("'", last_two_digits), las=2,cex.axis=1) # Rotate labels if needed
qqPlot(stand[,2],x="norm",main='EUR/GBP: QQ-plot',ylab='',xlab='',pch=20,cex=1.2,plot.it=TRUE,confidence=.95)
hist(stand[,2],main='EUR/GBP: rt and N(0,1)',freq=FALSE,ylab='',xlab='',xlim=c(-4,4),ylim=c(0,0.5))
lines(seq(-4,4,length=500),dnorm(seq(-4,4,length=500)), lwd=2)
hist(pnorm(stand[,2]),main='EUR/GBP: ut', ylab = '',xlab = '',freq=FALSE)
abline(h=1,lwd=2)
dev.off()

