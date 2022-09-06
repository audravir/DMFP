rm(list=ls(all=TRUE))
library(readr)
library(OpenMx)
library(sfsmisc)
library(moments)
library(xtable)
library(tseries)

##------------------------
## load log returns
##------------------------

impre = read_csv("data/10_dim_daily_return.csv")
T     = dim(impre)[1]
dm    = 10 # choose the dimension
cc    = 1 # =1 open to close; =0 close to close
          # ALWAYS use =1, after standardizing white noise
rets  = data.matrix(impre[(dim(impre)[1]-T+1):(dim(impre)[1]),(cc*10+2):(cc*10+1+dm)])
date  = as.Date(impre$Date[(dim(impre)[1]-T+1):(dim(impre)[1])],format="%Y-%m-%d")
names = c( "BAC", "JPM",  "IBM",  "MSFT",
           "XOM",  "AA",   "AXP",  "DD",  
           "GE",   "KO" )


## descriptive stats for crude returns

desc = matrix(NA,ncol=10,nrow=10)

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
             label = 'table:rets_desc',
             caption='Descriptive statistics for the return data and
              $p$-values for  Shapiro-Wilk (SW), Kolmogorov-Smirnov (KS) and 
             Jarque-Bera (JB) tests for Normality as well as Ljung-Box Q-test
             for autocorrelation for lags 5 and lags 10.'), 
      file = 'tables_and_figures/rets_desc.tex')



##------------------------
## load Realized Measures
##------------------------

impre = read_csv("data/10_dim_realized_covar.csv")
reaco = data.matrix(impre[(dim(impre)[1]-T+1):(dim(impre)[1]),-1])
RCov  = RCor = array(NA,c(dm,dm,T))
stand = RVs = matrix(NA,ncol=dm,nrow=T)

for(t in 1:T){
  x         = vech2full(reaco[t,])
  RCov[,,t] = x[1:dm,1:dm]
  stand[t,] = rets[t,]/sqrt(diag(RCov[,,t]))
  RCor[,,t] = nearcor(round(cov2cor(RCov[,,t]),8))$cor
  RVs[t,]   = sqrt(diag(RCov[,,t]))
}

apply(stand,2,mean)
apply(stand,2,sd)

#pdf('tables_and_figures/rets_desc.pdf',width = 11,height = 5)
par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {plot(date,rets[,i],type='l',
    main=names[i],ylab = '',xlab = '',col='gray60')
  lines(date,sqrt(RCov[i,i,]))}
#dev.off()


par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {plot(date,stand[,i],type='l',main=paste(colnames(rets)[i]))}


par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {
  acf(stand[,i],main=paste(colnames(rets)[i]))
}

par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {
  pacf(stand[,i],main=paste(colnames(rets)[i]))
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
save(marginals,file='temp/marginals.Rdata')

stand = stand-matrix(rep(apply(stand, 2, mean),T),nrow=T,byrow = TRUE)
stand = stand/matrix(rep(apply(stand, 2, sd),T),nrow=T,byrow = TRUE)

apply(stand,2,mean)
apply(stand,2,sd)


library(QTLRel)

#pdf('tables_and_figures/standrets_qq.pdf',width = 11,height = 5)
par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {qqPlot(stand[,i],x="norm",main=names[i],
                       ylab='',xlab='',pch=20,cex=1.2,
                       plot.it=TRUE,confidence=.95)}
#dev.off()


#pdf('tables_and_figures/standrets_hist.pdf',width = 11,height = 5)
par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {hist(stand[,i],main=names[i],
  freq=FALSE,ylab='',xlab='',xlim=c(-4,4),ylim=c(0,0.5))
lines(seq(-4,4,length=500),dnorm(seq(-4,4,length=500)),
      lwd=2) }
#dev.off()


#pdf('tables_and_figures/standrets_pnorm.pdf',width = 11,height = 5)
par(mfrow=c(2,ceiling(dm/2)))
for(i in 1:dm) {hist(pnorm(stand[,i]),main=names[i],
                     ylab = '',xlab = '')}
#dev.off()


desc = matrix(NA,ncol=10,nrow=10)

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
             label = 'table:standrets_desc',
             caption='Descriptive statistics for the standardized return data and
              $p$-values for  Shapiro-Wilk (SW), Kolmogorov-Smirnov (KS) and 
             Jarque-Bera (JB) tests for Normality as well as Ljung-Box Q-test
             for autocorrelation for lags 5 and lags 10..'), 
      file = 'tables_and_figures/standrets_desc.tex')

round(desc,4)

udata  = matrix(NA,ncol=dm,nrow=dim(stand)[1])
for(i in 1:dm){
  udata[,i] = pnorm(stand[,i])
}

library(plyr)
Sigma = alply(RCor,3) #

save(Sigma,rets,RVs,stand,udata,RCov,RCor,names,date,file='10data.Rdata')




