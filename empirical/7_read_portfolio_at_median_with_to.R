rm(list=ls(all=TRUE))
library(xtable)

load('empirical/temp/FX_portfolio_at_median_0619_15k.Rdata')

# compare weights for model 1
par(mfrow=c(3,1))

plot(ws_gmvr[1,,1],type='l',ylim=c(-1,1),main='GMVr')
for(i in 1:dm) lines(ws_gmvr[1,,i])

plot(ws_CVAR10r[1,,1],type='l',ylim=c(-1,1),main='CVAR10r')
for(i in 1:dm) lines(ws_CVAR10r[1,,i])

plot(ws_CVAR05r[1,,1],type='l',ylim=c(-1,1),main='CVAR05r')
for(i in 1:dm) lines(ws_CVAR05r[1,,i])

# how many assets are involved in the restricted portfolio
par(mfrow=c(3,2))
for(i in 1:length(models)){
  plot(apply((ws_gmvr[i,,]>0.01),1,sum),type='l',ylim=c(4,8),
       main=models[i])
}

############################################################

esfun=function(x,p){
  es=mean(x[which(x<quantile(x,p))])
  return(es)
}

load('data/rf.RData')
load('data/FXdata.Rdata')

## var-covar matrices

all.cov = list(cov1,cov2,cov3,cov4,cov5,cov6)
all.ws  = list(ws_gmvr,ws_CVAR05r,ws_CVAR10r)

names(all.ws) = c('gmvr','CVAR05r','CVAR10r')

all.prets = all.co = all.to = all.portsd = all.model.sd = list()
pret      = co = to = portsd = model.sd =array(NA,dim=c(length(models),K))

for(j in 1:length(all.ws)){
  for(i in 1:length(models)){
    for(t in 1:K){
      pret[i,t]   = sum(all.ws[[j]][i,t,]*rets[nn+t,])
      portsd[i,t] = sqrt(t(all.ws[[j]][i,t,])%*%RCov[,,nn+t]%*%all.ws[[j]][i,t,])
      co[i,t]     = (sum((all.ws[[j]][i,t,])^2))^(1/2)
      if(t<K)  to[i,t]   = sum(abs(all.ws[[j]][i,t+1,]-all.ws[[j]][i,t,]*((1+rets[nn+t,])/(1+pret[i,t]))))
      model.sd[i,t] = sqrt(t(all.ws[[j]][i,t,])%*%all.cov[[i]][,,t]%*%all.ws[[j]][i,t,])
    }
  }
  all.prets[[j]] = pret
  all.co[[j]]    = co
  all.to[[j]]    = to
  all.portsd[[j]]    = portsd
  all.model.sd[[j]] = model.sd
}

names(all.prets) = names(all.co) = names(all.to) = names(all.portsd) = names(all.model.sd) = names(all.ws)

## GMV
par(mfrow=c(2,3))
for(i in 1:length(models)) {
  plot(all.portsd[[1]][i,],type='l',main=models[i])
  lines(all.model.sd[[1]][i,],col=2)
  }

order(apply(abs(all.portsd[[1]]-all.model.sd[[1]]),1,mean),decreasing=FALSE)
order(apply((all.portsd[[1]]-all.model.sd[[1]])^2,1,mean),decreasing=FALSE)


## CVAR05
par(mfrow=c(2,3))
for(i in 1:length(models)) {
  plot(all.portsd[[2]][i,],type='l',main=models[i])
  lines(all.model.sd[[2]][i,],col=2)
  }

order(apply(abs(all.portsd[[2]]-all.model.sd[[2]]),1,mean),decreasing=FALSE)
order(apply((all.portsd[[2]]-all.model.sd[[2]])^2,1,mean),decreasing=FALSE)

## CVAR10
par(mfrow=c(2,3))
for(i in 1:length(models)) {
  plot(all.portsd[[3]][i,],type='l',main=models[i])
  lines(all.model.sd[[3]][i,],col=2)
  }

order(apply(abs(all.portsd[[3]]-all.model.sd[[3]]),1,mean),decreasing=FALSE)
order(apply((all.portsd[[3]]-all.model.sd[[3]])^2,1,mean),decreasing=FALSE)


# nrow = return+sd+SR+CVAR+CVAR+meansd+CO+TO (8)*3 (portfolios)

lapply(all.prets, function(x) apply(x,1,mean)*252)
lapply(all.prets, function(x) apply(x,1,sd)*sqrt(252))
lapply(all.prets, function(x) apply(sweep(x*252, 2, rf, '-'),1,mean)/(apply(x,1,sd)*sqrt(252)))
lapply(all.prets, function(x) apply(x*252,1,esfun,0.05))
lapply(all.prets, function(x) apply(x*252,1,esfun,0.10))
lapply(all.model.sd, function(x) apply(x*sqrt(252),1,mean))
lapply(all.co, function(x) apply(x,1,mean))
lapply(all.to, function(x) apply(x,1,mean,na.rm=TRUE))

# model ordering:
# 1 - Jore1
# 2 - Geweke
# 3 - EQual
# 4 - CAW
# 5 - DCC-t
# 6 - DCC-HEAVY-t

lapply(lapply(all.prets, function(x) apply(x,1,mean)*252), order,decreasing=FALSE)
lapply(lapply(all.prets, function(x) apply(x,1,sd)*sqrt(252)), order)
lapply(lapply(all.prets, function(x) apply(sweep(x*252, 2, rf, '-'),1,mean)/(apply(x,1,sd)*sqrt(252))), order, decreasing = TRUE)
lapply(lapply(all.prets, function(x) apply(x*252,1,esfun,0.05)), order, decreasing = TRUE)
lapply(lapply(all.prets, function(x) apply(x*252,1,esfun,0.10)), order, decreasing = TRUE)
lapply(lapply(all.model.sd, function(x) apply(x*sqrt(252),1,mean)), order)
lapply(lapply(all.co, function(x) apply(x,1,mean)), order, decreasing = TRUE)
lapply(lapply(all.to, function(x) apply(x,1,mean,na.rm=TRUE)), order, decreasing = TRUE)

# save portfolio return, pertfolio sd, and model produced sds

save(all.model.sd,all.portsd,all.prets,file='empirical/temp/data_mctest_2.Rdata')


# Group all results
# ncol = number of models
# nrow = return+sd+SR+CVAR+CVAR+meansd+TO+CO (8)*3 (portfolios)

ALL.res = matrix(NA,nrow=24,ncol=length(models))

# 1-8 GMV portfolio
# 9-16 CVAR05
# 17-24 CVAR10

# 1 return 2 sd 3 SR 4 CVAR05 5 CVAR10 6 meansd 7 TO 8 CO 

for(i in 1:3){
  ALL.res[i*8-7,] = lapply(all.prets, function(x) apply(x,1,mean)*252)[[i]]
  ALL.res[i*8-6,] = lapply(all.prets, function(x) apply(x,1,sd)*sqrt(252))[[i]]
  ALL.res[i*8-5,] = lapply(all.prets, function(x) apply(sweep(x*252, 2, rf, '-'),1,mean)/(apply(x,1,sd)*sqrt(252)))[[i]]
  ALL.res[i*8-4,] = lapply(all.prets, function(x) apply(x*252,1,esfun,0.05))[[i]]
  ALL.res[i*8-3,] = lapply(all.prets, function(x) apply(x*252,1,esfun,0.10))[[i]]
  ALL.res[i*8-2,] = lapply(all.model.sd, function(x) apply(x*sqrt(252),1,mean))[[i]]
  ALL.res[i*8-1,] = lapply(all.co, function(x) apply(x,1,mean))[[i]]
  ALL.res[i*8,] = lapply(all.to, function(x) apply(x,1,mean,na.rm=TRUE))[[i]]
}



# 1 return 2 sd 3 SR 4 CVAR05 5 CVAR10 6 meansd 7 TO 8 CO 

colnames(ALL.res) = models
rownames(ALL.res) = rep(c('mean','sdev','A.Sh.','CVaR05','CVaR10','avrg stdev','TO','CO'),3)


tableLines <- print(xtable(ALL.res,digits = 3,caption="Portfolio allocation results based on 1-step-ahead
              predictions for 2020/01/02 to 2023/01/31 out-of-sample period ($K$ = 797 observations).
              The three portfolios are:
             Global Minimum Variance (GMV), and minimum Conditional Value at Risk
             for lower 5 and 10 percentiles (CVaR05 and CVaR10), all with short-sale constraints.
             The table reports average portfolio return,
             overall standard deviation (in \\%),
             adjusted Sharpe ratio, Conditional Value at Risk
             for lower 5 and 10 percentiles,
             mean realized standard deviation,
             turnover, and concentration (all quantities annualized) for the
             pooled models (Geweke's, Jore's and
             equally weighted),
              two best individual models (CAW  and DCC-t) and a competitor model (DCC-HEAVY-t).",
                           align = "lccc|ccc",label='table:gmvfull_FX_new'),
                    scalebox=0.8,sanitize.text.function=function(x){x})
writeLines (tableLines, con = "tables_and_figures/gmvfull_FX_new.tex")



