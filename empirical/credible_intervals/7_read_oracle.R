
rm(list=ls(all=TRUE))
library(xtable)

load('empirical/temp/FX_portfolio_true.Rdata')

# compare weights for model 1
par(mfrow=c(3,1))

plot(ws_gmvr[1,,1],type='l',ylim=c(-1,1),main='GMVr')
for(i in 1:dm) lines(ws_gmvr[1,,i])

plot(ws_CVAR10r[1,,1],type='l',ylim=c(-1,1),main='CVAR10r')
for(i in 1:dm) lines(ws_CVAR10r[1,,i])

plot(ws_CVAR05r[1,,1],type='l',ylim=c(-1,1),main='CVAR05r')
for(i in 1:dm) lines(ws_CVAR05r[1,,i])

############################################################

esfun=function(x,p){
  es=mean(x[which(x<quantile(x,p))])
  return(es)
}

load('data/rf.RData')
load('data/FXdata.Rdata')

## var-covar matrices

all.ws  = list(ws_gmvr,ws_CVAR05r,ws_CVAR10r)

names(all.ws) = c('gmvr','CVAR05r','CVAR10r')

oracle.portsd = list()
portsd = array(NA,dim=c(1,K))

for(j in 1:length(all.ws)){
    for(t in 1:K){
      portsd[1,t] = sqrt(t(all.ws[[j]][1,t,])%*%RCov[,,nn+t]%*%all.ws[[j]][1,t,])
    }
  oracle.portsd[[j]]    = portsd
}

save(oracle.portsd,file='empirical/temp/oraclesd.Rdata')

