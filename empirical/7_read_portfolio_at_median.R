rm(list=ls(all=TRUE))

load('empirical/temp/FX_portfolio_at_median_0315.Rdata')

plot(ws_gmv[1,,1],type='l',ylim=c(-1,1),main='GMV')
for(i in 1:dm) lines(ws_gmv[1,,i])

plot(ws_minvar[1,,1],type='l',ylim=c(-1,1),main='GMVr')
for(i in 1:dm) lines(ws_minvar[1,,i])

plot(ws_CVAR10[1,,1],type='l',ylim=c(-1,1),main='CVAR10')
for(i in 1:dm) lines(ws_CVAR10[1,,i])

plot(ws_CVAR05[1,,1],type='l',ylim=c(-1,1),main='CVAR05')
for(i in 1:dm) lines(ws_CVAR05[1,,i])


############################################################

esfun=function(x,p){
  es=mean(x[which(x<quantile(x,p))])
  return(es)
}
load('data/rf.RData')

##

ALL.res = matrix(NA,ncol=16,nrow=length(models))

### portfolio ret

p.ret.gmv =p.ret.cvar05 =p.ret.cvar10 =p.ret.minvar= array(NA,dim=c(length(models),K))

for(i in 1:length(models)){
  for(t in 1:K){
    p.ret.gmv[i,t]    = sum(ws_gmv[i,t,]*rets[nn+t,])
    p.ret.cvar05[i,t] = sum(ws_CVAR05[i,t,]*rets[nn+t,])
    p.ret.cvar10[i,t] = sum(ws_CVAR10[i,t,]*rets[nn+t,])
    p.ret.minvar[i,t] = sum(ws_minvar[i,t,]*rets[nn+t,])
  }
}

### stdev annualized

ALL.res[,1] = apply(p.ret.gmv, 1,sd)*sqrt(252)
ALL.res[,5] = apply(p.ret.cvar05, 1,sd)*sqrt(252)
ALL.res[,9] = apply(p.ret.cvar10, 1,sd)*sqrt(252)
ALL.res[,13] = apply(p.ret.minvar, 1,sd)*sqrt(252)


### Adj.Sharpe annualized

ALL.res[,2]  = apply(sweep(p.ret.gmv*252, 2, rf, '-'),1,mean)/(apply(p.ret.gmv,1,sd)*sqrt(252))
ALL.res[,6]  = apply(sweep(p.ret.cvar05*252, 2, rf, '-'),1,mean)/(apply(p.ret.cvar05,1,sd)*sqrt(252))
ALL.res[,10] = apply(sweep(p.ret.cvar10*252, 2, rf, '-'),1,mean)/(apply(p.ret.cvar10,1,sd)*sqrt(252))
ALL.res[,14] = apply(sweep(p.ret.minvar*252, 2, rf, '-'),1,mean)/(apply(p.ret.minvar,1,sd)*sqrt(252))


### CVAR 5%

ALL.res[,3]  = apply(p.ret.gmv*252,1,esfun,0.05)
ALL.res[,7]  = apply(p.ret.cvar05*252,1,esfun,0.05)
ALL.res[,11] = apply(p.ret.cvar10*252,1,esfun,0.05)
ALL.res[,15] = apply(p.ret.minvar*252,1,esfun,0.05)


### CVAR 10%

ALL.res[,4]  = apply(p.ret.gmv*252,1,esfun,0.1)
ALL.res[,8]  = apply(p.ret.cvar05*252,1,esfun,0.1)
ALL.res[,12] = apply(p.ret.cvar10*252,1,esfun,0.1)
ALL.res[,16] = apply(p.ret.minvar*252,1,esfun,0.1)


rownames(ALL.res) = models

round(ALL.res,2)

# 1. jore's1
# 2. geweke's
# 3. equally w.
# 4. caw
# 5. dcct
# 6. dcc-heavy-t

####----------------
## Turnover
####----------------

# > dim(my.weights)
# [1]    6  797   14

TO = array(NA,dim=c(length(models),K))
# > dim(TO)
# [1]    6  797

for(i in 1:length(models)){
  for(t in 2:K){
    TO[i,t] = sum(abs(my.weights[i,t,]-my.weights[i,t-1,]*(1+rets[nn+t-1,])/(1+sum(my.weights[i,t-1,]*rets[nn+t-1,]))))
  }
  print(i)
}

apply(TO,1,mean,na.rm=TRUE)


####----------------
## Concentration
####----------------

# > dim(my.weights)
# [1]    6  797   14

CO = array(NA,dim=c(length(models),K))
# > dim(CO)
# [1]    6  797

for(i in 1:length(models)){
  for(t in 1:K){
    CO[i,t] = sqrt(sum(my.weights[i,t,]^2))
  }
  print(i)
}

apply(CO,1,mean,na.rm=TRUE)

####----------------
## Short Position
####----------------

# > dim(my.weights)
# [1]    6  797   14

SP = array(NA,dim=c(length(models),K))
# > dim(TO)
# [1]    6  797

for(i in 1:length(models)){
  for(t in 1:K){
    SP[i,t] = sum((my.weights[i,t,]<0)* my.weights[i,t,])
  }
  print(i)
}

apply(SP,1,mean,na.rm=TRUE)







# all.res.gmv = list()
# for(i in 1:length(models)){
#   # var05     = apply(gvm_ret[i,,],1,quantile,0.05)
#   # var10     = apply(gvm_ret[i,,],1,quantile,0.10)
#   # es05      = apply(gvm_ret[i,,],1,esfun,0.05)
#   # es10      = apply(gvm_ret[i,,],1,esfun,0.10)
#   psd       = apply(gvm_ret[i,,],1,sd)
#   GL        = (100*sqrt(252)*(apply(gvm_ret[4,,],1,sd)-apply(gvm_ret[i,,],1,sd))/
#                  apply(gvm_ret[i,,],1,sd))
#   shr       = apply(gvm_ret[i,,-1],1,sum)/(apply(gvm_ret[i,,-1],1,sd)*sqrt(252))
#   all.res.gmv   = c(all.res.gmv,list(data.frame(var05,var10,es05,es10,GL,shr,psd)))
# }
# 
# 
# 
# res.gmv = NULL
# 
# for(i in 1:length(models)){
#   res.gmv = rbind(res.gmv,c(quantile(all.res.gmv[[i]]$GL,0.05),
#                             median(all.res.gmv[[i]]$GL),
#                             quantile(all.res.gmv[[i]]$GL,0.95),
#                             quantile(all.res.gmv[[i]]$psd*100,0.05),
#                             median(all.res.gmv[[i]]$psd*100),
#                             quantile(all.res.gmv[[i]]$psd*100,0.95)))
# }
# 
# res.gmv
# 
# 
# res  = res.gmv[c(2,1,3,4,5,6),]
# 
# rownames(res ) = models[c(2,1,3,4,5,6)]
# colnames(res ) = c('P05','Median','P95','P05','Median','P95')
# 
# round(res ,3)
# 
# rm(reshf,Sig,Sigma,reslf,resH)
# 
# save.image('temp/portfolio_EX.Rdata')
# 
# 
# tableLines <- print(xtable(res,digits = 3,caption="GMV portfolio results based on 1-step-ahead predictions 
#              for  2021/01/04-2021/12/31 out-of-sample period 
#              ($K=252$ observations) for 5-variate dataset.
#              The table reports the posterior 5, 50 and 95 percentiles of G/L criteria as well as
#              portfolio standard deviation (in \\%) for the pooled models (Geweke's, Jore's and 
#              equally weighted), 
#               two best individual models (Additive Inverse Wishart  and 
#              Dynamic Conditional Correlation with $t$ copula) and a competitor model (DCC-HEAVY-t).",
#                            align = "lccc|ccc",label='table:gmvfull_EX'), 
#                     scalebox=0.8,sanitize.text.function=function(x){x})
# multicolumns <- "& \\\\multicolumn{3}{c}{G/L}
#                  & \\\\multicolumn{3}{c}{Portfolio stdev.}  \\\\\\\\"
# tableLines <- sub ("\\\\toprule\\n", paste0 ("\\\\toprule\n", multicolumns, "\n"), tableLines) ## booktabs = TRUE
# tableLines <- sub ("\\\\hline\\n",   paste0 ("\\\\hline\n",   multicolumns, "\n"), tableLines) ## booktabs = FALSE
# writeLines (tableLines, con = "tables_and_figures/gmvfull_EX.tex")
# 
# 
# 
