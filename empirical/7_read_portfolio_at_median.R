rm(list=ls(all=TRUE))
library(xtable)

load('empirical/temp/FX_portfolio_at_median_0315.Rdata')

par(mfrow=c(3,2))

plot(ws_gmv[1,,1],type='l',ylim=c(-1,1),main='GMV')
for(i in 1:dm) lines(ws_gmv[1,,i])

plot(ws_gmvr[1,,1],type='l',ylim=c(-1,1),main='GMVr')
for(i in 1:dm) lines(ws_gmvr[1,,i])

plot(ws_CVAR10[1,,1],type='l',ylim=c(-1,1),main='CVAR10')
for(i in 1:dm) lines(ws_CVAR10[1,,i])

plot(ws_CVAR10r[1,,1],type='l',ylim=c(-1,1),main='CVAR10')
for(i in 1:dm) lines(ws_CVAR10r[1,,i])

plot(ws_CVAR05[1,,1],type='l',ylim=c(-1,1),main='CVAR05')
for(i in 1:dm) lines(ws_CVAR05[1,,i])

plot(ws_CVAR05r[1,,1],type='l',ylim=c(-1,1),main='CVAR05')
for(i in 1:dm) lines(ws_CVAR05r[1,,i])

par(mfrow=c(2,3))
for(i in 1:length(models)){
  plot(apply((ws_gmv[i,,]<0)*ws_gmv[i,,],1,sum),type='l',ylim=c(-0.8,-0.3),
       main=models[i])
}

par(mfrow=c(2,3))
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

##

all.ws = list(ws_gmv,ws_gmvr,
              ws_CVAR05,ws_CVAR05r,ws_CVAR10,ws_CVAR10r)


names(all.ws) = c('gmv','gmvr','CVAR05','CVAR05r','CVAR10','CVAR10r')

all.prets = list()
pret= array(NA,dim=c(length(models),K))

for(j in 1:length(all.ws)){
  for(i in 1:length(models)){
    for(t in 1:K){
      pret[i,t] = sum(all.ws[[j]][i,t,]*rets[nn+t,])
    }
  }
  all.prets[[j]] = pret
}

names(all.prets) = c('gmv','gmvr','CVAR05','CVAR05r','CVAR10','CVAR10r')

lapply(all.prets, function(x) apply(x,1,sd)*sqrt(252))
lapply(all.prets, function(x) apply(sweep(x*252, 2, rf, '-'),1,mean)/(apply(x,1,sd)*sqrt(252)))
lapply(all.prets, function(x) apply(x*252,1,esfun,0.05))
lapply(all.prets, function(x) apply(x*252,1,esfun,0.10))





# 
# 
# ALL.res = matrix(NA,ncol=16,nrow=length(models))
# 
# ### portfolio ret
# 
# p.ret.gmv =p.ret.cvar05 =p.ret.cvar10 =p.ret.minvar= array(NA,dim=c(length(models),K))
# 
# for(i in 1:length(models)){
#   for(t in 1:K){
#     p.ret.gmv[i,t]    = sum(ws_gmv[i,t,]*rets[nn+t,])
#     p.ret.cvar05[i,t] = sum(ws_CVAR05[i,t,]*rets[nn+t,])
#     p.ret.cvar10[i,t] = sum(ws_CVAR10[i,t,]*rets[nn+t,])
#     p.ret.minvar[i,t] = sum(ws_gmvr[i,t,]*rets[nn+t,])
#   }
# }
# 
# ### stdev annualized
# 
# ALL.res[,1] = apply(p.ret.gmv, 1,sd)*sqrt(252)
# ALL.res[,5] = apply(p.ret.cvar05, 1,sd)*sqrt(252)
# ALL.res[,9] = apply(p.ret.cvar10, 1,sd)*sqrt(252)
# ALL.res[,13] = apply(p.ret.minvar, 1,sd)*sqrt(252)
# 
# 
# ### Adj.Sharpe annualized
# 
# ALL.res[,2]  = apply(sweep(p.ret.gmv*252, 2, rf, '-'),1,mean)/(apply(p.ret.gmv,1,sd)*sqrt(252))
# ALL.res[,6]  = apply(sweep(p.ret.cvar05*252, 2, rf, '-'),1,mean)/(apply(p.ret.cvar05,1,sd)*sqrt(252))
# ALL.res[,10] = apply(sweep(p.ret.cvar10*252, 2, rf, '-'),1,mean)/(apply(p.ret.cvar10,1,sd)*sqrt(252))
# ALL.res[,14] = apply(sweep(p.ret.minvar*252, 2, rf, '-'),1,mean)/(apply(p.ret.minvar,1,sd)*sqrt(252))
# 
# 
# ### CVAR 5%
# 
# ALL.res[,3]  = apply(p.ret.gmv*252,1,esfun,0.05)
# ALL.res[,7]  = apply(p.ret.cvar05*252,1,esfun,0.05)
# ALL.res[,11] = apply(p.ret.cvar10*252,1,esfun,0.05)
# ALL.res[,15] = apply(p.ret.minvar*252,1,esfun,0.05)
# 
# 
# ### CVAR 10%
# 
# ALL.res[,4]  = apply(p.ret.gmv*252,1,esfun,0.1)
# ALL.res[,8]  = apply(p.ret.cvar05*252,1,esfun,0.1)
# ALL.res[,12] = apply(p.ret.cvar10*252,1,esfun,0.1)
# ALL.res[,16] = apply(p.ret.minvar*252,1,esfun,0.1)
# 
# 
# rownames(ALL.res) = models
# colnames(ALL.res) = rep(c('sd','A.Sh.','CVaR05','CVaR10'),4)
# 
# round(ALL.res,2)

# 
# 
# tableLines <- print(xtable(ALL.res,digits = 3,caption="Portfolio allocation results based on 1-step-ahead 
#               predictions for 2020/01/02 to 2023/01/31 out-of-sample period ($K$ = 797 observations).
#               The four portfolios are:
#              the unrestricted Global Minimum Variance (GMV), minimum Conditional Value at Risk
#              for lower 5 and 10 percentiles (CVaR05 and CVaR10), and restricted GMV with a short-sale constraint.
#              The table reports portfolio standard deviation (in \\%),
#              adjusted Sharpe ratio, and realized Conditional Value at Risk
#              for lower 5 and 10 percentiles (all quantities annualized) for the 
#              pooled models (Geweke's, Jore's and
#              equally weighted),
#               two best individual models (CAW  and DCC-t) and a competitor model (DCC-HEAVY-t).",
#                            align = "lcccc|cccc|cccc|cccc",label='table:gmvfull_FX'),
#                     scalebox=0.6,sanitize.text.function=function(x){x})
# multicolumns <- "& \\\\multicolumn{4}{c}{GMV}
#                  & \\\\multicolumn{4}{c}{CVaR05}
#                  & \\\\multicolumn{4}{c}{CVaR10}
#                  & \\\\multicolumn{4}{c}{GMVrestr.}  \\\\\\\\"
# tableLines <- sub ("\\\\toprule\\n", paste0 ("\\\\toprule\n", multicolumns, "\n"), tableLines) ## booktabs = TRUE
# tableLines <- sub ("\\\\hline\\n",   paste0 ("\\\\hline\n",   multicolumns, "\n"), tableLines) ## booktabs = FALSE
# writeLines (tableLines, con = "tables_and_figures/gmvfull_FX.tex")
# 
# 


# # 1. jore's1
# # 2. geweke's
# # 3. equally w.
# # 4. caw
# # 5. dcct
# # 6. dcc-heavy-t
# 
# ####----------------
# ## Turnover
# ####----------------
# 
# # > dim(my.weights)
# # [1]    6  797   14
# 
# TO = array(NA,dim=c(length(models),K))
# # > dim(TO)
# # [1]    6  797
# 
# for(i in 1:length(models)){
#   for(t in 2:K){
#     TO[i,t] = sum(abs(my.weights[i,t,]-my.weights[i,t-1,]*(1+rets[nn+t-1,])/(1+sum(my.weights[i,t-1,]*rets[nn+t-1,]))))
#   }
#   print(i)
# }
# 
# apply(TO,1,mean,na.rm=TRUE)
# 
# 
# ####----------------
# ## Concentration
# ####----------------
# 
# # > dim(my.weights)
# # [1]    6  797   14
# 
# CO = array(NA,dim=c(length(models),K))
# # > dim(CO)
# # [1]    6  797
# 
# for(i in 1:length(models)){
#   for(t in 1:K){
#     CO[i,t] = sqrt(sum(my.weights[i,t,]^2))
#   }
#   print(i)
# }
# 
# apply(CO,1,mean,na.rm=TRUE)
# 
# ####----------------
# ## Short Position
# ####----------------
# 
# # > dim(my.weights)
# # [1]    6  797   14
# 
# SP = array(NA,dim=c(length(models),K))
# # > dim(TO)
# # [1]    6  797
# 
# for(i in 1:length(models)){
#   for(t in 1:K){
#     SP[i,t] = sum((my.weights[i,t,]<0)* my.weights[i,t,])
#   }
#   print(i)
# }
# 
# apply(SP,1,mean,na.rm=TRUE)




