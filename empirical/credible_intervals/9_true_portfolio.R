rm(list=ls(all=TRUE))
load('empirical/temp/res_FX.Rdata')
load('data/FXdata.Rdata')
library(xtable)
library(Rfast)
library(mvnfast)
library(profvis)
library(NMOF)

MCMCsize=10000

ws_gmvr = ws_CVAR05r = ws_CVAR10r = array(NA,c(1,K,dm))

for(t in 2:(nn+K)){
  
  t0 = Sys.time()
  
  if(t>nn){
    
    sample_rets = mvnfast::rmvn(MCMCsize,rep(0,dm),RCov[,,t])

    ##################

    res <- minvar(RCov[,,t], wmin = 0, wmax = 1)
    ws_gmvr[1,t-nn,]=c(res)
    
    res <- minCVaR(sample_rets, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10r[1,t-nn,]=c(res)
    
    res <- minCVaR(sample_rets, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05r[1,t-nn,]=c(res)
    
  }
  print(t)
  print(Sys.time()-t0)
}




rm(reshf,reslf,resH,Sigma,ws_DN,ws_gew,ws_jore1,Q,R,Rcaw,RCor,RCov,Rh,A.tmp,B.tmp,C.tmp,
   sample_retsdcct,sample_retsdccth,sample_retsxm,sample_standretdcct,sample_standretdccth,
   sample_standretxm,sample_udcct,sample_udccth,sample_uxm,selected,tmp)

save.image(file = 'empirical/temp/FX_portfolio_true.Rdata')
