rm(list=ls(all=TRUE))
load('empirical/temp/res_FX.Rdata')
load('data/FXdata.Rdata')
library(xtable)
library(Rfast)
library(mvnfast)
library(profvis)
library(NMOF)
library(PortfolioAnalytics)

load('empirical/temp/marginals.Rdata')
load('empirical/temp/RV_forc.Rdata')

# par(mfrow=c(3,5))
# for(i in 1:dm){
#   plot(tail(RVs[,i],K),type='l')
#   lines(RV_forc$'1sa'[,i],lwd=2,col=2)
# }

# marginals$rvs=tail(RVs,K)
marginals$rvs=RV_forc$'1sa'

load('empirical/temp/results_vectordcc_tcop.Rdata')
reslf = res
load('empirical/temp/results_caw_iw.Rdata')
reshf = res
load('empirical/temp/results_heavy_t.Rdata')
resH  = res

rm(res)

MCMCsize=1000

# 1. jore's1
# 2. geweke's
# 3. equally w.
# 4. caw
# 5. dcct
# 6. dcc-heavy-t
models = c('Jore1','Geweke','Equal','CAW','DCC-t','DCC-HEAVY-t')
ws_gmv = ws_CVAR05 = ws_CVAR10 = array(NA,c(length(models),K,dm))

# for Low Frequency (rank-1 DCC-t) model
Q  = R = array(NA,c(dm, dm, nn+K))
a  <- apply(reslf$r[ind,2:(dm+1)],2,median)
b  <- apply(reslf$r[ind,(dm+2):(2*dm+1)],2,median)
nu <- median(reslf$r[ind,1])
iota    = rep(1,dm)
Oiota   = Outer(iota,iota)

# for DCC-HEAVY-t model
Rh      = array(NA,c(dm, dm, nn+K))
nuh     <- median(resH$r[ind,1])
ah      <- median(resH$r[ind,2])
bh      <- median(resH$r[ind,3])
Pbar = Reduce('+',Sigma[1:nn])/nn

# for High Frequency (CAW) model
nux  = mean(reshf$r[ind,1])
b1   = apply(reshf$r[ind,2:(dm+1)],2,median)
b2   = apply(reshf$r[ind,(dm+2):(2*dm+1)],2,median)
Rcaw = array(NA,c(dm, dm, nn+K))
Rcaw[,,1] = Pbar

# some other initials, median weights
US  = matrix(runif(3*(nn+K)),ncol=nn+K,nrow=3)
A.tmp = B.tmp = C.tmp = matrix(NA,nrow=MCMCsize,ncol=dm)
class(A.tmp) <- "numeric"
class(B.tmp) <- "numeric"
class(C.tmp) <- "numeric"

wj1 = apply(ws_jore1,2,median)
wgw = apply(ws_gew,2,median)

par(mfrow=c(1,1))
plot(wj1,type='l')
lines(wgw,col=2,lwd=2)

# for DCC-t model
tdata <- qt(udata,nu)
Sbar  <- cova(tdata[1:nn,])
A     <- Outer(a,a)
B     <- Outer(b,b)
B0lf  <- (Oiota-A-B)*Sbar
INL   <- dt(tdata,df=nu,log=TRUE)
Q[,,1] <- cora(tdata[1:nn,])

# for DCC-HEAVY-t model
tdatah  <- qt(udata,nuh)
Rbar   <- cora(tdatah[1:nn,])
Rh[,,1] = Rbar

# CAW
B1  = Outer(b1,b1)
B2  = Outer(b2,b2)
B0  = (Oiota-B1-B2)*Pbar
mtdf = nux-dm

for(t in 2:(nn+K)){
  
  t0 = Sys.time()

  # for DCC-t model
  Q[,,t] <- B0lf+A*Outer(tdata[t-1,],tdata[t-1,])+B*Q[,,(t-1)]
  t.ma  = Q[,,t]
  t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
  t.R   = Outer(t.dv,t.dv)*t.ma
  R[,,t] = t.R
  
  # for DCC-HEAVY-t model
  Rh[,,t]   <- (1-bh)*Rbar-ah*Pbar+ah*Sigma[[t-1]]+bh*Rh[,,t-1]
  
  # for CAW model
  Rcaw[,,t] = B0+B1*Sigma[[t-1]]+B2*Rcaw[,,t-1]
  
  if(t>nn){

    mvnfast::rmvt(MCMCsize,rep(0,dm),(mtdf-1)/(mtdf+1)*Rcaw[,,t],df=mtdf+1,A=A.tmp)
    sample_uxm = pt(A.tmp,df = mtdf+1)

    mvnfast::rmvt(MCMCsize,rep(0,dm),R[,,t],df=nu,A=B.tmp)
    sample_udcct = pt(B.tmp,df = nu)

    mvnfast::rmvt(MCMCsize,rep(0,dm),Rh[,,t],df=nuh,A=C.tmp)
    sample_udccth = pt(C.tmp,df = nuh)

    sample_standretxm = qnorm(sample_uxm)
    sample_standretdcct = qnorm(sample_udcct)
    sample_standretdccth = qnorm(sample_udccth)

    sample_retsxm = ((t(sample_standretxm)*marginals$sd)+marginals$mean)*marginals$rvs[t-nn,]
    sample_retsdcct = ((t(sample_standretdcct)*marginals$sd)+marginals$mean)*marginals$rvs[t-nn,]
    sample_retsdccth = ((t(sample_standretdccth)*marginals$sd)+marginals$mean)*marginals$rvs[t-nn,]

    ##################

    invS = solve(cova(t(sample_retsxm)))
    nom  = invS%*%iota
    den  = as.vector(t(iota)%*%invS%*%iota)
    pws  = nom/den
    ws_gmv[4,t-nn,]= pws
    
    tmp = t(sample_retsxm)
    res <- minCVaR(tmp, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10[4,t-nn,]=c(res)
    res <- minCVaR(tmp, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05[4,t-nn,]=c(res)
    
    #
    invS = solve(cova(t(sample_retsdcct)))
    nom  = invS%*%iota
    den  = as.vector(t(iota)%*%invS%*%iota)
    pws  = nom/den
    ws_gmv[5,t-nn,]= pws
    
    tmp = t(sample_retsdcct)
    res <- minCVaR(tmp, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10[5,t-nn,]=c(res)
    res <- minCVaR(tmp, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05[5,t-nn,]=c(res)

    #
    invS = solve(cova(t(sample_retsdccth)))
    nom  = invS%*%iota
    den  = as.vector(t(iota)%*%invS%*%iota)
    pws  = nom/den
    ws_gmv[6,t-nn,]= pws

    tmp = t(sample_retsdccth)
    res <- minCVaR(tmp, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10[6,t-nn,]=c(res)
    res <- minCVaR(tmp, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05[6,t-nn,]=c(res)
    
    #
    if(wj1[t-nn]>US[1,t-nn]) selected = sample_retsxm else selected = sample_retsdcct

    invS = solve(cova(t(selected)))
    nom  = invS%*%iota
    den  = as.vector(t(iota)%*%invS%*%iota)
    pws  = nom/den
    ws_gmv[1,t-nn,]= pws
    
    tmp = t(selected)
    res <- minCVaR(tmp, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10[1,t-nn,]=c(res)
    res <- minCVaR(tmp, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05[1,t-nn,]=c(res)

    if(wgw[t-nn]>US[2,t-nn]) selected = sample_retsxm else selected = sample_retsdcct

    invS = solve(cova(t(selected)))
    nom  = invS%*%iota
    den  = as.vector(t(iota)%*%invS%*%iota)
    pws  = nom/den
    ws_gmv[2,t-nn,]= pws
    
    tmp = t(selected)
    res <- minCVaR(tmp, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10[2,t-nn,]=c(res)
    res <- minCVaR(tmp, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05[2,t-nn,]=c(res)

    if(0.5>US[3,t-nn]) selected = sample_retsxm else selected = sample_retsdcct

    invS = solve(cova(t(selected)))
    nom  = invS%*%iota
    den  = as.vector(t(iota)%*%invS%*%iota)
    pws  = nom/den
    ws_gmv[3,t-nn,]= pws
    
    
    ###
    tmp = data.frame(t(selected))
    colnames(tmp) = assets
    portfolio <- portfolio.spec(assets = assets)
    portfolio <- add.constraint(portfolio, type = "weight", min = 0)
    portfolio <- add.objective(portfolio, type = "risk", name = "var")
    opt_results <- optimize.portfolio(R = tmp, portfolio = portfolio)
    
    optimize.portfolio(tmp, portfolio = NULL, constraints = NULL,
                       objectives = NULL)
    
    
    ###
    
    tmp = t(selected)
    res <- minCVaR(tmp, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10[3,t-nn,]=c(res)
    res <- minCVaR(tmp, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05[3,t-nn,]=c(res)
  }
  print(t)
  print(Sys.time()-t0)
}

p1=1
p2=2

par(mfrow=c(1,1))
plot(tail(RCor[p1,p2,],K),col='gray80',lwd=1,type='l',ylim=c(-1,1))
lines(tail(R[p1,p2,],K),lwd=2)
lines(tail(Rh[p1,p2,],K),col=2,lwd=2)
lines(tail(Rcaw[p1,p2,],K),col=3,lwd=2)

save.image(file = 'empirical/temp/FX_portfolio_at_median.Rdata')

plot(ws_CVAR10[1,,1],type='l',ylim=c(-1,1))
for(i in 1:dm) lines(ws_CVAR10[1,,i])

plot(ws_CVAR05[1,,1],type='l',ylim=c(-1,1))
for(i in 1:dm) lines(ws_CVAR05[1,,i])

plot(ws_gmv[1,,1],type='l',ylim=c(-1,1))
for(i in 1:dm) lines(ws_gmv[1,,i])

############################################################

esfun=function(x,p){
  es=mean(x[which(x<quantile(x,p))])
  return(es)
}
load('data/rf.RData')

##

ALL.res = matrix(NA,ncol=12,nrow=length(models))

### portfolio ret

p.ret.gmv =p.ret.cvar05 =p.ret.cvar10 = array(NA,dim=c(length(models),K))

for(i in 1:length(models)){
  for(t in 1:K){
    p.ret.gmv[i,t]    = sum(ws_gmv[i,t,]*rets[nn+t,])
    p.ret.cvar05[i,t] = sum(ws_CVAR05[i,t,]*rets[nn+t,])
    p.ret.cvar10[i,t] = sum(ws_CVAR10[i,t,]*rets[nn+t,])
  }
}

### stdev annualized

ALL.res[,1] = apply(p.ret.gmv, 1,sd)*sqrt(252)
ALL.res[,5] = apply(p.ret.cvar05, 1,sd)*sqrt(252)
ALL.res[,9] = apply(p.ret.cvar10, 1,sd)*sqrt(252)

### Adj.Sharpe annualized

ALL.res[,2]  = apply(sweep(p.ret.gmv*252, 2, rf, '-'),1,mean)/(apply(p.ret.gmv,1,sd)*sqrt(252))
ALL.res[,6]  = apply(sweep(p.ret.cvar05*252, 2, rf, '-'),1,mean)/(apply(p.ret.cvar05,1,sd)*sqrt(252))
ALL.res[,10] = apply(sweep(p.ret.cvar10*252, 2, rf, '-'),1,mean)/(apply(p.ret.cvar10,1,sd)*sqrt(252))

### CVAR 5%

ALL.res[,3]  = apply(p.ret.gmv*252,1,esfun,0.05)
ALL.res[,7]  = apply(p.ret.cvar05*252,1,esfun,0.05)
ALL.res[,11] = apply(p.ret.cvar10*252,1,esfun,0.05)

### CVAR 10%

ALL.res[,4]  = apply(p.ret.gmv*252,1,esfun,0.1)
ALL.res[,8]  = apply(p.ret.cvar05*252,1,esfun,0.1)
ALL.res[,12] = apply(p.ret.cvar10*252,1,esfun,0.1)

rownames(ALL.res) = models

ALL.res

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
