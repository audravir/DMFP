rm(list=ls(all=TRUE))
load('empirical/temp/res_FX.Rdata')
load('data/FXdata.Rdata')
library(xtable)
library(Rfast)
library(mvnfast)
library(profvis)

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

MCMCsize=5000

# 1. jore's1
# 2. geweke's
# 3. equally w.
# 4. caw
# 5. dcct
# 6. dcc-heavy-t
models = c('Jore1','Geweke','Equal','CAW','DCC-t','DCC-HEAVY-t')
ws_gmv = array(NA,c(length(models),post.sample,K,dm))

# for Low Frequency (rank-1 DCC-t) model
Q  = R = array(NA,c(dm, dm, nn+K))
a  <- reslf$r[ind,2:(dm+1)]
b  <- reslf$r[ind,(dm+2):(2*dm+1)]
nu <- reslf$r[ind,1]
iota    = rep(1,dm)
Oiota   = Outer(iota,iota)

# for DCC-HEAVY-t model
Rh      = array(NA,c(dm, dm, nn+K))
nuh     <- resH$r[ind,1]
ah      <- resH$r[ind,2]
bh      <- resH$r[ind,3]
Pbar = Reduce('+',Sigma[1:nn])/nn

# for High Frequency (CAW) model
nux  = reshf$r[ind,1]
b1   = reshf$r[ind,2:(dm+1)]
b2   = reshf$r[ind,(dm+2):(2*dm+1)]
Rcaw = array(NA,c(dm, dm, nn+K))
Rcaw[,,1] = Pbar

US  = matrix(runif(3*(nn+K)),ncol=nn+K,nrow=3)
A.tmp = B.tmp = C.tmp = matrix(NA,nrow=MCMCsize,ncol=dm)
class(A.tmp) <- "numeric"
class(B.tmp) <- "numeric"
class(C.tmp) <- "numeric"

for(m in 1:post.sample){
  t0 = Sys.time()
  
  # for DCC-t model
  tdata <- qt(udata,nu[m])
  Sbar  <- cova(tdata[1:nn,])
  A     <- Outer(a[m,],a[m,])
  B     <- Outer(b[m,],b[m,])
  B0lf  <- (Oiota-A-B)*Sbar
  INL   <- dt(tdata,df=nu[m],log=TRUE)
  Q[,,1] <- cora(tdata[1:nn,])
  
  # for DCC-HEAVY-t model
  tdatah  <- qt(udata,nuh[m])
  Rbar   <- cora(tdatah[1:nn,])
  Rh[,,1] = Rbar
  
  # CAW
  B1  = Outer(b1[m,],b1[m,])
  B2  = Outer(b2[m,],b2[m,])
  B0  = (Oiota-B1-B2)*Pbar
  mtdf = nux[m]-dm

  for(t in 2:(nn+K)){
    # for DCC-t model
    Q[,,t] <- B0lf+A*Outer(tdata[t-1,],tdata[t-1,])+B*Q[,,(t-1)]
    t.ma  = Q[,,t]
    t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
    t.R   = Outer(t.dv,t.dv)*t.ma
    R[,,t] = t.R
    
    # for DCC-HEAVY-t model
    Rh[,,t]   <- (1-bh[m])*Rbar-ah[m]*Pbar+ah[m]*Sigma[[t-1]]+bh[m]*Rh[,,t-1]

    # for CAW model
    Rcaw[,,t] = B0+B1*Sigma[[t-1]]+B2*Rcaw[,,t-1]

    if(t>nn){

      mvnfast::rmvt(MCMCsize,rep(0,dm),(mtdf-1)/(mtdf+1)*Rcaw[,,t],df=mtdf+1,A=A.tmp)
      sample_uxm = pt(A.tmp,df = mtdf+1)

      mvnfast::rmvt(MCMCsize,rep(0,dm),R[,,t],df=nu[m],A=B.tmp)
      sample_udcct = pt(B.tmp,df = nu[m])

      mvnfast::rmvt(MCMCsize,rep(0,dm),Rh[,,t],df=nuh[m],A=C.tmp)
      sample_udccth = pt(C.tmp,df = nuh[m])

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
      ws_gmv[4,m,t-nn,]= pws

      invS = solve(cova(t(sample_retsdcct)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den
      ws_gmv[5,m,t-nn,]= pws

      invS = solve(cova(t(sample_retsdccth)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den
      ws_gmv[6,m,t-nn,]= pws

      if(ws_jore1[m,t-nn]>runif(1))  selected = sample_retsxm else selected = sample_retsdcct
      # if(ws_jore1[m,t-nn]>US[1,t-nn]) selected = sample_retsxm else selected = sample_retsdcct

      invS = solve(cova(t(selected)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den
      ws_gmv[1,m,t-nn,]= pws

      if(ws_gew[m,t-nn]>runif(1)) selected = sample_retsxm else selected = sample_retsdcct
      # if(ws_gew[m,t-nn]>US[2,t-nn]) selected = sample_retsxm else selected = sample_retsdcct

      invS = solve(cova(t(selected)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den
      ws_gmv[2,m,t-nn,]= pws

      if(0.5>runif(1)) selected = sample_retsxm else selected = sample_retsdcct
      # if(0.5>US[3,t-nn]) selected = sample_retsxm else selected = sample_retsdcct

      invS = solve(cova(t(selected)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den
      ws_gmv[3,m,t-nn,]= pws
    }
  }
  print(c(m,Sys.time()-t0))
}


save.image(file = 'empirical/temp/FX_portfolio_march12.Rdata')


###-----------
### GMV portfolio: mean, var, sharpe, TO,CO, SP, VaR, ES
###-----------
esfun=function(x,p){
  es=mean(x[which(x<quantile(x,p))])
  return(es)
}

gvm_ret = array(NA,dim=c(length(models),post.sample,K))

for(m in 1:length(ind)){
   for(i in 1:length(models)){
     for(t in 1:K){
       gvm_ret[i,m,t] = sum(ws_gmv[i,m,t,]*rets[nn+t,])
       }
   }
}

####----------------
## Turnover
####----------------

# > dim(ws_gmv)
# [1]    6 1000  797   14

TO = array(NA,dim=c(length(models),post.sample,K))
# > dim(TO)
# [1]    6 1000  797

for(i in 1:length(models)){
  for(m in 1:post.sample){
    for(t in 2:K){
      TO[i,m,t] = sum(abs(ws_gmv[i,m,t,]-ws_gmv[i,m,t-1,]*(1+rets[nn+t-1,])/(1+sum(ws_gmv[i,m,t-1,]*rets[nn+t-1,]))))
    }
  }
  print(i)
}

apply((apply(TO,c(1,2),mean,na.rm=TRUE)),1,median)


# at the mean
m_ws = apply(ws_gmv,c(1,3,4),median)
# > dim(m_ws)
# [1]   6 797  14

m_TO = array(NA,dim=c(length(models),K))

for(i in 1:length(models)){
  for(t in 1:(K-1)){
    m_TO[i,t+1]=sum(abs(m_ws[i,t+1,]-m_ws[i,t,]*(1+rets[nn+t,])/(1+sum(m_ws[i,t,]*rets[nn+t,])) ))
  }
}

# > dim(m_TO)
# [1]   6 797

apply(m_TO, 1,mean,na.rm=TRUE)


####----------------
## Concentration
####----------------

# > dim(ws_gmv)
# [1]    6 1000  797   14

CO = array(NA,dim=c(length(models),post.sample,K))
# > dim(TO)
# [1]    6 1000  797

for(i in 1:length(models)){
  for(m in 1:post.sample){
    for(t in 1:K){
      CO[i,m,t] = sqrt(sum(ws_gmv[i,m,t,]^2))
    }
  }
  print(i)
}

apply((apply(CO,c(1,2),mean,na.rm=TRUE)),1,median)


# at the mean

m_CO = array(NA,dim=c(length(models),K))

for(i in 1:length(models)){
  for(t in 1:(K)){
    m_CO[i,t]=sqrt(sum(m_ws[i,t,]^2))
  }
}

apply(m_CO, 1,mean,na.rm=TRUE)

####----------------
## Short Position
####----------------

# > dim(ws_gmv)
# [1]    6 1000  797   14

SP = array(NA,dim=c(length(models),post.sample,K))
# > dim(TO)
# [1]    6 1000  797

for(i in 1:length(models)){
  for(m in 1:post.sample){
    for(t in 1:K){
      SP[i,m,t] = sum((ws_gmv[i,m,t,]<0)* ws_gmv[i,m,t,])
    }
  }
  print(i)
}

apply((apply(SP,c(1,2),mean,na.rm=TRUE)),1,median)


# at the mean

m_SP = array(NA,dim=c(length(models),K))

for(i in 1:length(models)){
  for(t in 1:(K)){
    m_SP[i,t]=sum(m_ws[i,t,]*(m_ws[i,t,]<0))
  }
}

apply(m_SP, 1,mean,na.rm=TRUE)





# 1. jore's1
# 2. geweke's
# 3. equally w.
# 4. caw
# 5. dcct
# 6. dcc-heavy-t

apply(apply(gvm_ret,c(1,3),median),1,mean)
apply(apply(gvm_ret,c(1,3),median),1,sd)

load('data/rf.RData')

Sharpe = (apply(apply(gvm_ret,c(1,3),median),1,mean)*252)/(apply(apply(gvm_ret,c(1,3),median),1,sd)*sqrt(252))
Sharpe
AdjSharpe = (apply(sweep(apply(gvm_ret*252,c(1,3),median),2,rf,'-'),1,mean))/(apply(apply(gvm_ret,c(1,3),median),1,sd)*sqrt(252))
AdjSharpe


all.res.gmv = list()
for(i in 1:length(models)){
  var05     = apply(gvm_ret[i,,],1,quantile,0.05)
  var10     = apply(gvm_ret[i,,],1,quantile,0.10)
  es05      = apply(gvm_ret[i,,],1,esfun,0.05)
  es10      = apply(gvm_ret[i,,],1,esfun,0.10)
  psd       = apply(gvm_ret[i,,],1,sd)
  GL        = (100*(apply(gvm_ret[4,,],1,sd)-apply(gvm_ret[i,,],1,sd))/
                apply(gvm_ret[i,,],1,sd))
  shr       = apply(gvm_ret[i,,-1],1,sum)/(apply(gvm_ret[i,,-1],1,sd)*sqrt(252))
  all.res.gmv   = c(all.res.gmv,list(data.frame(var05,var10,es05,es10,GL,shr,psd)))
}



res.gmv = NULL

for(i in 1:length(models)){
  res.gmv = rbind(res.gmv,c(quantile(all.res.gmv[[i]]$GL,0.05),
                            median(all.res.gmv[[i]]$GL),
                            quantile(all.res.gmv[[i]]$GL,0.95),
                            quantile(all.res.gmv[[i]]$psd*100,0.05),
                            median(all.res.gmv[[i]]$psd*100),
                            quantile(all.res.gmv[[i]]$psd*100,0.95),
                            quantile(all.res.gmv[[i]]$es05,0.05),
                            median(all.res.gmv[[i]]$es05),
                            quantile(all.res.gmv[[i]]$es05,0.95)))
}

res.gmv


res  = res.gmv[c(2,1,3,4,5,6),]

rownames(res ) = models[c(2,1,3,4,5,6)]
colnames(res ) = c('P05','Median','P95','P05','Median','P95')

round(res ,3)

rm(reshf,Sig,Sigma,reslf,resH)
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



