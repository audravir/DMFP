rm(list=ls(all=TRUE))
load('temp/res_10.Rdata')
library(xtable)

# contains information about the marginal distributions
load('temp/RV_forc.Rdata')
load('temp/marginals.Rdata')
marginals$rvs=RV_forc$`1sa`
load('temp/results_scalar_dcc_t.Rdata')
resdcct = res
load('temp/results_xm1.Rdata')
resxm1  = res

rm(res)

MCMCsize=10000

# 1. jore's10
# 2. geweke's
# 3. equally w.
# 4. xm1
# 5. dcct
models = c('Jore1','Geweke','Eq.','AIW','DCC-t')
ws_gmv = array(NA,c(5,length(ind),K,dm))

# for DCC-t model
Q       = array(NA,c(dm, dm, T+K))
R       = array(NA,c(dm, dm, T+K))
a      <- resdcct$restdcc[ind,1]
b      <- resdcct$restdcc[ind,2]
nu     <- resdcct$restdcc[ind,3]
Q[,,1] <- cor(data[start:T,])
R[,,1] <- diag(diag(Q[,,1])^{-1/2})%*%Q[,,1]%*%diag(diag(Q[,,1])^{-1/2})

# for XM1 model
Sbar   = Reduce('+',Sig[start:T])/(T-start+1)
iota   = rep(1,dm)
M      = dim(resxm1$resc)[1]
lag = resxm1$resc[ind,1]
nux = resxm1$resc[ind,2]
b1  = resxm1$resc[ind,3:(dm+2)]
b2  = resxm1$resc[ind,(dm+3):(2*dm+2)]

iota = rep(1,dm)

for(m in 1:length(ind)){
  # for DCC-t model
  tdata  <- qt(udata,nu[m])
  S       = cov(tdata[start:T,])
  
  # xm1
  B0  = (iota%*%t(iota)-(b1[m,]%*%t(b1[m,]))-(b2[m,]%*%t(b2[m,])))*Sbar
  
  for(t in 2:(T+K)){
    # for DCC-t model
    Q[,,t]   <- S*(1-a[m]-b[m])+a[m]*(tdata[t-1,]%*%t(tdata[t-1,]))+b[m]*Q[,,(t-1)]
    R[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
    
    if(t>T){
      Vpred = B0+(b1[m,]%*%t(b1[m,]))*Sig[[t-1]]+
        (b2[m,]%*%t(b2[m,]))*Reduce('+',Sig[(t-lag[m]):(t-1)])/lag[m]
      mtdf = nux[m]-dm
      sample_uxm = pt(mvtnorm::rmvt(MCMCsize,(mtdf-1)/(mtdf+1)*Vpred,df=mtdf+1),df = mtdf+1)
      sample_udcct = pt(mvtnorm::rmvt(MCMCsize,R[,,t],df=nu[m]),df = nu[m])
      
      sample_standretxm = qnorm(sample_uxm)
      sample_standretdcct = qnorm(sample_udcct)
      
      sample_retsxm = ((t(sample_standretxm)*marginals$sd)+marginals$mean)*marginals$rvs[t-T,]
      sample_retsdcct = ((t(sample_standretdcct)*marginals$sd)+marginals$mean)*marginals$rvs[t-T,]
      
      
      
      
      ################## 
      #sample_retsxm=t(sample_standretxm )    ##################
      #sample_retsdcct=t(sample_standretdcct)  ###############
      
      invS = solve(cov(t(sample_retsxm)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den 
      ws_gmv[4,m,t-T,]= pws

      
      invS = solve(cov(t(sample_retsdcct)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den 
      ws_gmv[5,m,t-T,]= pws

      
      if(ws_jore1[m,t-T]>runif(1)) selected = sample_retsxm
      else selected = sample_retsdcct
      
      invS = solve(cov(t(selected)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den 
      ws_gmv[1,m,t-T,]= pws
      
      if(ws_gew[m,t-T]>runif(1)) selected = sample_retsxm
      else selected = sample_retsdcct
    
      invS = solve(cov(t(selected)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den 
      ws_gmv[2,m,t-T,]= pws
      
      if(0.5>runif(1)) selected = sample_retsxm
      else selected = sample_retsdcct
      
      invS = solve(cov(t(selected)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den 
      ws_gmv[3,m,t-T,]= pws
      
    }
    
  }
}



###-----------
### GMV portfolio: mean, var, sharpe, TO,CO, SP, VaR, ES
###-----------
esfun=function(x,p){
  es=mean(x[which(x<quantile(x,p))])
  return(es)
}

gvm_var = gvm_ret = gvm_sharpe = gmv_to = gmv_sp = gmv_co=
  array(NA,dim=c(5,length(ind),K))

for(m in 1:length(ind)){
   for(i in 1:5){
     for(t in 1:K){
       gvm_ret[i,m,t] = sum(ws_gmv[i,m,t,]*rets[T+t,])
       gvm_var[i,m,t] = t(ws_gmv[i,m,t,])%*%RCov[,,T+t]%*%(ws_gmv[i,m,t,])
       gvm_sharpe[i,m,t]=gvm_ret[i,m,t]/sqrt(gvm_var[i,m,t]) 
       gmv_co[i,m,t] = sqrt(sum(ws_gmv[i,m,t,]^2))
       gmv_sp[i,m,t] = sum(ws_gmv[i,m,t,]*(ws_gmv[i,m,t,]<0))
       }
   }
}

for(m in 1:length(ind)){
  for(i in 1:5){
    for(t in 2:K){
      parts = rep(NA,dm)
      for(d in 1:dm){
        parts[d] = abs(ws_gmv[i,m,t,d] - ws_gmv[i,m,t-1,d]*(1+rets[T+t-1,d])/
          (1+sum(ws_gmv[i,m,t-1,]*rets[T+t-1,])))
      }
      gmv_to[i,m,t]  = sum(parts)
    }
  }
}





res = matrix(NA,ncol=5,nrow=10)

for(i in 1:5){
  res[1,i] = median(apply(gvm_ret[i,,],1,quantile,0.05))
  res[2,i] = median(apply(gvm_ret[i,,],1,quantile,0.1))
  res[3,i] = median(apply(gvm_ret[i,,],1,esfun,0.05))
  res[4,i] = median(apply(gvm_ret[i,,],1,esfun,0.1))
  res[5,i] = median(apply(gvm_var[i,,],1,mean))
  res[6,i] = median(apply(gvm_ret[i,,],1,mean))
  res[7,i] = median(apply(gvm_sharpe[i,,],1,mean))
  res[8,i] = median(apply(gmv_to[i,,-1],1,mean))
  res[9,i] = median(apply(gmv_co[i,,],1,mean))
  res[10,i] = median(apply(gmv_sp[i,,],1,mean))
}

res = res[,c(2,1,3,4,5)]
colnames(res) = c( 'Geweke','Jore1','Equal','AIW','DCC-t')
rownames(res) = c('VaR5%','VaR10%','ES5%','ES10%',
                  'var','mean','Sharpe','TO','CO','SP')

print(xtable(res,
             caption = "GMV portfolio results based on 1-step-ahead predicitons 
             for  2009/01/02-2009/12/31 out-of-sample period ($K=252$ observations).
             The table reports the posterior medians of various Global Minimum
             Variance portfolio metrics for the pooled models: 
             Geweke's, Jore's and equally weighted, 
             as well as two best individual models, Additive Inverse Wishart (AIW) and 
             Dynamic Conditional Correlation with $t$ copula (DCC-t).",
             label = 'table:gmvfull3', digits = 3),
      file='tables_and_figures/gmvfull3.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE)


rm(resxm1,Sig,Sigma,resdcct)
save.image('portfolio.Rdata')
