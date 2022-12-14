rm(list=ls(all=TRUE))
load('temp/10variate/res_10.Rdata')
library(xtable)

# contains information about the marginal distributions
load('temp/10variate/RV_forc5.Rdata')
# forc 5 is from log-HAR model (logRV as dep.variable)
load('temp/10variate/marginals.Rdata')

marginals$rvs=RV_forc$`1sa`
load('temp/10variate/results_scalar_dcc_t.Rdata')
resdcct = res
load('temp/10variate/results_xm1.Rdata')
resxm1  = res
load('temp/10variate/results_dcch_t.Rdata')
resdccht  = res

rm(res)

MCMCsize=10000

# 1. jore's10
# 2. geweke's
# 3. equally w.
# 4. xm1
# 5. dcct
# 6. dcc-heavy-t
models = c('Jore1','Geweke','Equal','AIW','DCC-t','DCC-HEAVY-t')
ws_gmv = array(NA,c(length(models),length(ind),K,dm))

# for DCC-t model
Q       = array(NA,c(dm, dm, T+K))
R       = array(NA,c(dm, dm, T+K))
a      <- resdcct$restdcc[ind,1]
b      <- resdcct$restdcc[ind,2]
nu     <- resdcct$restdcc[ind,3]
Q[,,1] <- cor(data[start:T,])
R[,,1] <- diag(diag(Q[,,1])^{-1/2})%*%Q[,,1]%*%diag(diag(Q[,,1])^{-1/2})

# for DCC-HEAVY-t model
Rh      = array(NA,c(dm, dm, T+K))
ah      <- resdccht$restdcch[ind,1]
bh      <- resdccht$restdcch[ind,2]
nuh     <- resdccht$restdcch[ind,3]
Pbar = Reduce('+',Sig[1:T])/T
Rbar = cor(data[start:T,])
Rh[,,1] <- Rbar

# for XM1 model
Sbar   = Reduce('+',Sig[start:T])/(T-start+1)
iota   = rep(1,dm)
M      = dim(resxm1$resc)[1]
lag = resxm1$resc[ind,1]
nux = resxm1$resc[ind,2]
b1  = resxm1$resc[ind,3:(dm+2)]
b2  = resxm1$resc[ind,(dm+3):(2*dm+2)]

for(m in 1:length(ind)){
  # for DCC-t model
  tdata  <- qt(udata,nu[m])
  S       = cov(tdata[start:T,])
  
  # for DCC-HEAVY-t model
  tdatah  <- qt(udata,nuh[m])
  Rbar   <- cor(tdatah[start:T,])
  
  # xm1
  B0  = (iota%*%t(iota)-(b1[m,]%*%t(b1[m,]))-(b2[m,]%*%t(b2[m,])))*Sbar
  
  for(t in 2:(T+K)){
    # for DCC-t model
    Q[,,t]   <- S*(1-a[m]-b[m])+a[m]*(tdata[t-1,]%*%t(tdata[t-1,]))+b[m]*Q[,,(t-1)]
    R[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
    
    # for DCC-HEAVY-t model
    Rh[,,t]   <- Rbar+ah[m]*(Sig[[t-1]]-Pbar)+bh[m]*(Rh[,,t-1]-Rbar)
    
    if(t>T){
      Vpred = B0+(b1[m,]%*%t(b1[m,]))*Sig[[t-1]]+
        (b2[m,]%*%t(b2[m,]))*Reduce('+',Sig[(t-lag[m]):(t-1)])/lag[m]
      mtdf = nux[m]-dm
      sample_uxm = pt(mvtnorm::rmvt(MCMCsize,(mtdf-1)/(mtdf+1)*Vpred,df=mtdf+1),df = mtdf+1)
      sample_udcct = pt(mvtnorm::rmvt(MCMCsize,R[,,t],df=nu[m]),df = nu[m])
      sample_udccth = pt(mvtnorm::rmvt(MCMCsize,Rh[,,t],df=nuh[m]),df = nuh[m])
      
      sample_standretxm = qnorm(sample_uxm)
      sample_standretdcct = qnorm(sample_udcct)
      sample_standretdccth = qnorm(sample_udccth)
      
      sample_retsxm = ((t(sample_standretxm)*marginals$sd)+marginals$mean)*marginals$rvs[t-T,]
      sample_retsdcct = ((t(sample_standretdcct)*marginals$sd)+marginals$mean)*marginals$rvs[t-T,]
      sample_retsdccth = ((t(sample_standretdccth)*marginals$sd)+marginals$mean)*marginals$rvs[t-T,]
      
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

      invS = solve(cov(t(sample_retsdccth)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den 
      ws_gmv[6,m,t-T,]= pws
      
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
  array(NA,dim=c(length(models),length(ind),K))

for(m in 1:length(ind)){
   for(i in 1:length(models)){
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
  for(i in 1:length(models)){
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


resm =matrix(NA,ncol=length(models),nrow=12)
GL   = matrix(NA,ncol=length(models),nrow=1000)

for(i in 1:length(models)){
  resm[1,i] = median(apply(gvm_ret[i,,],1,quantile,0.05))
  resm[2,i] = median(apply(gvm_ret[i,,],1,quantile,0.1))
  resm[3,i] = median(apply(gvm_ret[i,,],1,esfun,0.05))
  resm[4,i] = median(apply(gvm_ret[i,,],1,esfun,0.1))
  resm[5,i] = median(apply(gmv_to[i,,-1],1,mean))
  resm[6,i] = median(apply(gmv_co[i,,],1,mean))
  resm[7,i] = median(apply(gmv_sp[i,,],1,mean))
  
  resm[8,i] = median(apply(gvm_ret[i,,],1,sum))#annualized
  resm[9,i] = median(apply(gvm_ret[i,,],1,sd)*sqrt(252))
  
  resm[10,i] = median(apply(gvm_ret[i,,-1]-0.01/252*gmv_to[i,,-1],1,sum))
  
  resm[11,i] = median(apply(gvm_ret[i,,-1]-0.01/252*gmv_to[i,,-1],1,sum)/
                        (apply(gvm_ret[i,,-1],1,sd)*sqrt(252)))
  
  resm[12,i] = median(100*sqrt(252)*(apply(gvm_ret[4,,],1,sd)-apply(gvm_ret[i,,],1,sd))/
                        apply(gvm_ret[i,,],1,sd))
  GL[,i]   = (100*sqrt(252)*(apply(gvm_ret[4,,],1,sd)-apply(gvm_ret[i,,],1,sd))/
                        apply(gvm_ret[i,,],1,sd))
}

res  = resm[,c(2,1,3,4,5,6)]

colnames(res ) = models[c(2,1,3,4,5,6)]
rownames(res ) = c('VaR5%','VaR10%','ES5%','ES10%','TO','CO','SP',
                   'return','sd','adj.return(1%)',
                   'adj.Sharpe','G/L')

round(res ,3)

rm(resxm1,Sig,Sigma,resdcct,resdccht)

save.image('temp/10variate/portfolio.Rdata')

print(xtable(res,align='ccccccc',
             caption = "GMV portfolio results based on 1-step-ahead predicitons 
             for  2009/01/02-2009/12/31 out-of-sample period ($K=252$ observations).
             The table reports the posterior medians of various Global Minimum
             Variance portfolio metrics for the pooled models: 
             Geweke's, Jore's and equally weighted, 
             as well as two best individual models, Additive Inverse Wishart (AIW) and 
             Dynamic Conditional Correlation with $t$ copula (DCC-t).",
             label = 'table:gmvfull', digits = 3),
      file='tables_and_figures/gmvfull.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,
      hline.after = getOption("xtable.hline.after", c(-1,0,7,nrow(resm))))



