rm(list=ls(all=TRUE))
load('temp/res_3_EX.Rdata')
library(xtable)
library(fPortfolio)


# contains information about the marginal distributions
# load('temp/RV_forc5.Rdata')
# forc 5 is from log-HAR model (logRV as dep.variable)
load('temp/marginals_EX.Rdata')

marginals$rvs=tail(RVs,K)
load('temp/results_scalar_dcc_t_EX.Rdata')
resdcct = res
load('temp/results_xm1_EX.Rdata')
resxm1  = res
load('temp/results_dcch_t_EX.Rdata')
resdccht  = res

rm(res)

MCMCsize=10000

cvarSpec <- portfolioSpec(
  model = list(type = "CVaR", optimize = "minRisk",
               estimator = "covEstimator", tailRisk = list(),
               params = list(alpha = 0.05)))

# 1. jore's1
# 2. geweke's
# 3. equally w.
# 4. xm1
# 5. dcct
# 6. dcc-heavy-t
models = c('Jore1','Geweke','Equal','AIW','DCC-t','DCC-HEAVY-t')
ws_gmv = ws_cvar = array(NA,c(length(models),length(ind),K,dm))

# for DCC-t model
Q       = array(NA,c(dm, dm, nn+K))
R       = array(NA,c(dm, dm, nn+K))
a      <- resdcct$restdcc[ind,1]
b      <- resdcct$restdcc[ind,2]
nu     <- resdcct$restdcc[ind,3]
Q[,,1] <- cor(data[start:nn,])
R[,,1] <- diag(diag(Q[,,1])^{-1/2})%*%Q[,,1]%*%diag(diag(Q[,,1])^{-1/2})

# for DCC-HEAVY-t model
Rh      = array(NA,c(dm, dm, nn+K))
ah      <- resdccht$restdcch[ind,1]
bh      <- resdccht$restdcch[ind,2]
nuh     <- resdccht$restdcch[ind,3]
Pbar = Reduce('+',Sig[1:nn])/nn
Rbar = cor(data[start:nn,])
Rh[,,1] <- Rbar

# for XM1 model
Sbar   = Reduce('+',Sig[start:nn])/(nn-start+1)
iota   = rep(1,dm)
M      = dim(resxm1$resc)[1]
lag = resxm1$resc[ind,1]
nux = resxm1$resc[ind,2]
b1  = resxm1$resc[ind,3:(dm+2)]
b2  = resxm1$resc[ind,(dm+3):(2*dm+2)]

US  = matrix(runif(3*(nn+K)),ncol=nn+K,nrow=3)

for(m in 1:length(ind)){
  t0 = Sys.time()
  
  # for DCC-t model
  tdata  <- qt(udata,nu[m])
  S       = cov(tdata[start:nn,])
  
  # for DCC-HEAVY-t model
  tdatah  <- qt(udata,nuh[m])
  Rbar   <- cor(tdatah[start:nn,])
  
  # xm1
  B0  = (iota%*%t(iota)-(b1[m,]%*%t(b1[m,]))-(b2[m,]%*%t(b2[m,])))*Sbar
  
  for(t in 2:(nn+K)){
    # for DCC-t model
    Q[,,t]   <- S*(1-a[m]-b[m])+a[m]*(tdata[t-1,]%*%t(tdata[t-1,]))+b[m]*Q[,,(t-1)]
    R[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
    
    # for DCC-HEAVY-t model
    Rh[,,t]   <- Rbar+ah[m]*(Sig[[t-1]]-Pbar)+bh[m]*(Rh[,,t-1]-Rbar)
    
    if(t>nn){
      Vpred = B0+(b1[m,]%*%t(b1[m,]))*Sig[[t-1]]+
        (b2[m,]%*%t(b2[m,]))*Reduce('+',Sig[(t-lag[m]):(t-1)])/lag[m]
      mtdf = nux[m]-dm
      sample_uxm = pt(mvtnorm::rmvt(MCMCsize,(mtdf-1)/(mtdf+1)*Vpred,df=mtdf+1),df = mtdf+1)
      sample_udcct = pt(mvtnorm::rmvt(MCMCsize,R[,,t],df=nu[m]),df = nu[m])
      sample_udccth = pt(mvtnorm::rmvt(MCMCsize,Rh[,,t],df=nuh[m]),df = nuh[m])
      
      sample_standretxm = qnorm(sample_uxm)
      sample_standretdcct = qnorm(sample_udcct)
      sample_standretdccth = qnorm(sample_udccth)
      
      sample_retsxm = ((t(sample_standretxm)*marginals$sd)+marginals$mean)*marginals$rvs[t-nn,]
      sample_retsdcct = ((t(sample_standretdcct)*marginals$sd)+marginals$mean)*marginals$rvs[t-nn,]
      sample_retsdccth = ((t(sample_standretdccth)*marginals$sd)+marginals$mean)*marginals$rvs[t-nn,]
      
      ################## 
      #sample_retsxm=t(sample_standretxm )    ##################
      #sample_retsdcct=t(sample_standretdcct)  ###############
      
      invS = solve(cov(t(sample_retsxm)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den 
      ws_gmv[4,m,t-nn,]= pws
      
      # tmp  = as.timeSeries(t(sample_retsxm))
      # port = minvariancePortfolio(tmp, spec = cvarSpec)
      # ws_cvar[4,m,t-nn,] = port@portfolio@portfolio$weights
      
      invS = solve(cov(t(sample_retsdcct)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den 
      ws_gmv[5,m,t-nn,]= pws
      
      # tmp  = as.timeSeries(t(sample_retsdcct))
      # port = minvariancePortfolio(tmp, spec = cvarSpec)
      # ws_cvar[5,m,t-nn,] = port@portfolio@portfolio$weights
      
      invS = solve(cov(t(sample_retsdccth)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den 
      ws_gmv[6,m,t-nn,]= pws
      
      # tmp  = as.timeSeries(t(sample_retsdccth))
      # port = minvariancePortfolio(tmp, spec = cvarSpec)
      # ws_cvar[6,m,t-nn,] = port@portfolio@portfolio$weights

      # if(ws_jore1[m,t-nn]>runif(1)) selected = sample_retsxm
      if(ws_jore1[m,t-nn]>US[1,t-nn]) selected = sample_retsxm
      else selected = sample_retsdcct
      
      invS = solve(cov(t(selected)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den 
      ws_gmv[1,m,t-nn,]= pws
      
      # tmp  = as.timeSeries(t(selected))
      # port = minvariancePortfolio(tmp, spec = cvarSpec)
      # ws_cvar[1,m,t-nn,] = port@portfolio@portfolio$weights
      
      # if(ws_gew[m,t-nn]>runif(1)) selected = sample_retsxm
      if(ws_gew[m,t-nn]>US[2,t-nn]) selected = sample_retsxm
      else selected = sample_retsdcct
    
      invS = solve(cov(t(selected)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den 
      ws_gmv[2,m,t-nn,]= pws
      
      # tmp  = as.timeSeries(t(selected))
      # port = minvariancePortfolio(tmp, spec = cvarSpec)
      # ws_cvar[2,m,t-nn,] = port@portfolio@portfolio$weights
      # 
      # if(0.5>runif(1)) selected = sample_retsxm
      if(0.5>US[3,t-nn]) selected = sample_retsxm
      else selected = sample_retsdcct
      
      invS = solve(cov(t(selected)))
      nom  = invS%*%iota
      den  = as.vector(t(iota)%*%invS%*%iota)
      pws  = nom/den 
      ws_gmv[3,m,t-nn,]= pws
      
      # tmp  = as.timeSeries(t(selected))
      # port = minvariancePortfolio(tmp, spec = cvarSpec)
      # ws_cvar[3,m,t-nn,] = port@portfolio@portfolio$weights
      # 
    }
    
  }
  print(c(m,Sys.time()-t0))
}

save.image(file = '5varportfolio.Rdata')


###-----------
### GMV portfolio: mean, var, sharpe, TO,CO, SP, VaR, ES
###-----------
esfun=function(x,p){
  es=mean(x[which(x<quantile(x,p))])
  return(es)
}

gvm_var = gvm_ret = gvm_sharpe = gmv_to = gmv_sp = gmv_co=
  array(NA,dim=c(length(models),length(ind),K))

# cvar_var = cvar_ret = cvar_sharpe = cvar_es=
#   array(NA,dim=c(length(models),length(ind),K))
# 
for(m in 1:length(ind)){
   for(i in 1:length(models)){
     for(t in 1:K){
       gvm_ret[i,m,t] = sum(ws_gmv[i,m,t,]*rets[nn+t,])
       # gvm_var[i,m,t] = t(ws_gmv[i,m,t,])%*%RCov[,,nn+t]%*%(ws_gmv[i,m,t,])
       # gvm_sharpe[i,m,t]=gvm_ret[i,m,t]/sqrt(gvm_var[i,m,t]) 
       # gmv_co[i,m,t] = sqrt(sum(ws_gmv[i,m,t,]^2))
       # gmv_sp[i,m,t] = sum(ws_gmv[i,m,t,]*(ws_gmv[i,m,t,]<0))
       
        # cvar_ret[i,m,t] = sum(ws_cvar[i,m,t,]*rets[nn+t,])
       # cvar_var[i,m,t] = t(ws_cvar[i,m,t,])%*%RCov[,,nn+t]%*%(ws_cvar[i,m,t,])
       # cvar_sharpe[i,m,t]=cvar_ret[i,m,t]/sqrt(cvar_var[i,m,t])
       # cvar_es[i,m,t] = esfun(cvar_ret[i,m,t],0.05)
       }
   }
}


# for(m in 1:length(ind)){
#   for(i in 1:length(models)){
#     for(t in 2:K){
#       parts = rep(NA,dm)
#       for(d in 1:dm){
#         parts[d] = abs(ws_gmv[i,m,t,d] - ws_gmv[i,m,t-1,d]*(1+rets[nn+t-1,d])/
#           (1+sum(ws_gmv[i,m,t-1,]*rets[nn+t-1,])))
#       }
#       gmv_to[i,m,t]  = sum(parts)
#     }
#   }
# }


# all.res.cvar = list()
# for(i in 1:length(models)){
#   var05     = apply(cvar_ret[i,,],1,quantile,0.05)
#   var10     = apply(cvar_ret[i,,],1,quantile,0.10)
#   es05      = apply(cvar_ret[i,,],1,esfun,0.05)
#   es10      = apply(cvar_ret[i,,],1,esfun,0.10)
#   GL        = (100*sqrt(252)*(apply(cvar_ret[4,,],1,sd)-apply(cvar_ret[i,,],1,sd))/
#                  apply(cvar_ret[i,,],1,sd))
#   shr       = apply(cvar_ret[i,,-1],1,sum)/(apply(cvar_ret[i,,-1],1,sd)*sqrt(252))
#   all.res.cvar   = c(all.res.cvar,list(data.frame(var05,var10,es05,es10,GL,shr)))
# }


all.res.gmv = list()
for(i in 1:length(models)){
  var05     = apply(gvm_ret[i,,],1,quantile,0.05)
  var10     = apply(gvm_ret[i,,],1,quantile,0.10)
  es05      = apply(gvm_ret[i,,],1,esfun,0.05)
  es10      = apply(gvm_ret[i,,],1,esfun,0.10)
  psd       = apply(gvm_ret[i,,],1,sd)
  GL        = (100*sqrt(252)*(apply(gvm_ret[4,,],1,sd)-apply(gvm_ret[i,,],1,sd))/
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
                            quantile(all.res.gmv[[i]]$psd*100,0.95)))
}

res.gmv


res  = res.gmv[c(2,1,3,4,5,6),]

rownames(res ) = models[c(2,1,3,4,5,6)]
colnames(res ) = c('P05','Median','P95','P05','Median','P95')

round(res ,3)

rm(resxm1,Sig,Sigma,resdcct,resdccht)

save.image('temp/portfolio_EX.Rdata')


tableLines <- print(xtable(res,digits = 3,caption="GMV portfolio results based on 1-step-ahead predictions 
             for  2021/01/04-2021/12/31 out-of-sample period 
             ($K=252$ observations) for 5-variate dataset.
             The table reports the posterior 5, 50 and 95 percentiles of G/L criteria as well as
             portfolio standard deviation (in \\%) for the pooled models (Geweke's, Jore's and 
             equally weighted), 
              two best individual models (Additive Inverse Wishart  and 
             Dynamic Conditional Correlation with $t$ copula) and a competitor model (DCC-HEAVY-t).",
                           align = "lccc|ccc",label='table:gmvfull_EX'), 
                    scalebox=0.8,sanitize.text.function=function(x){x})
multicolumns <- "& \\\\multicolumn{3}{c}{G/L}
                 & \\\\multicolumn{3}{c}{Portfolio stdev.}  \\\\\\\\"
tableLines <- sub ("\\\\toprule\\n", paste0 ("\\\\toprule\n", multicolumns, "\n"), tableLines) ## booktabs = TRUE
tableLines <- sub ("\\\\hline\\n",   paste0 ("\\\\hline\n",   multicolumns, "\n"), tableLines) ## booktabs = FALSE
writeLines (tableLines, con = "tables_and_figures/gmvfull_EX.tex")



