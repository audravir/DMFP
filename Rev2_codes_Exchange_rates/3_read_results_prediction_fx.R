rm(list=ls(all=TRUE))
load('data/3data_EX.Rdata')
library(xtable)
library(truncnorm)

data  = stand # Prediction etc is performed on STANDARDIZED returns
dm    = dim(data)[2] 
start = 1

nn.all       = length(date)
end.date = which(zoo::as.yearmon(date)=='Jan 2020')[1]-1
date[end.date]

nn    = end.date
Sig   = Sigma
K     = nn.all-end.date
c(nn,K,nn+K)
post.sample = 1000

##------
## Static
#  did not estimate separately
#  estimator is here, with some IWishart prior for the R matrix
#  the posterior predictive is multivariate t
##------

lL_static = matrix(NA,ncol=K,nrow=1)

for(t in (nn+1):(nn+K)){
  stdf = 10+(t-1)-dm+1
  scm  = ((10-dm-1)*diag(dm)+t(data[start:(t-1),])%*%data[start:(t-1),])/stdf 
  lL_static[(t-nn)] <- mvtnorm::dmvt(data[t,],delta=rep(0,dm),scm,df = stdf,log=FALSE)
}

sum(log(lL_static)) #the LPS for K=242


##------
## Rmf
#  Riskmetrics fixed
#  no estimation here, all is empirical
##------

Q       = array(NA,c(dm, dm, nn+K))
R       = array(NA,c(dm, dm, nn+K))
Q[,,1] <- cor(data[start:nn,])
R[,,1] <- diag(diag(Q[,,1])^{-1/2})%*%Q[,,1]%*%diag(diag(Q[,,1])^{-1/2})
lamRMf  = 0.999
Vpred  <- list()
lL_rmf = matrix(NA,ncol=K,nrow=1)

for(t in 2:(nn+K)){
  Q[,,t]   <- (1-lamRMf)*(data[t-1,]%*%t(data[t-1,]))+lamRMf*Q[,,(t-1)]
  R[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
  if(t>nn){
    lL_rmf[(t-nn)] <- mvtnorm::dmvnorm(data[t,], rep(0,dm), R[,,t], log=F)
  }}

sum(log(lL_static))
sum(log(lL_rmf))

##------
## scalar dcc Gaussian Copula
##------

load('temp/results_scalar_dcc_EX.Rdata')
M   = dim(res$resdcc)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

Q       = array(NA,c(dm, dm, nn+K))
R       = array(NA,c(dm, dm, nn+K))
Q[,,1] <- cor(data[start:nn,])
R[,,1] <- diag(diag(Q[,,1])^{-1/2})%*%Q[,,1]%*%diag(diag(Q[,,1])^{-1/2})
a      <- res$resdcc[ind,1]
b      <- res$resdcc[ind,2]
Vpred  <- list()
lL_dcc = matrix(NA,ncol=K,nrow=post.sample)

for(m in 1:post.sample){
  for(t in 2:(nn+K)){
    Q[,,t]   <- cov(data[start:nn,])*(1-a[m]-b[m])+a[m]*(data[t-1,]%*%t(data[t-1,]))+b[m]*Q[,,(t-1)]
    R[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
    if(t>nn){
      lL_dcc[m,(t-nn)] <- mvtnorm::dmvnorm(data[t,], rep(0,dm), R[,,t], log=F)
    }
  }}

sum(log(lL_static))
sum(log(lL_rmf))
sum(apply(log(lL_dcc),2,mean))

##------
## Rme
##------

load('temp/results_RMe_EX.Rdata')
M   = length(res$resRMe)
ind = round(seq(1,M,length=post.sample)) #thin every xth

Q       = array(NA,c(dm, dm, nn+K))
R       = array(NA,c(dm, dm, nn+K))
Q[,,1] <- cor(data[start:nn,])
R[,,1] <- diag(diag(Q[,,1])^{-1/2})%*%Q[,,1]%*%diag(diag(Q[,,1])^{-1/2})
lam = res$resRMe[ind]
Vpred  <- list()
lL_rme = matrix(NA,ncol=K,nrow=post.sample)

for(m in 1:post.sample){
  for(t in 2:(nn+K)){
    Q[,,t]   <- (1-lam[m])*(data[t-1,]%*%t(data[t-1,]))+lam[m]*Q[,,(t-1)]
    R[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
    if(t>nn){
      lL_rme[m,(t-nn)] <- mvtnorm::dmvnorm(data[t,], rep(0,dm), R[,,t], log=F)
    }
  }}

sum(log(lL_static))
sum(log(lL_rmf))
sum(apply(log(lL_dcc),2,mean))
sum(apply(log(lL_rme),2,mean))

##-----------------------
## scalar dcc t Copula
##-----------------------

load('temp/results_scalar_dcc_t_EX.Rdata')
M   = dim(res$restdcc)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth

Q       = array(NA,c(dm, dm, nn+K))
R       = array(NA,c(dm, dm, nn+K))
a      <- res$restdcc[ind,1]
b      <- res$restdcc[ind,2]
nu     <- res$restdcc[ind,3]

Q[,,1] <- cor(data[start:nn,])
R[,,1] <- diag(diag(Q[,,1])^{-1/2})%*%Q[,,1]%*%diag(diag(Q[,,1])^{-1/2})
lL_tdcc = matrix(NA,ncol=K,nrow=post.sample)

for(m in 1:post.sample){
  tdata  <- qt(udata,nu[m])
  S       = cov(tdata[start:nn,])
  
  for(t in 2:(nn+K)){
    Q[,,t]   <- S*(1-a[m]-b[m])+a[m]*(tdata[t-1,]%*%t(tdata[t-1,]))+b[m]*Q[,,(t-1)]
    R[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
    if(t>nn){
      lL_tdcc[m,(t-nn)] <- mvtnorm::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nu[m], log=F)*
        prod(dnorm(data[t,])/dt(tdata[t,],df=nu[m]))
    }
  }
}


sum(log(lL_static))
sum(log(lL_rmf))
sum(apply(log(lL_dcc),2,mean))
sum(apply(log(lL_rme),2,mean))
sum(apply(log(lL_tdcc),2,mean))

##------
## xm1
##------

Sbar   = Reduce('+',Sig[1:nn])/nn
iota   = rep(1,dm)

load('temp/results_xm1_EX.Rdata')
M   = dim(res$resc)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth
lL_xm1 = matrix(NA,ncol=K,nrow=post.sample)
lag = res$resc[ind,1]
nu  = res$resc[ind,2]
b1  = res$resc[ind,3:(dm+2)]
b2  = res$resc[ind,(dm+3):(2*dm+2)]

for(m in 1:post.sample){
    Vpred = list()
    B0  = (iota%*%t(iota)-(b1[m,]%*%t(b1[m,]))-(b2[m,]%*%t(b2[m,])))*Sbar
    for(t0 in 1:K){
      Vpred[[t0]] = B0+(b1[m,]%*%t(b1[m,]))*Sig[[nn+t0-1]]+
        (b2[m,]%*%t(b2[m,]))*Reduce('+',Sig[(nn+t0-lag[m]):(nn+t0-1)])/lag[m]
      mtdf = nu[m]-dm
      lL_xm1[m,t0] = mvtnorm::dmvt(data[(nn+t0),],delta=rep(0,dm),(mtdf-1)/(mtdf+1)*Vpred[[t0]],df = mtdf+1,log=FALSE)
}}


sum(log(lL_static))
sum(log(lL_rmf))
sum(apply(log(lL_dcc),2,mean))
sum(apply(log(lL_rme),2,mean))
sum(apply(log(lL_tdcc),2,mean))
sum(apply(log(lL_xm1),2,mean))


dcc  = cumsum(apply(log(lL_dcc),2,mean))
dcct = cumsum(apply(log(lL_tdcc),2,mean))
xm   = cumsum(apply(log(lL_xm1),2,mean))


standvol = matrix(NA,ncol=dm,nrow = K)

for(d in 1:dm){
  standvol[,d] = ((RCov[d,d,(nn+1):(nn+K)])-
                    mean((RCov[d,d,(nn+1):(nn+K)])))/
    sd((RCov[d,d,(nn+1):(nn+K)]))
}

mkvol    = rep(NA,K)
for(t0 in 1:K){
  mkvol[t0] = mean(standvol[t0,])
}




##------
## All models separately
##------

#pdf('tables_and_figures/all_bfs_EX.pdf',height=5,width=10)
par(mfrow=c(1,1), mar=c(3, 3, 1, 1) + 0.1)
plot(date[(nn+1):(nn+K)],mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='', xaxt="n")

par(new = TRUE)
plot(date[(nn+1):(nn+K)],cumsum(apply(log(lL_xm1[,]),2,median))-
       cumsum(log(lL_static[,])), xaxt="n",
     type='l',ylim=c(-1,50),lty=2,ylab='',xlab='',lwd=2)
abline(h=0)
lines(date[(nn+1):(nn+K)],cumsum(apply(log(lL_tdcc[,]),2,median))-
        cumsum(log(lL_static[,])),lty=2)
lines(date[(nn+1):(nn+K)],cumsum(apply(log(lL_rme[,]),2,median))-
        cumsum(log(lL_static[,])),col='gray40',lwd=2)
lines(date[(nn+1):(nn+K)],cumsum(log(lL_rmf[,]))-
        cumsum(log(lL_static[,])),col='gray60',lwd=2,lty=4)
lines(date[(nn+1):(nn+K)],cumsum(apply(log(lL_dcc[,]),2,median))-
        cumsum(log(lL_static[,])))
legend(x=date[(nn+1)]-10,y=30,col=c('gray60','gray40',1,1,1),
       lty=c(4,1,1,2,2),lwd=c(2,2,1,1,2),
       legend=c('RMf','RMe','DCC','DCC-t','AIW'))
atx <- seq(date[(nn+1)], date[(nn+K)], by=30)
axis(1, at=atx, labels=format(atx, "%Y/%m"))
#dev.off()

lpbf = c(sum(log(lL_static[,])),
sum(log(lL_rmf[,])),
sum(apply(log(lL_rme[,]),2,median)),
sum(apply(log(lL_dcc[,]),2,median)),
sum(apply(log(lL_tdcc[,]),2,median)),
sum(apply(log(lL_xm1[,]),2,median)))
names(lpbf) = c('Static','RMf','RMe','DCC','DCC-t','AIW')
lpbfdf = as.data.frame(t(lpbf))

max(lpbf)

print(xtable(lpbfdf,
             caption = '1-step-ahead log predictive scores ($LPS$) 
             for all individual models: Static, RiskMetrics fixed (RMf),
RiskMetrics estimated (RMe), 
Dynamic conditional correlation with Gaussian and t copulas (DCC
and DCC-t) and Additive Inverse Wishart (AIW) for 
2020/01/02 - 2021/12/31 out-of-sample period
(K = 504 observations).',
             label = 'table:lps_EX', digits = 2),
      file='tables_and_figures/lps_EX.tex',
      include.rownames = FALSE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE)



##------
## Density combination
##------

ws_gew    = lL_gew = matrix(NA,ncol=K,nrow=post.sample)
ws_jore1  = lL_jore1 = matrix(NA,ncol=K,nrow=post.sample)
ws_jore5  = lL_jore5 = matrix(NA,ncol=K,nrow=post.sample)
ws_jore10 = lL_jore10 = matrix(NA,ncol=K,nrow=post.sample)
ws_jore25 = lL_jore25 = matrix(NA,ncol=K,nrow=post.sample)
lL_equal  = lL_DN = matrix(NA,ncol=K,nrow=post.sample)

resDN= PMCMC_delNegro(lL_xm1,lL_tdcc,1000)

ws_DN = t(pnorm(resDN$weights_xs))

mean(resDN$acc)




for(m in 1:post.sample){
  for (t in 1:K){
    # weights geweke
    a1p = lL_xm1[m,1:t]
    a2p = lL_tdcc[m,1:t]
    opw = function(x){
      -sum(log(x*a1p+(1-x)*a2p))
    }
    ws_gew[m,t] = optim(0.5, opw, gr = NULL,
                      method = c("L-BFGS-B"),
                      lower = 0, upper = 1, hessian = FALSE)$par
    lL_gew[m,t] = log(ws_gew[m,t]*lL_xm1[m,t]+(1-ws_gew[m,t])*lL_tdcc[m,t])
    
    # DN
    
    lL_DN[m,t] = log(ws_DN[m,t]*lL_xm1[m,t]+(1-ws_DN[m,t])*lL_tdcc[m,t])
    
    # jore 1 past
    ws_jore1[m,t] = a1p[t]/(a1p[t]+a2p[t])
    lL_jore1[m,t] = log(ws_jore1[m,t]*lL_xm1[m,t]+(1-ws_jore1[m,t])*lL_tdcc[m,t])
    
    # jore 5 past
    pst = 5
    ws_jore5[m,t] = exp(sum(log(a1p[max(t-pst+1,1):t])))/
      (exp(sum(log(a1p[max(t-pst+1,1):t])))+exp(sum(log(a2p[max(t-pst+1,1):t]))))
    lL_jore5[m,t] = log(ws_jore5[m,t]*lL_xm1[m,t]+(1-ws_jore5[m,t])*lL_tdcc[m,t])
    
    # jore 10 past
    pst = 10
    ws_jore10[m,t] = exp(sum(log(a1p[max(t-pst+1,1):t])))/
      (exp(sum(log(a1p[max(t-pst+1,1):t])))+exp(sum(log(a2p[max(t-pst+1,1):t]))))
    lL_jore10[m,t] = log(ws_jore10[m,t]*lL_xm1[m,t]+(1-ws_jore10[m,t])*lL_tdcc[m,t])
    
    # jore 25 past
    pst = 25
    ws_jore25[m,t] = exp(sum(log(a1p[max(t-pst+1,1):t])))/
      (exp(sum(log(a1p[max(t-pst+1,1):t])))+exp(sum(log(a2p[max(t-pst+1,1):t]))))
    lL_jore25[m,t] = log(ws_jore25[m,t]*lL_xm1[m,t]+(1-ws_jore25[m,t])*lL_tdcc[m,t])
    
    # equally-weighted
    lL_equal[m,t] = log(0.5*lL_xm1[m,t]+0.5*lL_tdcc[m,t])
  }
}


par(mfrow=c(2,1))
plot(apply(ws_jore1,2,median),col=2,type='l',lwd=2)
lines(apply(ws_DN,2,median),lwd=2)

plot(cumsum(apply(lL_jore1,2,median))-cumsum(apply(log(lL_xm1),2,median)),
     col=2,type='l',lwd=2)
lines(cumsum(apply(lL_DN,2,median))-cumsum(apply(log(lL_xm1),2,median)),
      lwd=2)
lines(cumsum(apply(lL_jore1,2,median))-cumsum(apply(lL_DN,2,median)),col=4)
lines(cumsum(apply(log(resDN$likelihood),1,median))-
        cumsum(apply(log(lL_xm1),2,median)),col=3)



pdf('tables_and_figures/weights_EX.pdf',height=8,width=10)
par(mfrow=c(2,1), mar=c(3, 3, 1, 1) + 0.1)
plot(date[(nn+1):(nn+K)],mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(nn+1)]-30,date[(nn+K)]))
par(new = TRUE)
plot(date[(nn+1):(nn+K)],apply(ws_gew,2,median),ylim=c(0,1),
     type='l',ylab='',xlab='',xaxt="n",lwd=2,xlim=c(date[(nn+1)]-30,date[(nn+K)]))
lines(date[(nn+1):(nn+K)],apply(ws_jore1,2,median),col='gray40',lwd=2,lty=2)
lines(date[(nn+1):(nn+K)],apply(ws_jore5,2,median),lty=3)
lines(date[(nn+1):(nn+K)],apply(ws_jore10,2,median),col='gray60',lwd=2)
lines(date[(nn+1):(nn+K)],apply(ws_DN,2,median),lty=4,lwd=2)
abline(h=0.5)
axis(1, at=atx, labels=format(atx, "%Y/%m"))
legend(x=date[(nn+1)]-40,y=1,col=c(1,'gray40',1,'gray60',1),
       lty=c(1,2,3,1,4),lwd=c(2,2,1,2,2),
       legend=c('Geweke','Jore1','Jore5','Jore10','DelNegro'))

plot(date[(nn+1):(nn+K)],mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(nn+1)]-30,date[(nn+K)]))
par(new = TRUE)

best.model = lL_xm1

plot(date[(nn+1):(nn+K)],cumsum(apply(log(lL_xm1[,]),2,median))-
       cumsum(apply(log(best.model[,]),2,median)),
     type='l',ylim=c(-15,15),ylab='',xlab='',xaxt="n",
     xlim=c(date[(nn+1)]-30,date[(nn+K)]))
lines(date[(nn+1):(nn+K)],cumsum(apply(log(lL_tdcc[,]),2,median))-
        cumsum(apply(log(best.model[,]),2,median)))
lines(date[(nn+1):(nn+K)],cumsum(apply(lL_gew[,],2,median))-
        cumsum(apply(log(best.model[,]),2,median)),
      col=1,lwd=2)
lines(date[(nn+1):(nn+K)],cumsum(apply(lL_jore1[,],2,median))-
        cumsum(apply(log(best.model[,]),2,median)),
      col='gray40',lwd=2,lty=2)
lines(date[(nn+1):(nn+K)],cumsum(apply(lL_jore5[,],2,median))-
        cumsum(apply(log(best.model[,]),2,median)),col='gray40',lwd=2,lty=3)
lines(date[(nn+1):(nn+K)],cumsum(apply(lL_jore10[,],2,median))-
        cumsum(apply(log(best.model[,]),2,median)),col='gray60',lwd=2)
lines(date[(nn+1):(nn+K)],cumsum(apply(lL_DN,2,median))-
        cumsum(apply(log(best.model[,]),2,median)),lty=4,lwd=2)
lines(date[(nn+1):(nn+K)],cumsum(apply((lL_equal),2,median))-
        cumsum(apply(log(best.model[,]),2,median)),col='gray60',lwd=2,lty=5)
axis(1, at=atx, labels=format(atx, "%Y/%m"))
legend(x=date[(nn+1)]-40,y=15,col=c(1,1,'gray40',1,'gray60',1,'gray60'),
       lty=c(1,1,2,3,1,4,5),lwd=c(1,2,2,1,2,2,2),
       legend=c('DCC-t','Geweke','Jore1','Jore5','Jore10','DelNegro','Equal'))
dev.off()

save.image('temp/res_3_EX.Rdata')



lpbf.all = c(lpbf,sum(apply(lL_gew,2,median)),
    sum(apply(lL_jore1,2,median)),
  sum(apply(lL_jore5,2,median)),
  sum(apply(lL_jore10,2,median)),
  sum(apply(lL_DN,2,median)),
  sum(apply(lL_equal,2,median)))-sum(apply(log(best.model[,]),2,median))

names(lpbf.all) = c("Static", "RMf",    "RMe",    "DCC",    "DCC-t",  "AIW",
                    "Gew", "Jore1", "Jore5", "Jore10","DN", "Equal")

lpbf.all

print(xtable(t(data.frame(lpbf.all)),
             caption = 'Differences between the 1-step-ahead log predictive scores ($LPS$) 
             between the best fitting model (AIW) and the rest of models: Static, RiskMetrics fixed (RMf),
RiskMetrics estimated (RMe), 
Dynamic conditional correlation with Gaussian and t copulas (DCC
and DCC-t) and various AIW-DCC-t pools for 
2020/01/02 - 2021/12/31 out-of-sample period
(K = 504 observations).',
             label = 'table:lps_all_EX', digits = 2),
      file='tables_and_figures/lps_all_EX.tex',
      include.rownames = FALSE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE)

# -----------------------------
# Correlation between weights and avrg volatility
# -----------------------------

library(tidyquant)
datayahoo = getSymbols('^VIX', from = date[1],
                       to = tail(date,1)+1,warnings = FALSE,
                       auto.assign = FALSE,periodicity = "daily",
                       return.class = 'xts')

tmp = datayahoo$VIX.Adjusted
ind.date = match(date,index(tmp))
vix = tail(tmp[ind.date],K)

plot(vix)

ws = rep(1/dm,dm)

avrg2 =  rep(NA,K)
for(t in 1:K){
  avrg2[t] = t(ws)%*%RCov[,,nn+t]%*%(ws)
}

marketvols = cbind(mkvol,avrg2,vix)
cor(marketvols)
corrs      = array(NA,dim=c(dim(ws_gew)[1],6,dim(marketvols)[2]))

for(m in 1:length(ind)){
  for(i in 1:dim(marketvols)[2]){  
    corrs[m,1,i] = cor(ws_gew[m,],marketvols[,i])
    corrs[m,2,i] = cor(ws_jore1[m,],marketvols[,i])
    corrs[m,3,i] = cor(ws_jore5[m,],marketvols[,i])
    corrs[m,4,i] = cor(ws_jore10[m,],marketvols[,i])
    corrs[m,5,i] = cor(ws_DN[m,],marketvols[,i])
    corrs[m,6,i] = cor((log(lL_xm1)-log(lL_tdcc))[m,],marketvols[,i])
  }
}



corrs_res = apply(corrs,c(2,3),median)
apply(corrs,c(2,3),quantile,0.025) 
apply(corrs,c(2,3),quantile,0.975) 


corrs_res = data.frame(corrs_res)
colnames(corrs_res) = c('avrg RV','Mkt:eql','VIX')
rownames(corrs_res) = c('Geweke','Jore1','Jore5','Jore10','DelNegro',
                        'diff:logLik')

print(xtable::xtable(t(corrs_res),
                     caption = "Posterior medians of sample correlations 
             between the preference for high-frequency
             model  and four proxies for the   market volatility for 
             2009/01/02-2009/12/31 out-of-sample period ($K=252$ observations).
             The preference for the high-frequency model is measured as a 
             high-frequency 
             component weight in various  pooling schemes 
             as well as the difference between the daily log likelihood (diff:logLik) 
             between the AIW and DCC-t models.
             The proxies for the market volatility are: average standardized 
             realized volatility (avrg RV), 
             equally weighted market portfolio realized volatility (Mkt:eql), 
             MCap weighted market portfolio realized volaltity (Mkt:MCap) 
             and VIX index.",
                     label = 'table:corrs_EX', digits = 3),
      file='tables_and_figures/corrs_EX.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE)

