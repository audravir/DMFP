rm(list=ls(all=TRUE))
load('data/10data.Rdata')
library(xtable)
library(zoo)
library(truncnorm)
library(DMPF)


data  = stand # Prediction etc is performed on STANDARDIZED returns
dm    = dim(data)[2] 
start = 1
T     = 1990
Sig   = Sigma
K     = 252
post.sample = 1000

##------
## Static
#  did not estimate separately
#  estimator is here, with some IWishart prior for the R matrix
#  the posterior predictive is multivariate t
##------

lL_static = matrix(NA,ncol=K,nrow=1)

for(t in (T+1):(T+K)){
  stdf = 10+(t-1)-dm+1
  scm  = ((10-dm-1)*diag(dm)+t(data[start:(t-1),])%*%data[start:(t-1),])/stdf 
  lL_static[(t-T)] <- mvtnorm::dmvt(data[t,],delta=rep(0,dm),scm,df = stdf,log=FALSE)
}

sum(log(lL_static)) #the LPS for K=242


##------
## Rmf
#  Riskmetrics fixed
#  no estimation here, all is empirical
##------

Q       = array(NA,c(dm, dm, T+K))
R       = array(NA,c(dm, dm, T+K))
Q[,,1] <- cor(data[start:T,])
R[,,1] <- diag(diag(Q[,,1])^{-1/2})%*%Q[,,1]%*%diag(diag(Q[,,1])^{-1/2})
lamRMf  = 0.999
Vpred  <- list()
lL_rmf = matrix(NA,ncol=K,nrow=1)

for(t in 2:(T+K)){
  Q[,,t]   <- (1-lamRMf)*(data[t-1,]%*%t(data[t-1,]))+lamRMf*Q[,,(t-1)]
  R[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
  if(t>T){
    lL_rmf[(t-T)] <- mvtnorm::dmvnorm(data[t,], rep(0,dm), R[,,t], log=F)
  }}

sum(log(lL_static))
sum(log(lL_rmf))

##------
## scalar dcc Gaussian Copula
##------

load('temp/10variate/results_scalar_dcc.Rdata')
M   = dim(res$resdcc)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

Q       = array(NA,c(dm, dm, T+K))
R       = array(NA,c(dm, dm, T+K))
Q[,,1] <- cor(data[start:T,])
R[,,1] <- diag(diag(Q[,,1])^{-1/2})%*%Q[,,1]%*%diag(diag(Q[,,1])^{-1/2})
a      <- res$resdcc[ind,1]
b      <- res$resdcc[ind,2]
Vpred  <- list()
lL_dcc = matrix(NA,ncol=K,nrow=length(ind))

for(m in 1:length(ind)){
  for(t in 2:(T+K)){
    Q[,,t]   <- cov(data[start:T,])*(1-a[m]-b[m])+a[m]*(data[t-1,]%*%t(data[t-1,]))+b[m]*Q[,,(t-1)]
    R[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
    if(t>T){
      lL_dcc[m,(t-T)] <- mvtnorm::dmvnorm(data[t,], rep(0,dm), R[,,t], log=F)
    }
  }}

sum(log(lL_static))
sum(log(lL_rmf))
sum(apply(log(lL_dcc),2,mean))



postcorrsdcc = matrix(NA,ncol=dm*(dm-1)/2,nrow = K+T)

for(t in 2:(T+K)){
  tmp = 0
  for(m in 1:length(ind)){
    Q[,,t]   <- cov(data[start:T,])*(1-a[m]-b[m])+a[m]*(data[t-1,]%*%t(data[t-1,]))+b[m]*Q[,,(t-1)]
    R[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
    tmp = tmp+ R[,,t][lower.tri(R[,,t])]
  }
  postcorrsdcc[t,] <- tmp/length(ind)
}

##------
## Rme
##------

load('temp/10variate/results_RMe.Rdata')
M   = length(res$resRMe)
ind = round(seq(1,M,length=post.sample)) #thin every xth

Q       = array(NA,c(dm, dm, T+K))
R       = array(NA,c(dm, dm, T+K))
Q[,,1] <- cor(data[start:T,])
R[,,1] <- diag(diag(Q[,,1])^{-1/2})%*%Q[,,1]%*%diag(diag(Q[,,1])^{-1/2})
lam = res$resRMe[ind]
Vpred  <- list()
lL_rme = matrix(NA,ncol=K,nrow=length(ind))

for(m in 1:length(ind)){
  for(t in 2:(T+K)){
    Q[,,t]   <- (1-lam[m])*(data[t-1,]%*%t(data[t-1,]))+lam[m]*Q[,,(t-1)]
    R[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
    if(t>T){
      lL_rme[m,(t-T)] <- mvtnorm::dmvnorm(data[t,], rep(0,dm), R[,,t], log=F)
    }
  }}

sum(log(lL_static))
sum(log(lL_rmf))
sum(apply(log(lL_dcc),2,mean))
sum(apply(log(lL_rme),2,mean))

##-----------------------
## scalar dcc t Copula
##-----------------------

load('temp/10variate/results_scalar_dcc_t.Rdata')
M   = dim(res$restdcc)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth

Q       = array(NA,c(dm, dm, T+K))
R       = array(NA,c(dm, dm, T+K))
a      <- res$restdcc[ind,1]
b      <- res$restdcc[ind,2]
nu     <- res$restdcc[ind,3]

Q[,,1] <- cor(data[start:T,])
R[,,1] <- diag(diag(Q[,,1])^{-1/2})%*%Q[,,1]%*%diag(diag(Q[,,1])^{-1/2})
lL_tdcc = matrix(NA,ncol=K,nrow=length(ind))

for(m in 1:length(ind)){
  tdata  <- qt(udata,nu[m])
  S       = cov(tdata[start:T,])
  
  for(t in 2:(T+K)){
    Q[,,t]   <- S*(1-a[m]-b[m])+a[m]*(tdata[t-1,]%*%t(tdata[t-1,]))+b[m]*Q[,,(t-1)]
    R[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
    if(t>T){
      lL_tdcc[m,(t-T)] <- mvtnorm::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nu[m], log=F)*
        prod(dnorm(data[t,])/dt(tdata[t,],df=nu[m]))
    }
}}


sum(log(lL_static))
sum(log(lL_rmf))
sum(apply(log(lL_dcc),2,mean))
sum(apply(log(lL_rme),2,mean))
sum(apply(log(lL_tdcc),2,mean))


##-----------------------
## dcc-HEAVY t Copula
##-----------------------

load('temp/10variate/results_dcch_t.Rdata')
M   = dim(res$restdcch)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth

R       = array(NA,c(dm, dm, T+K))
a      <- res$restdcch[ind,1]
b      <- res$restdcch[ind,2]
nu     <- res$restdcch[ind,3]

Pbar = Reduce('+',Sig[1:T])/T
Rbar = cor(data[start:T,])

R[,,1] <- Rbar

lL_tdcch = matrix(NA,ncol=K,nrow=length(ind))

for(m in 1:length(ind)){
  tdata  <- qt(udata,nu[m])
  Rbar   <- cor(tdata[start:T,])

  for(t in 2:(T+K)){
    R[,,t]   <- Rbar+a[m]*(Sig[[t-1]]-Pbar)+b[m]*(R[,,t-1]-Rbar)
    
    if(t>T){
      lL_tdcch[m,(t-T)] <- mvtnorm::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nu[m], log=F)*
        prod(dnorm(data[t,])/dt(tdata[t,],df=nu[m]))
    }
  }
}

sum(log(lL_static))
sum(log(lL_rmf))
sum(apply(log(lL_dcc),2,mean))
sum(apply(log(lL_rme),2,mean))
sum(apply(log(lL_tdcc),2,mean))
sum(apply(log(lL_tdcch),2,mean))


##------
## xm1
##------

Sbar   = Reduce('+',Sig[1:T])/T
iota   = rep(1,dm)

load('temp/10variate/results_xm1.Rdata')
M   = dim(res$resc)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth
lL_xm1 = matrix(NA,ncol=K,nrow=length(ind))
lag = res$resc[ind,1]
nu  = res$resc[ind,2]
b1  = res$resc[ind,3:(dm+2)]
b2  = res$resc[ind,(dm+3):(2*dm+2)]

for(m in 1:length(ind)){
    Vpred = list()
    B0  = (iota%*%t(iota)-(b1[m,]%*%t(b1[m,]))-(b2[m,]%*%t(b2[m,])))*Sbar
    for(t0 in 1:K){
      Vpred[[t0]] = B0+(b1[m,]%*%t(b1[m,]))*Sig[[T+t0-1]]+
        (b2[m,]%*%t(b2[m,]))*Reduce('+',Sig[(T+t0-lag[m]):(T+t0-1)])/lag[m]
      mtdf = nu[m]-dm
      lL_xm1[m,t0] = mvtnorm::dmvt(data[(T+t0),],delta=rep(0,dm),(mtdf-1)/(mtdf+1)*Vpred[[t0]],df = mtdf+1,log=FALSE)
}}


sum(log(lL_static))
sum(log(lL_rmf))
sum(apply(log(lL_dcc),2,mean))
sum(apply(log(lL_rme),2,mean))
sum(apply(log(lL_tdcc),2,mean))
sum(apply(log(lL_tdcch),2,mean))
sum(apply(log(lL_xm1),2,mean))

##------
## dcc-HEAVY Gaussian Copula
##------

load('temp/10variate/results_dcch.Rdata')
M   = dim(res$resdcch)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

R       = array(NA,c(dm, dm, T+K))
R[,,1] <- cor(data[start:T,])
a      <- res$resdcch[ind,1]
b      <- res$resdcch[ind,2]
Rbar = cor(data[start:T,])
Pbar = Reduce('+',Sigma[1:T])/T
lL_dcch = matrix(NA,ncol=K,nrow=length(ind))

for(m in 1:length(ind)){
  for(t in 2:(T+K)){
    R[,,t] <- Rbar+a[m]*(Sigma[[t-1]]-Pbar)+b[m]*(R[,,t-1]-Rbar)
    if(t>T){
      lL_dcch[m,(t-T)] <- mvtnorm::dmvnorm(data[t,], rep(0,dm), R[,,t], log=F)
    }
  }
}

sum(log(lL_static))
sum(log(lL_rmf))
sum(apply(log(lL_dcc),2,mean))
sum(apply(log(lL_rme),2,mean))
sum(apply(log(lL_tdcc),2,mean))
sum(apply(log(lL_tdcch),2,mean))
sum(apply(log(lL_xm1),2,mean))
sum(apply(log(lL_dcch),2,mean))


postcorrsh = matrix(NA,ncol=dm*(dm-1)/2,nrow = K+T)

for(t in 2:(T+K)){
  tmp = 0
  for(m in 1:length(ind)){
    R[,,t] <- Rbar+a[m]*(Sigma[[t-1]]-Pbar)+b[m]*(R[,,t-1]-Rbar)
    tmp = tmp+ R[,,t][lower.tri(R[,,t])]
  }
    postcorrsh[t,] <- tmp/length(ind)
}

Z = matrix(NA,ncol=dm*(dm-1)/2,nrow = K+T)
for(t in 1:(K+T)) Z[t,] = Sigma[[t]][upper.tri(Sigma[[t]])]

# rolling window corrs
C = combn(1:dm,2)
Cp = C[,order(C[2,])]
rollcorr = Zroll = matrix(NA,ncol=dm*(dm-1)/2,nrow = K+T)

for(i in 1:(dm*(dm-1)/2)) {
  rollcorr[,i] = rollapply(data[,Cp[,i]],width=252,function(x) cor(x[,1],x[,2]), by.column=FALSE,align = c("center"),fill=NA) 
  Zroll[,i] = rollapply(Z[,i],width=5,mean,fill=NA,align='center')
}
  
pdf('tables_and_figures/comp_heavy.pdf',width=20,height=20)
par(mfrow=c(9,5))
for(i in 1:(dm*(dm-1)/2)){
  plot(date,Zroll[,i],type='l',col='gray70',main=paste(Cp[1,i],',',Cp[2,i],sep=''),ylab='',
       xlab='',ylim=c(-0.3,1),lwd=3)
  lines(date,postcorrsh[,i],col=2,lwd=2)
  lines(date,postcorrsdcc[,i],col=4,lwd=2)
  lines(date,rollcorr[,i],col=1,lwd=2)
  # legend(x=date[1],y=1,col=c('gray70',1,2,4),
  #        lwd=c(2,2,2,2),legend=c('6m-rolling RCor',
  #         '6m-rolling empirical corr','DCC-HEAVY','DCC'))
}
dev.off()


##################

dcc  =cumsum(apply(log(lL_dcc),2,mean))
dcct  =cumsum(apply(log(lL_tdcc),2,mean))
xm   = cumsum(apply(log(lL_xm1),2,mean))


standvol = matrix(NA,ncol=dm,nrow = K)

for(d in 1:dm){
  standvol[,d] = ((RCov[d,d,(T+1):(T+K)])-
                    mean((RCov[d,d,(T+1):(T+K)])))/
    sd((RCov[d,d,(T+1):(T+K)]))
}

mkvol    = rep(NA,K)
for(t0 in 1:K){
  mkvol[t0] = mean(standvol[t0,])
}


##------
## All models separately
##------

pdf('tables_and_figures/all_bfs.pdf',height=5,width=10)
par(mfrow=c(1,1), mar=c(3, 3, 1, 1) + 0.1)
plot(tail(date,K),mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='', xaxt="n")

par(new = TRUE)
plot(tail(date,K),cumsum(apply(log(lL_xm1[,]),2,median))-
       cumsum(log(lL_static[,])), xaxt="n",
     type='l',ylim=c(-5,30),lty=2,ylab='',xlab='',lwd=2)
abline(h=0)
lines(tail(date,K),cumsum(apply(log(lL_tdcc[,]),2,median))-
        cumsum(log(lL_static[,])),lty=2)
lines(tail(date,K),cumsum(apply(log(lL_rme[,]),2,median))-
        cumsum(log(lL_static[,])),col='gray40',lwd=2)
lines(tail(date,K),cumsum(log(lL_rmf[,]))-
        cumsum(log(lL_static[,])),col='gray60',lwd=2,lty=4)
lines(tail(date,K),cumsum(apply(log(lL_dcc[,]),2,median))-
        cumsum(log(lL_static[,])))
lines(tail(date,K),cumsum(apply(log(lL_tdcch),2,median))-
        cumsum(log(lL_static[,])),col='gray40',lty=6,lwd=3)
legend(x=date[(T+1)]-10,y=30,col=c('gray60','gray40',1,1,1,'gray40'),
       lty=c(4,1,1,2,2,6),lwd=c(2,2,1,1,2,2),
       legend=c('RMf','RMe','DCC','DCC-t','AIW','DCC-HEAVY-t'))
atx <- seq(date[(T+1)], date[(T+K)], by=30)
axis(1, at=atx, labels=format(atx, "%Y/%m"))
dev.off()

lpbf = c(sum(log(lL_static[,])),
sum(log(lL_rmf[,])),
sum(apply(log(lL_rme[,]),2,median)),
sum(apply(log(lL_dcc[,]),2,median)),
sum(apply(log(lL_tdcc[,]),2,median)),
sum(apply(log(lL_xm1[,]),2,median)),
sum(apply(log(lL_tdcch),2,median)))
names(lpbf) = c('Static','RMf','RMe','DCC','DCC-t','AIW','DCC-HEAVY-t')
lpbfdf = as.data.frame(t(lpbf))

max(lpbf)

print(xtable(lpbfdf,align= 'ccccccc|c',
             caption = '1-step-ahead log predictive scores ($LPS$) 
             for all individual models: Static, RiskMetrics fixed (RMf),
RiskMetrics estimated (RMe), 
Dynamic conditional correlation with Gaussian and $t$ copulas (DCC
and DCC-t), Additive Inverse Wishart (AIW) and DCC-HEAVY model with $t$ copula for 
2009/01/02-2009/12/31 out-of-sample period
(K = 252 observations).',
             label = 'table:lps', digits = 2),
      file='tables_and_figures/lps.tex',
      include.rownames = FALSE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,scalebox=1)



##------
## Density combination
##------

ws_gew    = lL_gew = matrix(NA,ncol=K,nrow=length(ind))
ws_jore1  = lL_jore1 = matrix(NA,ncol=K,nrow=length(ind))
ws_jore5  = lL_jore5 = matrix(NA,ncol=K,nrow=length(ind))
ws_jore10 = lL_jore10 = matrix(NA,ncol=K,nrow=length(ind))
ws_jore25 = lL_jore25 = matrix(NA,ncol=K,nrow=length(ind))
lL_equal  = lL_DN = matrix(NA,ncol=K,nrow=length(ind))


resDN= PMCMC_delNegro(lL_xm1,lL_tdcc,1000,c(0,2),10000,propsd = 0.5)

mean(resDN$acc)
hist(resDN$beta,freq=FALSE)
grid = seq(-1,1,length=500)
lines(grid,dtruncnorm(grid,-1,1,0,1),lwd=2)
inddn = seq(1,length(resDN$beta),length=1000)
bp = resDN$beta[inddn]
lines(density(bp),col=2,lwd=3)
c(quantile(resDN$beta,0.025),median(resDN$beta),quantile(resDN$beta,0.975))

median(bp)
pacf(bp)

ws_DN = t(pnorm(resDN$weights_xs))


for(m in 1:length(ind)){
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


pdf('tables_and_figures/weights.pdf',height=8,width=10)
par(mfrow=c(2,1), mar=c(3, 3, 1, 1) + 0.1)
plot(date[(T+1):(T+K)],mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(T+1)]-60,date[(T+K)]))
par(new = TRUE)
plot(date[(T+1):(T+K)],apply(ws_gew,2,median),ylim=c(0,1),
     type='l',ylab='',xlab='',xaxt="n",lwd=2,xlim=c(date[(T+1)]-60,date[(T+K)]))
lines(date[(T+1):(T+K)],apply(ws_jore1,2,median),col='gray40',lwd=2,lty=2)
lines(date[(T+1):(T+K)],apply(ws_jore5,2,median),lty=3)
lines(date[(T+1):(T+K)],apply(ws_jore10,2,median),col='gray60',lwd=2)
lines(date[(T+1):(T+K)],apply(ws_DN,2,median),lty=4,lwd=2)
abline(h=0.5)
axis(1, at=atx, labels=format(atx, "%Y/%m"))
legend(x=date[(T+1)]-65,y=1,col=c(1,'gray40',1,'gray60',1),
       lty=c(1,2,3,1,4),lwd=c(2,2,1,2,2),
       legend=c('Geweke','Jore1','Jore5','Jore10','DelNegro'))

plot(date[(T+1):(T+K)],mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(T+1)]-60,date[(T+K)]))
par(new = TRUE)

plot(date[(T+1):(T+K)],cumsum(apply(log(lL_xm1[,]),2,median))-
       cumsum(apply(log(lL_xm1[,]),2,median)),
     type='l',ylim=c(-20,15),ylab='',xlab='',xaxt="n",
     xlim=c(date[(T+1)]-60,date[(T+K)]))
lines(date[(T+1):(T+K)],cumsum(apply(log(lL_tdcc[,]),2,median))-
        cumsum(apply(log(lL_xm1[,]),2,median)))
lines(date[(T+1):(T+K)],cumsum(apply(lL_gew[,],2,median))-
        cumsum(apply(log(lL_xm1[,]),2,median)),
      col=1,lwd=2)
lines(date[(T+1):(T+K)],cumsum(apply(lL_jore1[,],2,median))-
        cumsum(apply(log(lL_xm1[,]),2,median)),
      col='gray40',lwd=2,lty=2)
lines(date[(T+1):(T+K)],cumsum(apply(lL_jore5[,],2,median))-
        cumsum(apply(log(lL_xm1[,]),2,median)),col='gray40',lwd=2,lty=3)
lines(date[(T+1):(T+K)],cumsum(apply(lL_jore10[,],2,median))-
        cumsum(apply(log(lL_xm1[,]),2,median)),col='gray60',lwd=2)
lines(date[(T+1):(T+K)],cumsum(apply(lL_DN,2,median))-
        cumsum(apply(log(lL_xm1[,]),2,median)),lty=4,lwd=2)
lines(date[(T+1):(T+K)],cumsum(apply((lL_equal),2,median))-
        cumsum(apply(log(lL_xm1[,]),2,median)),col='gray60',lwd=2,lty=5)
lines(tail(date,K),cumsum(apply(log(lL_tdcch),2,median))-
        cumsum(apply(log(lL_xm1[,]),2,median)),col='gray40',lty=6,lwd=2)

axis(1, at=atx, labels=format(atx, "%Y/%m"))
legend(x=date[(T+1)]-65,y=15,col=c(1,1,'gray40',1,'gray60',1,'gray60','gray40'),
       lty=c(1,1,2,3,1,4,5,6),lwd=c(1,2,2,1,2,2,2,2),
       legend=c('DCC-t','Geweke','Jore1','Jore5',
                'Jore10','DelNegro','Equal','DCC-HEAVY-t'))
dev.off()

save.image('temp/10variate/res_10.Rdata')

# -----------------------------
# Jore1 vs Dn
# -----------------------------

library(DMFP)

plot(resDN$beta,type='l')
hist(resDN$beta,freq = FALSE)


pdf('tables_and_figures/dn_vs_jore.pdf',height=7,width=10)
par(mfrow=c(2,1))
plot(tail(date,K),apply(ws_jore1,2,median),col=2,type='l',lwd=2,
     main='HF component weight, DN in black, Jore1 in red',ylab='',xlab='')
lines(tail(date,K),apply(ws_DN,2,median),lwd=2)
abline(h=c(0,0.5,1),lty=3)
#legend(x=date[T+1],y=0.8,col=c(1,2),lwd=c(2,2),legend = c('DelNegro','Jore1'))

plot(tail(date,K),cumsum(apply(lL_jore1,2,median))-cumsum(apply(log(lL_xm1),2,median)),
     col=2,type='l',lwd=2,main='BF against AIW model, DN in black, Jore1 in red',ylab='',xlab='')
lines(tail(date,K),cumsum(apply(lL_DN,2,median))-cumsum(apply(log(lL_xm1),2,median)),
      lwd=2)
lines(tail(date,K),cumsum(apply(lL_jore1,2,median))-cumsum(apply(lL_DN,2,median)),
      col='gray60',lty=2,lwd=2)
abline(h=0)
dev.off()






# -----------------------------
# Correlation between weights and avrg volatility
# -----------------------------
VIX <- read.csv("data/VIX.csv")
vix = VIX$Adj.Close

# avrg is the avrg squared returns across 10 assets
# avrg1 is the avrg RV across 10 assets

# names 

library(quantmod)

# live mcap data

mcap = rep(NA,dm)
for(i in 1:dm){
  mcap[i] = getQuote(names[i], what = yahooQF(c("Market Capitalization")))[[2]]/10^9
}

# mcap = c(353.613,468.6,133.403,1945,270.604,6.666,133.525,
#          43.275,117.853,238.393) #June 2021

ws = rep(1/dm,dm)
ws3 = mcap/sum(mcap)

avrg2 = avrg3 = rep(NA,K)
for(t in 1:K){
  avrg2[t] = t(ws)%*%RCov[,,T+t]%*%(ws)
  avrg3[t] = t(ws3)%*%RCov[,,T+t]%*%(ws3)
}

marketvols = cbind(mkvol,avrg2,avrg3,vix)
cor(marketvols)
corrs      = array(NA,dim=c(dim(ws_gew)[1],6,4))

for(m in 1:length(ind)){
  for(i in 1:4){  
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
colnames(corrs_res) = c('avrg RV','Mkt:eql','Mkt:Mcap','VIX')
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
             label = 'table:corrs', digits = 3),
      file='tables_and_figures/corrs.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE)


