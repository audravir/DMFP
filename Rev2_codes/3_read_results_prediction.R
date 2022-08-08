rm(list=ls(all=TRUE))
load('data/10data.Rdata')
library(xtable)
data  = stand # Prediction etc is performed on STANDARDIZED returns
dm    = dim(data)[2] 
start = 1
T     = 1990
Sig   = Sigma
K     = 252

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

load('temp/results_scalar_dcc.Rdata')
M   = dim(res$resdcc)[1] # size of MCMC
ind = round(seq(1,M,length=M/25)) #thin every 5th

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

##------
## Rme
##------

load('temp/results_RMe.Rdata')
M   = length(res$resRMe)
ind = round(seq(1,M,length=M/25)) #thin every 5th

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

load('temp/results_scalar_dcc_t.Rdata')
M   = dim(res$restdcc)[1]
ind = round(seq(1,M,length=M/25)) #thin every 5th

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

##------
## xm1
##------

Sbar   = Reduce('+',Sig[1:T])/T
iota   = rep(1,dm)

load('temp/results_xm1.Rdata')
M   = dim(res$resc)[1]
ind = round(seq(1,M,length=M/25)) #thin every 5th
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
sum(apply(log(lL_xm1),2,mean))


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
plot(date[(T+1):(T+K)],mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='', xaxt="n")

par(new = TRUE)
plot(date[(T+1):(T+K)],cumsum(apply(log(lL_xm1[,]),2,median))-
       cumsum(log(lL_static[,])), xaxt="n",
     type='l',ylim=c(-15,30),lty=2,ylab='',xlab='',lwd=2)
abline(h=0)
lines(date[(T+1):(T+K)],cumsum(apply(log(lL_tdcc[,]),2,median))-
        cumsum(log(lL_static[,])),lty=2)
lines(date[(T+1):(T+K)],cumsum(apply(log(lL_rme[,]),2,median))-
        cumsum(log(lL_static[,])),col='gray40',lwd=2)
lines(date[(T+1):(T+K)],cumsum(log(lL_rmf[,]))-
        cumsum(log(lL_static[,])),col='gray60',lwd=2,lty=4)
lines(date[(T+1):(T+K)],cumsum(apply(log(lL_dcc[,]),2,median))-
        cumsum(log(lL_static[,])))
legend(x=date[(T+1)]-10,y=30,col=c('gray60','gray40',1,1,1),
       lty=c(4,1,1,2,2),lwd=c(2,2,1,1,2),
       legend=c('RMf','RMe','DCC','DCC-t','AIW'))
atx <- seq(date[(T+1)], date[(T+K)], by=30)
axis(1, at=atx, labels=format(atx, "%Y/%m"))
dev.off()

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
2009/01/02-2009/12/31 out-of-sample period
(K = 252 observations).',
             label = 'table:lps', digits = 2),
      file='tables_and_figures/lps.tex',
      include.rownames = FALSE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE)



##------
## Density combination
##------

ws_gew    = lL_gew = matrix(NA,ncol=K,nrow=length(ind))
ws_jore1  = lL_jore1 = matrix(NA,ncol=K,nrow=length(ind))
ws_jore5  = lL_jore5 = matrix(NA,ncol=K,nrow=length(ind))
ws_jore10 = lL_jore10 = matrix(NA,ncol=K,nrow=length(ind))
ws_jore25 = lL_jore25 = matrix(NA,ncol=K,nrow=length(ind))
lL_equal  = matrix(NA,ncol=K,nrow=length(ind))

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

resdn = PMCMC_delNegro(lL_xm1,lL_tdcc,500)

pdf('tables_and_figures/weights.pdf',height=8,width=10)
par(mfrow=c(2,1), mar=c(3, 3, 1, 1) + 0.1)
plot(date[(T+1):(T+K)],mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(T+1)]-30,date[(T+K)]))
par(new = TRUE)
plot(date[(T+1):(T+K)],apply(ws_gew,2,median),ylim=c(0,1),
     type='l',ylab='',xlab='',xaxt="n",lwd=2,xlim=c(date[(T+1)]-30,date[(T+K)]))
lines(date[(T+1):(T+K)],apply(ws_jore1,2,median),col='gray40',lwd=2,lty=2)
lines(date[(T+1):(T+K)],apply(ws_jore5,2,median),lty=3)
lines(date[(T+1):(T+K)],apply(ws_jore10,2,median),col='gray60',lwd=2)
lines(date[(T+1):(T+K)],apply(pnorm(resdn$weights_xs),1,median),lty=4,lwd=2)
abline(h=0.5)
axis(1, at=atx, labels=format(atx, "%Y/%m"))
legend(x=date[(T+1)]-40,y=1,col=c(1,'gray40',1,'gray60',1),
       lty=c(1,2,3,1,4),lwd=c(2,2,1,2,2),
       legend=c('Geweke','Jore1','Jore5','Jore10','DelNegro'))

plot(date[(T+1):(T+K)],mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(T+1)]-30,date[(T+K)]))
par(new = TRUE)

plot(date[(T+1):(T+K)],cumsum(apply(log(lL_xm1[,]),2,median))-
       cumsum(apply(log(lL_xm1[,]),2,median)),
     type='l',ylim=c(-15,15),ylab='',xlab='',xaxt="n",
     xlim=c(date[(T+1)]-30,date[(T+K)]))
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
lines(date[(T+1):(T+K)],cumsum(apply(log(resdn$likelihood),1,median))-
        cumsum(apply(log(lL_xm1[,]),2,median)),lty=4,lwd=2)
lines(date[(T+1):(T+K)],cumsum(apply((lL_equal),2,median))-
        cumsum(apply(log(lL_xm1[,]),2,median)),col='gray60',lwd=2,lty=5)
axis(1, at=atx, labels=format(atx, "%Y/%m"))
legend(x=date[(T+1)]-40,y=15,col=c(1,1,'gray40',1,'gray60',1,'gray60'),
       lty=c(1,1,2,3,1,4,5),lwd=c(1,2,2,1,2,2,2),
       legend=c('DCC-t','Geweke','Jore1','Jore5','Jore10','DelNegro','Equal'))

dev.off()

save.image('temp/res_10.Rdata')

# -----------------------------
# Correlation between weights and avrg volatility
# -----------------------------
VIX <- read.csv("data/VIX.csv")
vix = VIX$Adj.Close

# avrg is the avrg squared returns across 10 assets
# avrg1 is the avrg RV across 10 assets

ws = rep(1/dm,dm)
mcap = c(353.613,468.6,133.403,1945,270.604,6.666,133.525,
         43.275,117.853,238.393)

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
    corrs[m,5,i] = cor(resdn$weights_xs[,m],marketvols[,i])
    corrs[m,6,i] = cor((log(lL_xm1)-log(lL_tdcc))[m,],marketvols[,i])
  }
}



corrs_res = apply(corrs,c(2,3),median) 
apply(corrs,c(2,3),quantile,0.975) 


corrs_res = data.frame(corrs_res)
colnames(corrs_res) = c('avrg RV','Mkt:eql','Mkt:Mcap','VIX')
rownames(corrs_res) = c('Geweke','Jore1','Jore5','Jore10','DelNegro',
                        'diff:logLik')

print(xtable(t(corrs_res),
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


