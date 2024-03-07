rm(list=ls(all=TRUE))
load('data/FXdata.Rdata')
library(xtable)
# library(truncnorm)
library(Rfast)
library(mvnfast)
library(zoo)

data  = stand # Prediction etc is performed on STANDARDIZED returns
dm    = dim(data)[2] 
start = 1

nn.all       = length(date)
end.date = which(zoo::as.yearmon(date)=="ene 2020")[1]-1
if(is.na(date[end.date])){end.date = which(zoo::as.yearmon(date)=="jan 2020")[1]-1}
date[end.date]

nn    = end.date
Sig   = Sigma
K     = nn.all-end.date
c(nn,K,nn+K)
post.sample = 1000
p1 = 5
p2 = 6

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
  lL_static[(t-nn)] <- mvtnorm::dmvt(data[t,],delta=rep(0,dm),scm,df = stdf,log=TRUE)
}

sum(lL_static) #the LPS for K=797 

##------
## Rmf
#  Riskmetrics fixed
#  no estimation here
##------

Q       = array(NA,c(dm, dm, nn+K))
R.rmf   = array(NA,c(dm, dm, nn+K))
Q[,,1] <- cor(data[start:nn,])
R.rmf[,,1] <- diag(diag(Q[,,1])^{-1/2})%*%Q[,,1]%*%diag(diag(Q[,,1])^{-1/2})
lamRMf  = 0.999
Vpred  <- list()
lL_rmf = matrix(NA,ncol=K,nrow=1)

for(t in 2:(nn+K)){
  Q[,,t]   <- (1-lamRMf)*(data[t-1,]%*%t(data[t-1,]))+lamRMf*Q[,,(t-1)]
  R.rmf[,,t]   <- diag(diag(Q[,,t])^{-1/2})%*%Q[,,t]%*%diag(diag(Q[,,t])^{-1/2})
  if(t>nn){
    lL_rmf[(t-nn)] <- mvtnorm::dmvnorm(data[t,], rep(0,dm), R.rmf[,,t], log=TRUE)
  }
}

sum(lL_static)
sum(lL_rmf)

##------
## vector dcc Gaussian Copula
##------

load('empirical/temp/results_vectordcc.Rdata')
M   = dim(res$resdcc)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

Q       = R = array(NA,c(dm, dm, nn+K))
Q[,,1] <- cor(data[start:nn,])
a      <- res$resdcc[ind,1:dm]
b      <- res$resdcc[ind,(dm+1):(dm*2)]
iota   = rep(1,dm)
Oiota  = Outer(iota,iota)
Sbar   = cov(data[start:nn,])

Vpred  <- list()
lL_dcc = matrix(NA,ncol=K,nrow=post.sample)

for(m in 1:post.sample){
  A      = Outer(a[m,],a[m,])
  B      = Outer(b[m,],b[m,])

  for(t in 2:(nn+K)){
    Q[,,t]   <- (Oiota-A-B)*Sbar+A*Outer(data[t-1,],data[t-1,])+B*Q[,,(t-1)]
    t.ma  = Q[,,t]
    t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
    t.R   = Outer(t.dv,t.dv)*t.ma
    R[,,t] = t.R
    if(t>nn){
      lL_dcc[m,(t-nn)] <- mvnfast::dmvn(data[t,], rep(0,dm), t.R, log=TRUE)
    }
  }
}

# should be the same!
res$Vpred[ind][[post.sample]][1:5]
R[,,nn+1][1:5]

sum(lL_static)
sum(lL_rmf)
sum(apply(lL_dcc,2,median))

# at the median of estimated parameters

A      = Outer(apply(a,2,median),apply(a,2,median))
B      = Outer(apply(b,2,median),apply(b,2,median))
R.dcc  = array(NA,c(dm, dm, nn+K))
tmp.dcc = rep(NA,nn+K)

for(t in 2:(nn+K)){
  Q[,,t]   <- (Oiota-A-B)*Sbar+A*Outer(data[t-1,],data[t-1,])+B*Q[,,(t-1)]
  t.ma  = Q[,,t]
  t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
  t.R   = Outer(t.dv,t.dv)*t.ma
  R.dcc[,,t] = t.R
  tmp.dcc[t] = mvnfast::dmvn(data[t,], rep(0,dm), t.R, log=TRUE)
}

sum(tail(tmp.dcc,K))
sum(apply(lL_dcc,2,median))

##-----------------------
## vector dcc t Copula
##-----------------------

load('empirical/temp/results_vectordcc_tcop.Rdata')

M   = dim(res$r)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth

Q  = R = array(NA,c(dm, dm, nn+K))
a  <- res$r[ind,2:(dm+1)]
b  <- res$r[ind,(dm+2):(2*dm+1)]
nu <- res$r[ind,1]

lL_tdcc = matrix(NA,ncol=K,nrow=post.sample)
iota    = rep(1,dm)
Oiota   = Outer(iota,iota)
INL.N   = dnorm(data,log=TRUE)

for(m in 1:post.sample){
  tdata <- qt(udata,nu[m])
  Sbar  <- cova(tdata[start:nn,])
  A     <- Outer(a[m,],a[m,])
  B     <- Outer(b[m,],b[m,])
  B0    <- (Oiota-A-B)*Sbar
  INL   <- dt(tdata,df=nu[m],log=TRUE)
  Q[,,1] <- cora(tdata[start:nn,])
  
  for(t in 2:(nn+K)){
    Q[,,t] <- B0+A*Outer(tdata[t-1,],tdata[t-1,])+B*Q[,,(t-1)]
    t.ma  = Q[,,t]
    t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
    t.R   = Outer(t.dv,t.dv)*t.ma
    R[,,t] = t.R
    if(t>nn){
      lL_tdcc[m,(t-nn)] = mvnfast::dmvt(tdata[t,], rep(0,dm),t.R, nu[m], log=TRUE)+sum(INL.N[t,])-sum(INL[t,])
    }
  }
}

# should be the same!
res$Vpred[ind][[post.sample]][1:5]
R[,,nn+1][1:5]

zoomin = c(1:100)

par(mfrow=c(1,1))
plot(apply(lL_dcc,2,median)[zoomin],type='l')
lines(apply(lL_tdcc,2,median)[zoomin],col=2)

sum(lL_static)
sum(lL_rmf)
sum(apply(lL_dcc,2,median))
sum(apply(lL_tdcc,2,median))

# at the median of estimated parameters

A      = Outer(apply(a,2,median),apply(a,2,median))
B      = Outer(apply(b,2,median),apply(b,2,median))
tdata <- qt(udata,median(nu))
Sbar  <- cova(tdata[start:nn,])
B0    <- (Oiota-A-B)*Sbar
R.dcct = array(NA,c(dm, dm, nn+K))
tmp.dcct = rep(NA,nn+K)
INL   <- dt(tdata,df=median(nu),log=TRUE)

for(t in 2:(nn+K)){
  Q[,,t] <- B0+A*Outer(tdata[t-1,],tdata[t-1,])+B*Q[,,(t-1)]
  t.ma  = Q[,,t]
  t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
  t.R   = Outer(t.dv,t.dv)*t.ma
  R.dcct[,,t] = t.R
  tmp.dcct[t] = mvnfast::dmvt(tdata[t,], rep(0,dm),t.R, median(nu), log=TRUE)+sum(INL.N[t,])-sum(INL[t,])
}

sum(tail(tmp.dcct,K))
sum(apply(lL_tdcc,2,median))


roll_corr <- rollapply(data = cbind(data[,p1], data[,p2]), width = 100,
                       function(z) cor(z[,1], z[,2]), by.column = FALSE,
                       align = "right",fill=NA)
plot(roll_corr,type='l',ylim=c(-1,1))
lines(R.rmf[p1,p2,],type='l',ylim=c(0,1))
lines(R.dcc[p1,p2,],col=2)
lines(R.dcct[p1,p2,],col=4)

#-----------------------
# dcc-HEAVY-scalar t Copula
#-----------------------

load('empirical/temp/results_heavy_t.Rdata')
M     = dim(res$r)[1]
ind   = round(seq(1,M,length=post.sample)) #thin every xth
lL_ht = matrix(NA,ncol=K,nrow=post.sample)

nu = res$r[ind,1]
a  = res$r[ind,2]
b  = res$r[ind,3]

Rpred = res$Rpred[ind]
Pbar  = Reduce('+',Sigma[1:nn])/nn
R     = array(NA,c(dm, dm, nn+K))
INL.N = dnorm(data,log=TRUE)

for(m in 1:post.sample){
  
  tdata  = qt(udata,nu[m])
  Rbar   = cor(tdata[1:nn,])
  R[,,1] = Rbar
  INL.T  = dt(tdata,df=nu[m],log=TRUE)
  
  for(t in 2:(nn+K)){
    R[,,t] = (1-b[m])*Rbar-a[m]*Pbar+a[m]*Sigma[[t-1]]+b[m]*R[,,t-1]
    if(t>nn){
      lL_ht[m,(t-nn)]= mvnfast::dmvt(tdata[t,], rep(0,dm),R[,,t],nu[m],log=TRUE)+sum(INL.N[t,])-sum(INL.T[t,])
    }
  }
}

# should be the same
Rpred[[post.sample]][1:5]
R[,,nn+1][1:5]

sum(lL_static)
sum(lL_rmf)
sum(apply(lL_dcc,2,median))
sum(apply(lL_tdcc,2,median))
sum(apply(lL_ht,2,median))

# at the median of estimated parameters

a  = median(a)
b  = median(b)
nu = median(nu)
tdata  <- qt(udata,nu)
Rbar   <- cor(tdata[1:nn,])

R.ht = array(NA,c(dm, dm, nn+K))
R.ht[,,1] <- Rbar
tmp.ht = rep(NA,nn+K)
INL   <- dt(tdata,df=nu,log=TRUE)

for(t in 2:(nn+K)){
  R.ht[,,t] = (1-b)*Rbar-a*Pbar+a*Sigma[[t-1]]+b*R.ht[,,t-1]
  tmp.ht[t]= mvnfast::dmvt(tdata[t,], rep(0,dm),R.ht[,,t], nu, log=TRUE)+sum(INL.N[t,])-sum(INL[t,])
}

sum(tail(tmp.ht,K))
sum(apply(lL_ht,2,median))

roll_corr <- rollapply(data = cbind(data[,p1], data[,p2]), width = 100,
                       function(z) cor(z[,1], z[,2]), by.column = FALSE,
                       align = "right",fill=NA)
plot(roll_corr,type='l',ylim=c(-1,1))
lines(R.rmf[p1,p2,],type='l',ylim=c(0,1))
lines(R.dcc[p1,p2,],col=2)
lines(R.dcct[p1,p2,],col=4)
lines(R.ht[p1,p2,],col=3)

##------
## XM
##------

Sbar   = Reduce('+',Sig[1:nn])/nn
iota   = rep(1,dm)

load('empirical/temp/results_xm.Rdata')
M   = dim(res$r)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth
lL_xm = matrix(NA,ncol=K,nrow=post.sample)

lag = res$r[ind,1]
nu  = res$r[ind,2]
b1  = res$r[ind,3:(dm+2)]
b2  = res$r[ind,(dm+3):(2*dm+2)]

for(m in 1:post.sample){
  Vpred = vector(mode = "list", length = K)
  B1  = Outer(b1[m,],b1[m,])
  B2  = Outer(b2[m,],b2[m,])
  B0  = (Oiota-B1-B2)*Sbar
  
  for(t0 in 1:K){
    Vpred[[t0]] = B0+B1*Sig[[nn+t0-1]]+B2*(Reduce('+',Sig[(nn+t0-lag[m]):(nn+t0-1)])/lag[m])
    mtdf = nu[m]-dm
    lL_xm[m,t0] = mvnfast::dmvt(data[nn+t0,], rep(0,dm), (mtdf-1)/(mtdf+1)*Vpred[[t0]], df = mtdf+1, log=TRUE)
  }
}

# should be the same! ok
res$Vpred[ind][[post.sample]][1:5]
Vpred[[1]][1:5]

sum(lL_static)
sum(lL_rmf)
sum(apply(lL_dcc,2,median))
sum(apply(lL_tdcc,2,median))
sum(apply(lL_ht,2,median))
sum(apply(lL_xm,2,median))

# at the median of estimated parameters

B1  = Outer(apply(b1,2,median),apply(b1,2,median))
B2  = Outer(apply(b2,2,median),apply(b2,2,median))
B0  = (Oiota-B1-B2)*Sbar
R.xm = array(NA,c(dm, dm, nn+K))
tmp.xm  = rep(NA,nn+K)
tmp.true = rep(NA,nn+K)
mtdf = median(nu)-dm

for(t in 2:(nn+K)){
  R.xm[,,t]   <- B0+B1*Sig[[t-1]]+B2*(Reduce('+',Sig[max(c(t-median(lag)),1):(t-1)])/min(c(t-1,median(lag))))
  tmp.xm[t] = mvnfast::dmvt(data[t,], rep(0,dm), (mtdf-1)/(mtdf+1)*R.xm[,,t], df = mtdf+1, log=TRUE)
  tmp.true[t] = mvnfast::dmvn(data[t,], rep(0,dm),RCor[,,t], log=TRUE)
}

sum(tail(tmp.xm,K))
sum(apply(lL_xm,2,median))

roll_corr <- rollapply(data = cbind(data[,p1], data[,p2]), width = 100,
                       function(z) cor(z[,1], z[,2]), by.column = FALSE,
                       align = "right",fill=NA)
plot(tail(RCor[p1,p2,],K),col='gray80',type='l',ylim=c(-1,1))
lines(tail(roll_corr,K),lwd=3)
lines(tail(R.dcct[p1,p2,],K),col=2,lwd=2)
lines(tail(R.xm[p1,p2,],K),col=4,lwd=2)
lines(tail(R.ht[p1,p2,],K),col=3,lwd=2)

##------
## CAW-iw
##------

Sbar   = Reduce('+',Sig[1:nn])/nn
iota   = rep(1,dm)

load('empirical/temp/results_caw_iw.Rdata')
M   = dim(res$r)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth
lL_caw = matrix(NA,ncol=K,nrow=post.sample)

nu  = res$r[ind,1]
b1  = res$r[ind,2:(dm+1)]
b2  = res$r[ind,(dm+2):(2*dm+1)]
Vlast = res$Vpred[ind]
R = vector(mode = "list", length = K+nn)
R[[1]] = Sbar

for(m in 1:post.sample){
  B1  = Outer(b1[m,],b1[m,])
  B2  = Outer(b2[m,],b2[m,])
  B0  = (Oiota-B1-B2)*Sbar
  mtdf = nu[m]-dm
  
  for(t in 2:(K+nn)){
    R[[t]] = B0+B1*Sig[[t-1]]+B2*R[[t-1]]
    if(t>nn){
      lL_caw[m,t-nn] = mvnfast::dmvt(data[t,], rep(0,dm), (mtdf-1)/(mtdf+1)*R[[t]], df = mtdf+1, log=TRUE)
    }
  }
}

# should be the same!
R[[nn+1]][1:5]
res$Vpred[ind][[post.sample]][1:5]

sum(lL_static)
sum(lL_rmf)
sum(apply(lL_dcc,2,median))
sum(apply(lL_tdcc,2,median))
sum(apply(lL_ht,2,median))
sum(apply(lL_xm,2,median))
sum(apply(lL_caw,2,median))

# at the median of estimated parameters

B1  = Outer(apply(b1,2,median),apply(b1,2,median))
B2  = Outer(apply(b2,2,median),apply(b2,2,median))
B0  = (Oiota-B1-B2)*Sbar
R.caw = array(NA,c(dm, dm, nn+K))
tmp.caw  = rep(NA,nn+K)
mtdf = median(nu)-dm
R.caw[,,1] = Sbar

for(t in 2:(nn+K)){
  R.caw[,,t]   <- B0+B1*Sig[[t-1]]+B2*R.caw[,,t-1]
  tmp.caw[t] = mvnfast::dmvt(data[t,], rep(0,dm), (mtdf-1)/(mtdf+1)*R.caw[,,t], df = mtdf+1, log=TRUE)
}

sum(tail(tmp.caw,K))
sum(apply(lL_caw,2,median))

roll_corr <- rollapply(data = cbind(data[,p1], data[,p2]), width = 100,
                       function(z) cor(z[,1], z[,2]), by.column = FALSE,
                       align = "right",fill=NA)
plot(tail(RCor[p1,p2,],K),col='gray80',type='l',ylim=c(-1,1))
lines(tail(roll_corr,K),lwd=3)
lines(tail(R.dcct[p1,p2,],K),col=2,lwd=2)
lines(tail(R.xm[p1,p2,],K),col=4,lwd=2)
lines(tail(R.ht[p1,p2,],K),col=3,lwd=2)
lines(tail(R.caw[p1,p2,],K),col=6,lwd=2)

##----------------------
## save here
##----------------------
rm(res,Sig,R,Q)
save.image("empirical/temp/individual_llh.RData")
# load("empirical/temp/individual_llh.RData")


ALLRC = matrix(NA,nrow=50,ncol=nn+K)
for(j in 1:50){
  ALLRC[j,] = rollapply(data = cbind(data[,p1], data[,p2]), width = 25+j*2,
                       function(z) cor(z[,1], z[,2]), by.column = FALSE,
                       align = "right",fill=NA)
}

roll_corr <- rollapply(data = cbind(data[,p1], data[,p2]), width = 100,
                       function(z) cor(z[,1], z[,2]), by.column = FALSE,
                       align = "right",fill=NA)

plot(tail(date,K),tail(ALLRC[1,],K),ylim=c(-1,1),col='gray90',type='l')
for(j in 2:50) lines(tail(date,K),tail(ALLRC[j,],K),col='gray90')


lines(tail(date,K),tail(R.caw[p1,p2,],K),col='pink',lwd=2)
lines(tail(date,K),tail(R.xm[p1,p2,],K),col=4)
lines(tail(date,K),tail(R.dcct[p1,p2,],K),col=6,lwd=2)
lines(tail(date,K),tail(R.ht[p1,p2,],K),col=2,lwd=2)


####

mkvol = apply(tail(apply(RVs^2,2,scale),K),1,mean)

par(mfrow=c(1,1))
par(mfrow=c(1,1), mar=c(3, 3, 1, 1) + 0.1)
plot(tail(date,K),mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='', xaxt="n")

par(new = TRUE)
plot(tail(date,K),cumsum(apply(lL_tdcc,2,median))-
       cumsum(apply(lL_dcc,2,median)),type='l')


##------
## All models separately
##------

pdf('tables_and_figures/all_bfs_FX.pdf',height=4,width=10)
par(mfrow=c(1,1), mar=c(4, 3, 1, 1) + 0.1)
plot(tail(date,K),mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='', xaxt="n")

par(new = TRUE)
plot(tail(date,K),rep(0,K), xaxt="n",ylim=c(-5,2500),ylab='',xlab='',type='l')
lines(tail(date,K),cumsum(lL_rmf[,])-cumsum(lL_static[,]),lwd=2,lty=1,col=1)
lines(tail(date,K),cumsum(apply(lL_dcc[,],2,median))-cumsum(lL_static[,]),lwd=2,lty=2,col='coral')
lines(tail(date,K),cumsum(apply(lL_tdcc[,],2,median))-cumsum(lL_static[,]),lwd=2,col='coral3')
lines(tail(date,K),cumsum(apply(lL_xm,2,median))-cumsum(lL_static[,]),lwd=2,lty=2,col='royalblue')
lines(tail(date,K),cumsum(apply(lL_caw,2,median))-cumsum(lL_static[,]),lwd=2,col='royalblue4')
lines(tail(date,K),cumsum(apply(lL_ht,2,median))-cumsum(lL_static[,]),lwd=2,lty=1,col='violet')

legend(x=date[nn]-5,y=2500,col=c(1,'coral','coral3','royalblue','royalblue4','violet','gray90'),
       lty=c(1,2,1,2,1,1,1),lwd=c(2,2,2,2,2,2,2),
       legend=c('RMf','DCC','DCC-t','AIW','CAW','DCC-HEAVY-t','avrg.stand.RV'))
atx <- seq(date[(nn+1)], date[(nn+K)], by=25)
axis(1, at=atx, labels=format(atx, "%Y/%m"),las=2,cex.axis=0.75)
dev.off()

lpbf = c(sum(lL_static[,]),
         sum(lL_rmf[,]),
         sum(apply(lL_dcc[,],2,median)),
         sum(apply(lL_tdcc[,],2,median)),
         sum(apply(lL_xm[,],2,median)),
         sum(apply(lL_caw[,],2,median)),
         sum(apply(lL_ht,2,median)))
names(lpbf) = c('Static','RMf','DCC','DCC-t','AIW','CAW','DCC-HEAVY-t')
lpbfdf = as.data.frame(t(lpbf))


print(xtable(lpbfdf,align= 'cccccccc',
             caption = '1-step-ahead log predictive scores ($LPS$) 
             for all individual models: Static, RiskMetrics fixed (RMf),
rank-1 Dynamic conditional correlation with Gaussian and $t$ copulas (DCC
and DCC-t), Additive Inverse Wishart (AIW), Conditional Autoregressive Wishart (CAW)
and DCC-HEAVY model with $t$ copula for 
2020/01/02 - 2023/01/31 out-of-sample period
(K = 797 observations).',
             label = 'table:lps_FX', digits = 2),
      file='tables_and_figures/lps_FX.tex',
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

library(DMFP)
resDN= PMCMC_delNegro(lL_caw,lL_tdcc,1000,c(0,2),10000,propsd = 0.5)

ws_DN = t(pnorm(resDN$weights_xs))

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

if (dim(ws_DN)[1] > post.sample){
  dnind = round(seq(1,dim(ws_DN)[1],length=post.sample))
  ws_DN = ws_DN[dnind,]
}

resDNxs = list()

resDNxs$beta = resDN$beta[dnind]
resDNxs$weights_xs=resDN$weights_xs[,dnind]
resDNxs$likelihood=resDN$likelihood[,dnind]
resDNxs$acc=resDN$acc[dnind]


##----------------------
## save here
##----------------------
rm(resDN,R.caw,R.dcc,R.dcct,R.ht,R.rmf,R.xm,ws_jore5,ws_jore25)
save.image("empirical/temp/DNweights.RData")
# load("empirical/temp/DNweights.RData")

lhf = exp(lL_caw)
llf = exp(lL_tdcc)

for(m in 1:post.sample){
  for (t in 1:K){
    # weights geweke
    a1p = lhf[m,1:t]
    a2p = llf[m,1:t]
    opw = function(x){
      -sum(log(x*a1p+(1-x)*a2p))
    }
    ws_gew[m,t] = optim(0.5, opw, gr = NULL,
                      method = c("L-BFGS-B"),
                      lower = 0, upper = 1, hessian = FALSE)$par
    lL_gew[m,t] = log(ws_gew[m,t]*lhf[m,t]+(1-ws_gew[m,t])*llf[m,t])
    
    # DN
    
    lL_DN[m,t] = log(ws_DN[m,t]*lhf[m,t]+(1-ws_DN[m,t])*llf[m,t])
    
    # jore 1 past
    ws_jore1[m,t] = a1p[t]/(a1p[t]+a2p[t])
    lL_jore1[m,t] = log(ws_jore1[m,t]*lhf[m,t]+(1-ws_jore1[m,t])*llf[m,t])
    
    # # jore 5 past
    # pst = 5
    # ws_jore5[m,t] = exp(sum(log(a1p[max(t-pst+1,1):t])))/
    #   (exp(sum(log(a1p[max(t-pst+1,1):t])))+exp(sum(log(a2p[max(t-pst+1,1):t]))))
    # lL_jore5[m,t] = log(ws_jore5[m,t]*lhf[m,t]+(1-ws_jore5[m,t])*llf[m,t])
    
    # jore 10 past
    pst = 10
    ws_jore10[m,t] = exp(sum(log(a1p[max(t-pst+1,1):t])))/
      (exp(sum(log(a1p[max(t-pst+1,1):t])))+exp(sum(log(a2p[max(t-pst+1,1):t]))))
    lL_jore10[m,t] = log(ws_jore10[m,t]*lhf[m,t]+(1-ws_jore10[m,t])*llf[m,t])
    
    # # jore 25 past
    # pst = 25
    # ws_jore25[m,t] = exp(sum(log(a1p[max(t-pst+1,1):t])))/
    #   (exp(sum(log(a1p[max(t-pst+1,1):t])))+exp(sum(log(a2p[max(t-pst+1,1):t]))))
    # lL_jore25[m,t] = log(ws_jore25[m,t]*lhf[m,t]+(1-ws_jore25[m,t])*llf[m,t])
    
    # equally-weighted
    lL_equal[m,t] = log(0.5*lhf[m,t]+0.5*llf[m,t])
  }
}

move.axis = 100

pdf('tables_and_figures/weights_FX.pdf',height=8,width=14)
par(mfrow=c(2,1), mar=c(4, 3, 1, 1) + 0.1)
plot(date[(nn+1):(nn+K)],mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]))
par(new = TRUE)
plot(date[(nn+1):(nn+K)],apply(ws_jore1,2,median),ylim=c(0,1),
     type='l',ylab='',xlab='',xaxt="n",lwd=3,xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]),col='pink')
lines(date[(nn+1):(nn+K)],apply(ws_jore10,2,median),col='pink4',lwd=2,lty=2)
lines(date[(nn+1):(nn+K)],apply(ws_DN,2,median),lty=1,lwd=2,col='blue')
lines(date[(nn+1):(nn+K)],apply(ws_gew,2,median),col=1,lwd=2,lty=1)
abline(h=0.5)
axis(1, at=atx, labels=format(atx, "%Y/%m"),las=2,cex.axis=0.75)
legend(x=date[(nn+1)]-move.axis-35,y=1,col=c('pink','pink4','blue',1),
       lty=c(1,2,1,1),lwd=c(2,2,2,2),
       legend=c('Jore1','Jore10','DelNegro','Geweke'))


plot(tail(date,K),mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]))
par(new = TRUE)

best.model = lL_caw

plot(tail(date,K),cumsum(apply(best.model,2,median))-cumsum(apply(best.model,2,median)),
     type='l',ylim=c(-500,500),ylab='',xlab='',xaxt="n",
     xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]))
lines(tail(date,K),cumsum(apply(lL_tdcc,2,median))-cumsum(apply(best.model,2,median)),lwd=2,
      col='coral3')
lines(tail(date,K),cumsum(apply(lL_gew,2,median))-cumsum(apply(best.model,2,median)),
      col=1,lwd=2)
lines(tail(date,K),cumsum(apply(lL_jore1,2,median))-cumsum(apply(best.model,2,median)),
      col='pink',lwd=2)
lines(tail(date,K),cumsum(apply(lL_jore10,2,median))-cumsum(apply(best.model,2,median)),
      col='pink4',lwd=2,lty=2)
lines(tail(date,K),cumsum(apply(lL_DN,2,median))-cumsum(apply(best.model,2,median)),lwd=2,
      col='blue')
lines(tail(date,K),cumsum(apply((lL_equal),2,median))-cumsum(apply(best.model,2,median)),
      col='coral',lwd=2,lty=2)
lines(tail(date,K),cumsum(apply(lL_ht,2,median))-cumsum(apply(best.model,2,median)),
      col='violet',lwd=2)
axis(1, at=atx, labels=format(atx, "%Y/%m"),las=2,cex.axis=0.75)
legend(x=date[(nn+1)]-move.axis-35,y=500,col=c('coral3',1,'pink','pink4','blue','coral','violet'),
       lty=c(1,1,1,2,1,2,1),lwd=c(2,2,2,2,2,2,2),
       legend=c('DCC-t','Geweke','Jore1','Jore10','DelNegro','Equal','DCC-HEAVY-t'))
dev.off()


# -----------------------------
# LPTS
# -----------------------------

# first, find the quantiles Q=50,25,10,5

cuts=seq(0.6,1,length=5000)
perc = rep(NA,length(cuts))
for(i in 1:length(cuts)){
  perc[i] = mean(apply(tail(udata,K)<cuts[i],1,all))
}

ind50 = which(apply(tail(udata,K)<cuts[which.min(abs(perc - 0.5))],1,all))
ind25 = which(apply(tail(udata,K)<cuts[which.min(abs(perc - 0.25))],1,all))
ind10 = which(apply(tail(udata,K)<cuts[which.min(abs(perc - 0.10))],1,all))
ind05 = which(apply(tail(udata,K)<cuts[which.min(abs(perc - 0.05))],1,all))

pdf('tables_and_figures/lpts_FX.pdf',height=8,width=12)
par(mfrow=c(2,2))
plot(density(apply(lL_caw, 1,mean)),xlim=c(-4,-2.9),ylab='',xlab='',main='LPS',
     lwd=2,lty=1,col='royalblue')
lines(density(apply(lL_tdcc, 1,mean)),lwd=2,col='coral3')
lines(density(apply(lL_ht, 1,mean)),lwd=2,col='violet')
lines(density(apply(lL_gew, 1,mean)),lwd=2,col='black')
lines(density(apply(lL_jore1, 1,mean)),lwd=2,col='pink')
lines(density(apply(lL_equal, 1,mean)),lwd=2,col='coral',lty=2)
legend(x=-4,y=110,lwd=c(2,2,2,2,2,2),lty=c(1,1,1,1,1,2),
       col=c('coral3','royalblue','violet','black','pink','coral'),
       legend=c('DCC-t','CAW','DCC-HEAVY-t','Geweke','Jore1','Equal'))

plot(density(apply(lL_caw[,ind25], 1,mean)),xlim=c(-0.65,0.25),ylab='',xlab='',main='LPTS(25%)',
     lwd=2,lty=1,col='royalblue')
lines(density(apply(lL_tdcc[,ind25], 1,mean)),lwd=2,col='coral3')
lines(density(apply(lL_ht[,ind25], 1,mean)),lwd=2,col='violet')
lines(density(apply(lL_gew[,ind25], 1,mean)),lwd=2,col='black')
lines(density(apply(lL_jore1[,ind25], 1,mean)),lwd=2,col='pink')
lines(density(apply(lL_equal[,ind25], 1,mean)),lwd=2,col='coral',lty=2)
legend(x=-0.65,y=90,lwd=c(2,2,2,2,2,2),lty=c(1,1,1,1,1,2),
       col=c('coral3','royalblue','violet','black','pink','coral'),
       legend=c('DCC-t','CAW','DCC-HEAVY-t','Geweke','Jore1','Equal'))

plot(density(apply(lL_caw[,ind10], 1,mean)),xlim=c(0.2,1),ylab='',xlab='',main='LPTS(10%)',
     lwd=2,lty=1,col='royalblue')
lines(density(apply(lL_tdcc[,ind10], 1,mean)),lwd=2,col='coral3')
lines(density(apply(lL_ht[,ind10], 1,mean)),lwd=2,col='violet')
lines(density(apply(lL_gew[,ind10], 1,mean)),lwd=2,col='black')
lines(density(apply(lL_jore1[,ind10], 1,mean)),lwd=2,col='pink')
lines(density(apply(lL_equal[,ind10], 1,mean)),lwd=2,col='coral',lty=2)
legend(x=0.2,y=90,lwd=c(2,2,2,2,2,2),lty=c(1,1,1,1,1,2),
       col=c('coral3','royalblue','violet','black','pink','coral'),
       legend=c('DCC-t','CAW','DCC-HEAVY-t','Geweke','Jore1','Equal'))

plot(density(apply(lL_caw[,ind05], 1,mean)),xlim=c(0.4,1.2),ylab='',xlab='',main='LPTS(5%)',
     lwd=2,lty=1,col='royalblue')
lines(density(apply(lL_tdcc[,ind05], 1,mean)),lwd=2,col='coral3')
lines(density(apply(lL_ht[,ind05], 1,mean)),lwd=2,col='violet')
lines(density(apply(lL_gew[,ind05], 1,mean)),lwd=2,col='black')
lines(density(apply(lL_jore1[,ind05], 1,mean)),lwd=2,col='pink')
lines(density(apply(lL_equal[,ind05], 1,mean)),lwd=2,col='coral',lty=2)
legend(x=0.4,y=90,lwd=c(2,2,2,2,2,2),lty=c(1,1,1,1,1,2),
       col=c('coral3','royalblue','violet','black','pink','coral'),
       legend=c('DCC-t','CAW','DCC-HEAVY-t','Geweke','Jore1','Equal'))
dev.off()

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

marketvols = cbind(mkvol,avrg2,coredata(vix))
cor(marketvols,use = "pairwise.complete.obs")
par(mfrow=c(1,1))
plot(tail(date,K),marketvols[,1],type='l')
lines(tail(date,K),marketvols[,2],col=2)
lines(tail(date,K),sqrt(marketvols[,3]),col=3)



vols.ws = cbind(marketvols,apply(ws_gew,2,median),
                apply(ws_jore1,2,median),apply(ws_DN,2,median),apply(lL_caw-lL_tdcc,2,median))
dim(vols.ws)
colnames(vols.ws) = c('avrg RV','Mkt:eql','VIX','Geweke','Jore1','DelNegro','logLik diff')
cor(vols.ws,use = "pairwise.complete.obs")


corrs_res = cor(vols.ws,use = "pairwise.complete.obs")
corrs_res[upper.tri(corrs_res)] = NA

print(xtable::xtable(corrs_res,
                     caption = "Sample correlations 
             between three proxies for the   market volatility and
             the preference for high-frequency model for
             2020/01/02 to 2023/01/31 out-of-sample period ($K=797$ observations).
             The proxies for the market volatility are: average standardized 
             realized volatility (avrg RV), 
             equally weighted market portfolio realized volatility (Mkt:eql) 
             and VIX index. The preference for the high-frequency model is measured as a 
             high-frequency 
             component weight in various  pooling schemes 
             as well as the difference between the daily log likelihood (logLik diff) 
             between the CAW and DCC-t models.",
                     label = 'table:corrs_FX', digits = 3),
      file='tables_and_figures/corrs_FX.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,scalebox = 0.8)

##
save(ws_DN,ws_gew,ws_jore1,K,dm,post.sample,ind,nn,file="empirical/temp/res_FX.Rdata")


