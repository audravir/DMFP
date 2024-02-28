rm(list=ls(all=TRUE))
load('data/FXdata.Rdata')

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
post.sample = 100
p1 = 1
p2 = 2


#-----------------------
# dcc-HEAVY-scalar Normal Copula
#-----------------------

load('empirical/temp/results_heavy_normal.Rdata')
M     = dim(res$r)[1]
ind   = round(seq(1,M,length=post.sample)) #thin every xth
lL_hn = matrix(NA,ncol=K,nrow=post.sample)

a = res$r[ind,1]
b = res$r[ind,2]

Rpred = res$Rpred[ind]
Pbar  = Reduce('+',Sig[1:nn])/nn
Rbar  = cor(data[1:nn,])
R     = array(NA,c(dm, dm, nn+K))
R[,,1] = Rbar

for(m in 1:post.sample){
  for(t in 2:(nn+K)){
    R[,,t] = (1-b[m])*Rbar-a[m]*Pbar+a[m]*Sig[[t-1]]+b[m]*R[,,t-1]
      if(t>nn){
        lL_hn[m,(t-nn)]= mvnfast::dmvn(data[t,], rep(0,dm),R[,,t],log=TRUE)
    }
  }
}

# should be the same
Rpred[[post.sample]][1:5]
R[,,nn+1][1:5]

sum(apply(lL_hn,2,median))


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
Pbar  = Reduce('+',Sig[1:nn])/nn
R     = array(NA,c(dm, dm, nn+K))
INL.N = dnorm(data,log=TRUE)

for(m in 1:post.sample){
  
  tdata  = qt(udata,nu[m])
  Rbar   = cor(tdata[1:nn,])
  R[,,1] = Rbar
  INL.T  = dt(tdata,df=nu[m],log=TRUE)
  
  for(t in 2:(nn+K)){
    R[,,t] = (1-b[m])*Rbar-a[m]*Pbar+a[m]*Sig[[t-1]]+b[m]*R[,,t-1]
    if(t>nn){
      lL_ht[m,(t-nn)]= mvnfast::dmvt(tdata[t,], rep(0,dm),R[,,t],nu[m],log=TRUE)+sum(INL.N[t,])-sum(INL.T[t,])
    }
  }
}

# should be the same
Rpred[[post.sample]][1:5]
R[,,nn+1][1:5]

sum(apply(lL_hn,2,median))
sum(apply(lL_ht,2,median))




#-----------------------
# dcc-HEAVY-scalar (joint) t Copula
#-----------------------

load('empirical/temp/results_heavy_scalar_joint.Rdata')
M   = dim(res$r)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth
lL_hsj = matrix(NA,ncol=K,nrow=post.sample)

nu = res$r[ind,1]
a  = res$r[ind,2]
b  = res$r[ind,3]
Rpred = res$Rpred[ind]
Pbar = Reduce('+',Sig[1:nn])/nn

# one way
for(m in 1:post.sample){
  Vpred = vector(mode = "list", length = K)
  tdata  <- qt(udata,nu[m])
  Rbar   <- cor(tdata[1:nn,])

  Vpred[[1]] = Rpred[[m]]
  lL_hsj[m,1] = mvnfast::dmvt(tdata[nn+1,], rep(0,dm), Vpred[[1]], df = nu[m], log=TRUE)-
    sum(dt(tdata[nn+1,],df=nu[m],log=TRUE))+sum(dnorm(data[nn+1,],log=TRUE))

  for(t0 in 2:K){
    Vpred[[t0]] = Rbar+a[m]*(Sig[[nn+t0-1]]-Pbar)+b[m]*(Vpred[[t0-1]]-Rbar)
    lL_hsj[m,t0] = mvnfast::dmvt(tdata[nn+t0,], rep(0,dm), Vpred[[t0]], df = nu[m], log=TRUE)-
      sum(dt(tdata[nn+t0,],df=nu[m],log=TRUE))+sum(dnorm(data[nn+t0,],log=TRUE))
  }
}

# second way
lL_hsj2 = matrix(NA,ncol=K,nrow=post.sample)
Pbar = Reduce('+',Sig[1:nn])/nn
R    = array(NA,c(dm, dm, nn+K))
INL.N    = dnorm(data,log=TRUE)

for(m in 1:post.sample){
  tdata  <- qt(udata,nu[m])
  Rbar   <- cor(tdata[1:nn,])
  R[,,1] <- Rbar
  INL   <- dt(tdata,df=nu[m],log=TRUE)

  for(t in 2:(nn+K)){
    R[,,t] = Rbar+a[m]*(Sig[[t-1]]-Pbar)+b[m]*(R[,,t-1]-Rbar)
    if(t>nn){
      lL_hsj2[m,(t-nn)]= mvnfast::dmvt(tdata[t,], rep(0,dm),R[,,t], nu[m], log=TRUE)+sum(INL.N[t,])-sum(INL[t,])
    }
  }
}

# compare
R[,,nn+1][1:3,1:3]
Vpred[[1]][1:3,1:3]

sum(apply(lL_hsj,2,median))
sum(apply(lL_hsj2,2,median))

par(mfrow=c(2,2))
plot(nu,type='l')
plot(a,type='l')
plot(b,type='l')

# at the median of estimated parameters

nu = median(res$r[ind,1])
a  = median(res$r[ind,2])
b  = median(res$r[ind,3])
R_hsj = array(NA,c(dm, dm, nn+K))
R_hsj[,,1] <- Rbar
tmp.hsj <- rep(NA,nn+K)
tdata   <- qt(udata,nu)
Rbar    <- cor(tdata[1:nn,])
INL     <- dt(tdata,df=nu,log=TRUE)

for(t in 2:(nn+K)){
  R_hsj[,,t] = Rbar+a*(Sig[[t-1]]-Pbar)+b*(R_hsj[,,t-1]-Rbar)
  tmp.hsj[t] = mvnfast::dmvt(tdata[t,], rep(0,dm),R_hsj[,,t], nu, log=TRUE)+sum(INL.N[t,])-sum(INL[t,])
}

sum(tail(tmp.hsj,K))
sum(apply(lL_hsj,2,median))

avrgmatrix = Reduce('+',res$Rpred[ind])/post.sample
avrgmatrix[1:3,1:3]
R_hsj[1:3,1:3,nn+1]

##-----------------------
## dcc-HEAVY-scalar (separate) t Copula
##-----------------------

load('empirical/temp/results_heavy_scalar_separate.Rdata')
M   = dim(res$r)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth
lL_hss = matrix(NA,ncol=K,nrow=post.sample)

nu = res$r[ind,1]
a  = res$r[ind,2]
b  = res$r[ind,3]
Rpred = res$Rpred[ind]
Pbar = Reduce('+',Sig[1:nn])/nn

# one way
for(m in 1:post.sample){
  Vpred = vector(mode = "list", length = K)
  tdata  <- qt(udata,nu[m])
  Rbar   <- cor(tdata[1:nn,])
  
  Vpred[[1]] = Rpred[[m]]
  lL_hss[m,1] = mvnfast::dmvt(tdata[nn+1,], rep(0,dm), Vpred[[1]], df = nu[m], log=TRUE)-
    sum(dt(tdata[nn+1,],df=nu[m],log=TRUE))+sum(dnorm(data[nn+1,],log=TRUE))
  
  for(t0 in 2:K){
    Vpred[[t0]] = Rbar+a[m]*(Sig[[nn+t0-1]]-Pbar)+b[m]*(Vpred[[t0-1]]-Rbar)
    lL_hss[m,t0] = mvnfast::dmvt(tdata[nn+t0,], rep(0,dm), Vpred[[t0]], df = nu[m], log=TRUE)-
      sum(dt(tdata[nn+t0,],df=nu[m],log=TRUE))+sum(dnorm(data[nn+t0,],log=TRUE))
  }
}

# second way
lL_hss2 = matrix(NA,ncol=K,nrow=post.sample)
Pbar = Reduce('+',Sig[1:nn])/nn
R    = array(NA,c(dm, dm, nn+K))
INL.N = dnorm(data,log=TRUE)

for(m in 1:post.sample){
  tdata  <- qt(udata,nu[m])
  Rbar   <- cor(tdata[1:nn,])
  R[,,1] <- Rbar
  INL   <- dt(tdata,df=nu[m],log=TRUE)
  
  for(t in 2:(nn+K)){
    R[,,t] = Rbar+a[m]*(Sig[[t-1]]-Pbar)+b[m]*(R[,,t-1]-Rbar)
    if(t>nn){
      lL_hss2[m,(t-nn)]= mvnfast::dmvt(tdata[t,], rep(0,dm),R[,,t], nu[m], log=TRUE)+sum(INL.N[t,])-sum(INL[t,])
    }
  }
}

# compare
R[,,nn+1][1:3,1:3]
Vpred[[1]][1:3,1:3]

sum(apply(lL_hss,2,median))
sum(apply(lL_hss2,2,median))

par(mfrow=c(2,2))
plot(nu,type='l')
plot(a,type='l')
plot(b,type='l')

# at the median of estimated parameters

nu = median(res$r[ind,1])
a  = median(res$r[ind,2])
b  = median(res$r[ind,3])
R_hss = array(NA,c(dm, dm, nn+K))
R_hss[,,1] <- Rbar
tmp.hss <- rep(NA,nn+K)
tdata   <- qt(udata,nu)
Rbar    <- cor(tdata[1:nn,])
INL     <- dt(tdata,df=nu,log=TRUE)

for(t in 2:(nn+K)){
  R_hss[,,t] = Rbar+a*(Sig[[t-1]]-Pbar)+b*(R_hss[,,t-1]-Rbar)
  tmp.hss[t] = mvnfast::dmvt(tdata[t,], rep(0,dm),R_hss[,,t], nu, log=TRUE)+sum(INL.N[t,])-sum(INL[t,])
}

sum(tail(tmp.hss,K))
sum(apply(lL_hss,2,median))

avrgmatrix = Reduce('+',res$Rpred[ind])/post.sample
avrgmatrix[1:3,1:3]
R_hss[1:3,1:3,nn+1]

##-----------------------
## dcc-HEAVY-matrix (true) t Copula
##-----------------------

load("empirical/temp/results_heavy_true.Rdata")

M   = dim(res$r)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth
lL_hm = matrix(NA,ncol=K,nrow=post.sample)

nu = res$r[ind,1]
a  = res$r[ind,2:(dm+1)]
b  = res$r[ind,(dm+2):(dm*2+1)]
Rpred = res$Rpred[ind]
Pbar  = Reduce('+',Sig[1:nn])/nn
iota  = rep(1,dm)
Oiota = Outer(iota,iota)

# one way
for(m in 1:post.sample){
  Vpred = vector(mode = "list", length = K)
  tdata  <- qt(udata,nu[m])
  Rbar   <- cor(tdata[1:nn,])
  A      = Outer(a[m,],a[m,])
  B      = Outer(b[m,],b[m,])
  Rtilde = (Oiota-B)*Rbar-A*Pbar
  
  Vpred[[1]] = Rpred[[m]]
  lL_hm[m,1] = mvnfast::dmvt(tdata[nn+1,], rep(0,dm), Vpred[[1]], df = nu[m], log=TRUE)-
    sum(dt(tdata[nn+1,],df=nu[m],log=TRUE))+sum(dnorm(data[nn+1,],log=TRUE))
  
  for(t0 in 2:K){
    Vpred[[t0]] = Rtilde+A*Sig[[nn+t0-1]]+B*Vpred[[t0-1]]
    lL_hm[m,t0] = mvnfast::dmvt(tdata[nn+t0,], rep(0,dm), Vpred[[t0]], df = nu[m], log=TRUE)-
      sum(dt(tdata[nn+t0,],df=nu[m],log=TRUE))+sum(dnorm(data[nn+t0,],log=TRUE))
  }
}

# second way
lL_hm2 = matrix(NA,ncol=K,nrow=post.sample)
Pbar  = Reduce('+',Sig[1:nn])/nn
R     = array(NA,c(dm, dm, nn+K))
INL.N = dnorm(data,log=TRUE)

for(m in 1:post.sample){
  tdata  <- qt(udata,nu[m])
  Rbar   <- cor(tdata[1:nn,])
  R[,,1] <- Rbar
  INL   <- dt(tdata,df=nu[m],log=TRUE)
  A      = Outer(a[m,],a[m,])
  B      = Outer(b[m,],b[m,])
  Rtilde = (Oiota-B)*Rbar-A*Pbar
  
  for(t in 2:(nn+K)){
    R[,,t] = Rtilde+A*Sig[[t-1]]+B*R[,,t-1]
    if(t>nn){
      lL_hm2[m,(t-nn)]= mvnfast::dmvt(tdata[t,], rep(0,dm),R[,,t], nu[m], log=TRUE)+sum(INL.N[t,])-sum(INL[t,])
    }
  }
}

# compare
R[,,nn+1][1:3,1:3]
Vpred[[1]][1:3,1:3]

sum(apply(lL_hm,2,median))
sum(apply(lL_hm2,2,median))

# at the median of estimated parameters

A  = Outer(apply(a,2,median),apply(a,2,median))
B  = Outer(apply(b,2,median),apply(b,2,median))
nu = median(nu)
tdata  <- qt(udata,nu)
Rbar   <- cor(tdata[1:nn,])
Rtilde = (Oiota-B)*Rbar-A*Pbar
R_hm = array(NA,c(dm, dm, nn+K))
R_hm[,,1] <- Rbar
tmp.hm  = rep(NA,nn+K)
INL   <- dt(tdata,df=nu,log=TRUE)

for(t in 2:(nn+K)){
  R_hm[,,t] = Rtilde+A*Sig[[t-1]]+B*R_hm[,,t-1]
  tmp.hm[t]= mvnfast::dmvt(tdata[t,], rep(0,dm),R_hm[,,t], nu, log=TRUE)+sum(INL.N[t,])-sum(INL[t,])
}

sum(tail(tmp.hm,K))
sum(apply(lL_hm,2,median))

avrgmatrix = Reduce('+',res$Rpred[ind])/post.sample
avrgmatrix[1:3,1:3]
R_hm[1:3,1:3,nn+1]

##-----------------------
## Volatilities
##-----------------------

par(mfrow=c(1,1))
plot(tail(R_hsj[p1,p2,],K),type='l',ylim=c(-1,1))
lines(tail(R_hss[p1,p2,],K),col=2)
lines(tail(R_hm[p1,p2,],K),col=4,lwd=2)

# ##-----------------------
# ## dcc-HEAVY t Copula
# ##-----------------------
# 
# load('empirical/temp/results_heavy.Rdata')
# M   = dim(res$restdcch)[1]
# ind = round(seq(1,M,length=post.sample)) #thin every xth
# 
# R       = array(NA,c(dm, dm, nn+K))
# a      <- res$restdcch[ind,2:(dm+1)]
# b      <- res$restdcch[ind,(dm+2):(dm*2+1)]
# nu     <- res$restdcch[ind,1]
# 
# Pbar = Reduce('+',Sig[1:nn])/nn
# Rbar = cor(data[start:nn,])
# 
# R[,,1] <- Rbar
# 
# lL_tdcch = matrix(NA,ncol=K,nrow=length(ind))
# 
# for(m in 1:length(ind)){
#   tdata <- qt(udata,nu[m])
#   Rbar  <- cor(tdata[start:nn,])
#   A     <- diag (a[m,])
#   B     <- diag (b[m,])
#   INL   <- dt(tdata,nu[m],log=TRUE)
#   
#   for(t in 2:(nn+K)){
#     R[,,t]   <- Rbar+A*(Sig[[t-1]]-Pbar)+B*(R[,,t-1]-Rbar)
#     
#     if(t>nn){
#       lL_tdcch[m,(t-nn)] <- mvnfast::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nu[m], log=TRUE)+
#         sum(INL.N[t,])-sum(INL[t,])
#     }
#   }
# }
# 
# dim(lL_tdcch)
# 
# zoomin = c(1:100)
# 
# par(mfrow=c(1,1))
# plot(apply(lL_dcc,2,median)[zoomin],type='l')
# lines(apply(lL_tdcch,2,median)[zoomin],col=2)
# 
# sum(lL_static)
# sum(lL_rmf)
# sum(apply(lL_dcc,2,median))
# sum(apply(lL_tdcc,2,median))
# sum(apply(lL_tdcch,2,median))
# 
# ###############################
# 
# # at the median of estimated parameters
# 
# A       = diag(apply(a,2,median))
# B       = diag(apply(b,2,median))
# tdata   <- qt(udata,median(nu))
# Rbar    = cor(tdata[start:nn,])
# R.heavy = array(NA,c(dm, dm, nn+K))
# tmp.ll  = rep(NA,nn+K)
# R.heavy[,,1] <- Rbar
# INL   <- dt(tdata,median(nu),log=TRUE)
# 
# microbenchmark(A%*%(Sig[[t-1]]-Pbar),mat.mult(A,(Sig[[t-1]]-Pbar)),unit = 'relative',times = 1000)
# 
# mm1=A%*%(Sig[[t-1]]-Pbar)
# mm2=diag(A)*(Sig[[t-1]]-Pbar)
# 
# mm1-mm2
# 
# for(t in 2:(nn+K)){
#   R.heavy[,,t]   <- Rbar+A*(Sig[[t-1]]-Pbar)+B*(R.heavy[,,t-1]-Rbar)
#   tmp.ll[t] = mvnfast::dmvt(tdata[t,], rep(0,dm), R.heavy[,,t], df = median(nu), log=TRUE)+
#     sum(INL.N[t,])-sum(INL[t,])
# }
# 
# sum(tail(tmp.ll,K))
# 
# roll_corr <- rollapply(data = cbind(data[,p1], data[,p2]), width = 100,
#                        function(z) cor(z[,1], z[,2]), by.column = FALSE,
#                        align = "center",fill=NA)
# plot(roll_corr,type='l',ylim=c(-1,1))
# lines(R.rmf[p1,p2,],type='l',ylim=c(0,1))
# lines(R.dcc[p1,p2,],col=2)
# lines(R.heavy[p1,p2,],col=4)

##------
## XM
##------

Sbar   = Reduce('+',Sig[1:nn])/nn
iota   = rep(1,dm)

load('empirical/temp/results_xm.Rdata')
M   = dim(res$resc)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth
lL_xm1 = matrix(NA,ncol=K,nrow=post.sample)

lag = res$resc[ind,1]
nu  = res$resc[ind,2]
b1  = res$resc[ind,3:(dm+2)]
b2  = res$resc[ind,(dm+3):(2*dm+2)]

for(m in 1:post.sample){
  Vpred = vector(mode = "list", length = K)
  B1  = Outer(b1[m,],b1[m,])
  B2  = Outer(b2[m,],b2[m,])
  B0  = (Oiota-B1-B2)*Sbar
  
  for(t0 in 1:K){
    Vpred[[t0]] = B0+B1*Sig[[nn+t0-1]]+B2*(Reduce('+',Sig[(nn+t0-lag[m]):(nn+t0-1)])/lag[m])
    mtdf = nu[m]-dm
    lL_xm1[m,t0] = mvnfast::dmvt(data[nn+t0,], rep(0,dm), (mtdf-1)/(mtdf+1)*Vpred[[t0]], df = mtdf+1, log=TRUE)
  }
}

sum(lL_static)
sum(lL_rmf)
sum(apply(lL_dcc,2,median))
sum(apply(lL_tdcc,2,median))
sum(apply(lL_xm1,2,median))

# at the median of estimated parameters

B1  = Outer(apply(b1,2,median),apply(b1,2,median))
B2  = Outer(apply(b2,2,median),apply(b2,2,median))
B0  = (Oiota-B1-B2)*Sbar
R.xm = array(NA,c(dm, dm, nn+K))
tmp.xm  = rep(NA,nn+K)
tmp.true = rep(NA,nn+K)
mtdf = median(nu)-dm
tmp.dcc = rep(NA,nn+K)

for(t in 2:(nn+K)){
  R.xm[,,t]   <- B0+B1*Sig[[t-1]]+B2*(Reduce('+',Sig[max(c(t-median(lag)),1):(t-1)])/min(c(t-1,median(lag))))
  tmp.xm[t] = mvnfast::dmvt(data[t,], rep(0,dm), (mtdf-1)/(mtdf+1)*R.xm[,,t], df = mtdf+1, log=TRUE)
  tmp.true[t] = mvnfast::dmvn(data[t,], rep(0,dm),RCor[,,t], log=TRUE)
  tmp.dcc[t] = mvnfast::dmvn(data[t,], rep(0,dm), R.dcc[,,t], log=TRUE)
}

sum(tail(tmp.xm,K))
sum(tail(tmp.dcct,K))
sum(tail(tmp.dcc,K))
sum(tail(tmp.true,K))


roll_corr <- rollapply(data = cbind(data[,p1], data[,p2]), width = 100,
                       function(z) cor(z[,1], z[,2]), by.column = FALSE,
                       align = "center",fill=NA)
plot(RCor[p1,p2,],col=3,type='l',ylim=c(-1,1))
lines(roll_corr,lwd=3)
lines(R.dcct[p1,p2,],col=2)
lines(R.xm[p1,p2,],col=4)

##------
## CAW
##------

llNW = function(y,S,nu,iter=1000){
  fun=function(x) mvnfast::dmvn(y, rep(0,dm), x, log=TRUE)
  V = rWishart(iter,nu,S/nu)
  res = mean(apply(V,3,fun))
  return(res)
}

load('empirical/temp/results_caw.Rdata')
M   = dim(res$resc)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth
lL_caw = lL_cawN =matrix(NA,ncol=K,nrow=post.sample)

nu=res$resc[ind,1]
b1=res$resc[ind,2:(dm+1)]
b2=res$resc[ind,(dm+2):(dm*2+1)]
Vlast = res$Vpred[ind]

for(m in 1:post.sample){
  Vpred = vector(mode = "list", length = K)
  B1  = Outer(b1[m,],b1[m,])
  B2  = Outer(b2[m,],b2[m,])
  B0  = (Oiota-B1-B2)*Sbar
  
  Vpred[[1]] = B0+B1*Vlast[[m]]+B2*Sig[[nn]]
  # lL_caw[m,1] = llNW(data[nn+1,],Vpred[[1]],nu[m])
  lL_cawN[m,1] = mvnfast::dmvn(data[nn+1,], rep(0,dm), Vpred[[1]], log=TRUE)
  
  
  for(t0 in 2:K){
    Vpred[[t0]] = B0+B1*Vpred[[t0-1]]+B2*Sig[[nn+t0-1]]
    #    lL_caw[m,t0] = llNW(data[nn+t0,],Vpred[[t0]],nu[m])
    lL_cawN[m,t0] = mvnfast::dmvn(data[nn+t0,], rep(0,dm), Vpred[[t0]], log=TRUE)
  }
  
}

sum(lL_static)
sum(lL_rmf)
sum(apply(lL_dcc,2,median))
sum(apply(lL_tdcc,2,median))
sum(apply(lL_xm1,2,median))
sum(apply(lL_cawN,2,median))

# at the median of estimated parameters

B1  = Outer(apply(b1,2,median),apply(b1,2,median))
B2  = Outer(apply(b2,2,median),apply(b2,2,median))
B0  = (Oiota-B1-B2)*Sbar
R.caw = array(NA,c(dm, dm, nn+K))
tmp.caw  = rep(NA,nn+K)
R.caw[,,1] = Sbar

for(t in 2:(nn+K)){
  R.caw[,,t]   <- B0+B1*R.caw[,,t-1]+B2*Sig[[t-1]]
  tmp.caw[t] = mvnfast::dmvn(data[t,], rep(0,dm), R.caw[,,t], log=TRUE)
}

sum(tail(tmp.caw,K))

roll_corr <- rollapply(data = cbind(data[,p1], data[,p2]), width = 100,
                       function(z) cor(z[,1], z[,2]), by.column = FALSE,
                       align = "right",fill=NA)
plot(RCor[p1,p2,],col=3,type='l',ylim=c(-1,1))
lines(roll_corr,lwd=3)
lines(R.caw[p1,p2,],col=2)
lines(R.xm[p1,p2,],col=4)

p1=4
p2=6


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

#### COMPARE

dcct = cumsum(apply(lL_tdcc,2,mean))
xm   = cumsum(apply(lL_xm1,2,mean))

plot(tail(date,K),xm-dcct,type='l')

mean((lL_xm1-lL_tdcc)<0)

####


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

pdf('tables_and_figures/all_bfs_EX.pdf',height=4,width=10)
par(mfrow=c(1,1), mar=c(3, 3, 1, 1) + 0.1)
plot(tail(date,K),mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='', xaxt="n")

par(new = TRUE)
plot(tail(date,K),cumsum(apply((lL_xm1[,]),2,median))-
       cumsum((lL_static[,])), xaxt="n",
     type='l',ylim=c(-5,2500),lty=2,ylab='',xlab='',lwd=2)
abline(h=0)
lines(tail(date,K),cumsum(apply((lL_tdcc[,]),2,median))-
        cumsum((lL_static[,])),lty=2)
lines(tail(date,K),cumsum(apply((lL_dcc[,]),2,median))-
        cumsum((lL_static[,])))
lines(tail(date,K),cumsum(apply((lL_tdcch),2,median))-
        cumsum((lL_static[,])),col='gray40',lty=6,lwd=3)
legend(x=date[nn]-5,y=2500,col=c(1,1,1,'gray40','gray90'),
       lty=c(1,2,2,6,1),lwd=c(1,1,2,2,2),
       legend=c('DCC','DCC-t','AIW','DCC-HEAVY-t','avrg.stand.RV'))
atx <- seq(date[(nn+1)], date[(nn+K)], by=30)
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


print(xtable(lpbfdf,align= 'ccccccc|c',
             caption = '1-step-ahead log predictive scores ($LPS$) 
             for all individual models: Static, RiskMetrics fixed (RMf),
RiskMetrics estimated (RMe), 
Dynamic conditional correlation with Gaussian and $t$ copulas (DCC
and DCC-t), Additive Inverse Wishart (AIW) and DCC-HEAVY model with $t$ copula for 
2021/01/04 - 2021/12/31 out-of-sample period
(K = 252 observations).',
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

library(DMFP)
resDN= PMCMC_delNegro(lL_xm1,lL_tdcc,1000,c(0,2),10000,propsd = 0.5)

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




pdf('tables_and_figures/weights_EX.pdf',height=8,width=10)
par(mfrow=c(2,1), mar=c(3, 3, 1, 1) + 0.1)
plot(date[(nn+1):(nn+K)],mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(nn+1)]-60,date[(nn+K)]))
par(new = TRUE)
plot(date[(nn+1):(nn+K)],apply(ws_gew,2,median),ylim=c(0,1),
     type='l',ylab='',xlab='',xaxt="n",lwd=2,xlim=c(date[(nn+1)]-60,date[(nn+K)]))
lines(date[(nn+1):(nn+K)],apply(ws_jore1,2,median),col='gray40',lwd=2,lty=2)
lines(date[(nn+1):(nn+K)],apply(ws_jore5,2,median),lty=3)
lines(date[(nn+1):(nn+K)],apply(ws_jore10,2,median),col='gray60',lwd=2)
lines(date[(nn+1):(nn+K)],apply(ws_DN,2,median),lty=4,lwd=2)
abline(h=0.5)
axis(1, at=atx, labels=format(atx, "%Y/%m"))

legend(x=date[(T+1)]-65,y=1,col=c(1,'gray40',1,'gray60',1),
       lty=c(1,2,3,1,4),lwd=c(2,2,1,2,2),
       legend=c('Geweke','Jore1','Jore5','Jore10','DelNegro'))

plot(tail(date,K),mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(nn+1)]-60,date[(nn+K)]))
par(new = TRUE)

best.model = lL_xm1

plot(tail(date,K),cumsum(apply(log(lL_xm1[,]),2,median))-
       cumsum(apply(log(best.model[,]),2,median)),
     type='l',ylim=c(-15,15),ylab='',xlab='',xaxt="n",
     xlim=c(date[(nn+1)]-60,date[(nn+K)]))
lines(tail(date,K),cumsum(apply(log(lL_tdcc[,]),2,median))-
        cumsum(apply(log(best.model[,]),2,median)))
lines(tail(date,K),cumsum(apply(lL_gew[,],2,median))-
        cumsum(apply(log(best.model[,]),2,median)),
      col=1,lwd=2)
lines(tail(date,K),cumsum(apply(lL_jore1[,],2,median))-
        cumsum(apply(log(best.model[,]),2,median)),
      col='gray40',lwd=2,lty=2)
lines(tail(date,K),cumsum(apply(lL_jore5[,],2,median))-
        cumsum(apply(log(best.model[,]),2,median)),col='gray40',lwd=2,lty=3)
lines(tail(date,K),cumsum(apply(lL_jore10[,],2,median))-
        cumsum(apply(log(best.model[,]),2,median)),col='gray60',lwd=2)
lines(tail(date,K),cumsum(apply(lL_DN,2,median))-
        cumsum(apply(log(best.model[,]),2,median)),lty=4,lwd=2)
lines(tail(date,K),cumsum(apply((lL_equal),2,median))-
        cumsum(apply(log(best.model[,]),2,median)),col='gray60',lwd=2,lty=5)
lines(tail(date,K),cumsum(apply(log(lL_tdcch),2,median))-
        cumsum(apply(log(best.model[,]),2,median)),col='gray40',lty=6,lwd=2)
axis(1, at=atx, labels=format(atx, "%Y/%m"))
legend(x=date[(nn+1)]-65,y=15,col=c(1,1,'gray40',1,'gray60',1,'gray60','gray40'),
       lty=c(1,1,2,3,1,4,5,6),lwd=c(1,2,2,1,2,2,2,2),
       legend=c('DCC-t','Geweke','Jore1','Jore5',
                'Jore10','DelNegro','Equal','DCC-HEAVY-t'))
dev.off()

save.image('temp/res_3_EX.Rdata')



lpbf.all = c(lpbf,sum(apply(lL_gew,2,median)),
             sum(apply(lL_jore1,2,median)),
             sum(apply(lL_jore5,2,median)),
             sum(apply(lL_jore10,2,median)),
             sum(apply(lL_DN,2,median)),
             sum(apply(lL_equal,2,median)))-sum(apply(log(best.model[,]),2,median))

names(lpbf.all) = c(names(lpbf),"Gew", "Jore1", "Jore5", "Jore10","DN", "Equal")

lpbf.all

print(xtable(t(data.frame(lpbf.all)),align='cccccccccccccc',
             caption = 'Differences in the 1-step-ahead log predictive scores ($LPS$) 
             between the best fitting individual model (AIW) and the rest of the models: Static, RiskMetrics fixed (RMf),
RiskMetrics estimated (RMe), 
Dynamic conditional correlation with Gaussian and $t$ copulas (DCC
and DCC-t), DCC-HEAVY with $t$ copula and various AIW-DCC-t pools for 
2021/01/04 - 2021/12/31 out-of-sample period
(K = 252 observations).',
             label = 'table:lps_all_EX', digits = 2),
      file='tables_and_figures/lps_all_EX.tex',
      include.rownames = FALSE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,scalebox = 0.8)

