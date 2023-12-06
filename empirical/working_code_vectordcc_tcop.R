rm(list=ls(all=TRUE))
load('data/FXdata.Rdata')
library(matrixcalc)
library(mvnfast)

library(profvis)

profvis({
nn       = length(date)
end.date = which(zoo::as.yearmon(date)=="ene 2020")[1]-1
if(is.na(date[end.date])){end.date = which(zoo::as.yearmon(date)=="jan 2020")[1]-1}

date[end.date]

T0  = end.date  #
T0 # for estimation
K = nn-end.date
K # for oos evaluation

T0+K
nn
# the same

T0=100
data = stand[1:T0,1:5]
M=1000
propsd=0.001
propsdnu = 0.01

###

t0   = Sys.time()
TT   = dim(data)[1]
dm   = dim(data)[2]
bi   = min(M,10^4)

udata = pnorm(data)*TT/(TT+1) 

Qold = array(NA,c(dm, dm, TT))
R    = array(NA,c(dm, dm, TT))
aold   <- rep(0.1,dm)
bold   <- rep(0.99,dm)
nuold  <- 20
tdata  <- qt(udata,nuold)
llold  <- rep(0,TT)

Qold[,,1] = cor(tdata)
R[,,1] <- diag(diag(Qold[,,1])^{-1/2})%*%Qold[,,1]%*%diag(diag(Qold[,,1])^{-1/2})

resdcc <- matrix(NA,ncol=dm*2+1,nrow=(bi+M))
iota   = rep(1,dm)
Sbar   = cov(tdata)
B0     = (iota%*%t(iota)-aold%*%t(aold)-bold%*%t(bold))*Sbar

accdcc <- rep(0,bi+M)
accnu  <- rep(0,bi+M)
Vpred  = Qpred = list()

for(t in 2:TT){
  Qold[,,t]   <- B0+(aold%*%t(aold))*(tdata[t-1,]%*%t(tdata[t-1,]))+(bold%*%t(bold))*Qold[,,(t-1)]
  R[,,t]   <- diag(diag(Qold[,,t])^{-1/2})%*%Qold[,,t]%*%diag(diag(Qold[,,t])^{-1/2})
  inlik=sum(dt(tdata[t,],df=nuold,log=TRUE))
  llold[t] <- dmvt(tdata[t,], rep(0,dm), R[,,t], df=nuold, log=TRUE)-inlik
}

for(m in 1:(M+bi)){
  
  t1 = Sys.time()
  ##-----
  ## bs
  ##-----
  repeat{
    bn  = rnorm(dm*2,c(aold,bold),sd=propsd)
    anew = bn[1:dm]
    bnew = bn[(dm+1):(2*dm)]
    B0  = (iota%*%t(iota)-anew%*%t(anew)-bnew%*%t(bnew))*Sbar
    if(anew[1]>0 && bnew[1]>0 && is.positive.definite(B0) && (sum(abs(anew%*%t(anew)+bnew%*%t(bnew))<1)==dm^2)) break
  }
  
  llnew <- rep(0,TT)
  Qnew  = Qold
  
  for(t in 2:TT){
    Qnew[,,t] <- B0+(anew%*%t(anew))*(tdata[t-1,]%*%t(tdata[t-1,]))+(bnew%*%t(bnew))*Qnew[,,(t-1)]
    R[,,t]    <- diag(diag(Qnew[,,t])^{-1/2})%*%Qnew[,,t]%*%diag(diag(Qnew[,,t])^{-1/2})
    inlik = sum(dt(tdata[t,],df=nuold,log=TRUE))
    llnew[t]  <- dmvt(tdata[t,], rep(0,dm), R[,,t], df=nuold, log=TRUE)-inlik
  }
  
  if((sum(llnew)-sum(llold)+
      sum(dnorm(anew,0,sqrt(10),log=T))-sum(dnorm(aold,0,sqrt(10),log=T))+
      sum(dnorm(bnew,0,sqrt(10),log=T))-sum(dnorm(bold,0,sqrt(10),log=T)))>log(runif(1)))
  {
    llold  = llnew
    aold   = anew
    bold   = bnew
    Qold   = Qnew
    accdcc[m] = 1
  }
  
  
  ##-----
  ## nu
  ##-----
  repeat{
    nunew = rnorm(1,nuold,propsdnu)
    if(nunew>dm) break
  }
  
  tdata  = qt(udata,nunew)
  Sbar   = cov(tdata)
  B0     = (iota%*%t(iota)-aold%*%t(aold)-bold%*%t(bold))*Sbar
  
  llnew <- rep(0,TT)
  Qnew  = Qold
  
  for(t in 2:TT){
    Qnew[,,t] <- B0+(aold%*%t(aold))*(tdata[t-1,]%*%t(tdata[t-1,]))+(bold%*%t(bold))*Qnew[,,(t-1)]
    R[,,t]    <- diag(diag(Qnew[,,t])^{-1/2})%*%Qnew[,,t]%*%diag(diag(Qnew[,,t])^{-1/2})
    inlik=sum(dt(tdata[t,],df=nunew,log=TRUE))
    llnew[t]  <- dmvt(tdata[t,], rep(0,dm), R[,,t], df=nunew, log=TRUE)-inlik
      
  }
  
  if((sum(llnew)-sum(llold)+
      dexp(nunew,0.1,log=TRUE)-dexp(nuold,0.1,log=TRUE))>log(runif(1)))
  {
    llold  = llnew
    nuold  = nunew
    Qold   = Qnew
    accnu[m] = 1
  }
  
  tdata  = qt(udata,nuold)
  Sbar   = cov(tdata)
  B0     = (iota%*%t(iota)-aold%*%t(aold)-bold%*%t(bold))*Sbar
  
  resdcc[m,] <- c(nuold,aold,bold) 
  Qpred[[m]] <- B0+aold*(tdata[TT,]%*%t(tdata[TT,]))+bold*Qold[,,TT]
  Vpred[[m]] <- diag(diag(Qpred[[m]])^{-1/2})%*%Qpred[[m]]%*%diag(diag(Qpred[[m]])^{-1/2})
  
  if(!m%%100){
    print(paste(round(m/(M+bi)*100),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
  }
}  
})

# res = list(Vpred[(bi+1):(bi+M)],Qpred[(bi+1):(bi+M)],resdcc[(bi+1):(bi+M),],
#            accdcc[(bi+1):(bi+M)],accnu[(bi+1):(bi+M)])
# names(res) = c('Vpred','Qpred','resdcc','accdcc','accnu')
# save(res,file=paste('empirical/temp/results_vectordcc_tcop.Rdata',sep=''))