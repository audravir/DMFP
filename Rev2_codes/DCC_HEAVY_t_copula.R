rm(list=ls(all=TRUE))
library(sfsmisc)
library(truncnorm)
library(DMFP)
library(matrixcalc)
library(mixAK)
library(LaplacesDemon)

load('data/10data.Rdata')
T0  = 1990
standrets = stand
M = 25000

data  = standrets[1:T0,]
udata = pnorm(data)*T0/(T0+1) 
Sig   = Sigma[1:T0]
  
TT   = dim(data[1:T0,])[1]
dm   = dim(data[1:T0,])[2]
bi   = M



  R    = array(NA,c(dm, dm, TT))
  aold   <- 0.01
  bold   <- 0.98
  nuold  <- 10
  tdata  <- qt(udata,nuold)
  
  Pbar = Reduce('+',Sig)/T0
  
  R[,,1] <- cor(tdata)
  
  llold   <- rep(0,TT)
  restdcch <- matrix(NA,ncol=3,nrow=(bi+M))
  acctdcch <- rep(0,bi+M)

  for(t in 2:TT){
    Rbar     <- cor(tdata)
    R[,,t]   <- Rbar+aold*(Sig[[t-1]]-Pbar)+bold*(R[,,t-1]-Rbar)
    llold[t] <- mvtnorm::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nuold, log=T)-
      sum(dt(tdata[t,],df=nuold,log=TRUE))
  }
  
  for(m in 1:(M+bi)){
    parnew = rnorm(2,c(aold,bold),0.02)
    nunew  = truncnorm::rtruncnorm(1,a = 4,mean = nuold,sd = 2) # 2 for 3variate
    anew   = parnew[1]
    bnew   = parnew[2]
    llnew  = rep(0,TT)
    tdata  = qt(udata,nunew)
    Rbar   = cor(tdata)
    R[,,1] = Rbar
    
    for(t in 2:TT){
      R[,,t]   <- Rbar+anew*(Sig[[t-1]]-Pbar)+bnew*(R[,,t-1]-Rbar)
      llnew[t] <- mvtnorm::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nunew, log=TRUE)-
        sum(dt(tdata[t,],df=nunew,log=TRUE))
    }
    
    if((sum(llnew)-sum(llold)+
        dbeta(anew,3,10,log=TRUE)-dbeta(aold,3,10,log=TRUE)+
        dbeta(bnew,10,3,log=TRUE)-dbeta(bold,10,3,log=TRUE)+
        dexp(nunew,0.1,log=TRUE)-dexp(nuold,0.1,log=TRUE)+
        log(truncnorm::dtruncnorm(nunew,a = 4,mean = nuold,sd = 2))-
        log(truncnorm::dtruncnorm(nuold,a = 4,mean = nunew,sd = 2)))>log(runif(1)))
    {
      llold  = llnew
      aold   = anew
      bold   = bnew
      nuold  = nunew
      acctdcch[m] = 1
    }
    restdcch[m,] <- c(aold,bold,nuold) 
  }  
  
  
  res = list(restdcch[(bi+1):(bi+M),],acctdcch[(bi+1):(bi+M)])
  names(res) = c('restdcch','acctdcch')
  
  par(mfrow=c(1,3))
  for(i in 1:3) plot(res$restdcch[,i],type='l')
  
  mean(res$acctdcch)
  
  save(res,file='temp/10variate/results_dcch_t.Rdata')
  
  
  