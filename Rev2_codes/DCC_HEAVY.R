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
# function(standrets,Sigma,M)

data  = standrets[1:T0,]
Sig   = Sigma[1:T0]

TT   = dim(data[1:T0,])[1]
dm   = dim(data[1:T0,])[2]
bi   = M

Rbar = cor(data)
Pbar = Reduce('+',Sig)/T0

R    = array(NA,c(dm, dm, TT))
aold   <- 0.01
bold   <- 0.98
R[,,1] <- Rbar

llold   <- rep(0,TT)
resdcch <- matrix(NA,ncol=2,nrow=(bi+M))
accdcch <- rep(0,bi+M)


for(t in 2:TT){
  R[,,t] <- Rbar+aold*(Sig[[t-1]]-Pbar)+bold*(R[,,t-1]-Rbar)
  llold[t] <- mvtnorm::dmvnorm(data[t,], rep(0,dm), R[,,t], log=T)
}

for(m in 1:(M+bi)){
  parnew = rnorm(2,c(aold,bold),0.02) #0.01 is too small
  anew   = parnew[1]
  bnew   = parnew[2]
  llnew  = rep(0,TT)

  for(t in 2:TT){
    R[,,t] <- Rbar+anew*(Sig[[t-1]]-Pbar)+bnew*(R[,,t-1]-Rbar)
    llnew[t] <- mvtnorm::dmvnorm(data[t,], rep(0,dm), R[,,t], log=T)
  }
  
  if((sum(llnew)-sum(llold)+
      dbeta(anew,3,10,log=TRUE)-dbeta(aold,3,10,log=TRUE)+
      dbeta(bnew,10,3,log=TRUE)-dbeta(bold,10,3,log=TRUE))>log(runif(1)))
  {
    llold  = llnew
    aold   = anew
    bold   = bnew
    accdcch[m] = 1
  }
  resdcch[m,] <- c(aold,bold) 
}  

res = list(resdcch[(bi+1):(bi+M),],accdcch[(bi+1):(bi+M)])
names(res) = c('resdcch','accdcch')
mean(accdcch[(bi+1):(bi+M)])
plot(resdcch[,1],type='l')
plot(resdcch[,2],type='l')

save(res,file='temp/10variate/results_dcch.Rdata')
