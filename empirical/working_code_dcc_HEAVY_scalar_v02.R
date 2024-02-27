rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')

library(Rfast)
library(mvnfast)
library(matrixcalc)

nn       = length(date)
end.date = which(zoo::as.yearmon(date)=="ene 2020")[1]-1
if(is.na(date[end.date])){end.date = which(zoo::as.yearmon(date)=="jan 2020")[1]-1}

date[end.date]

T0  = end.date  #
T0
K = nn-end.date
K

T0+K
nn


data  = stand[1:T0,]
Sig   = Sigma[1:T0]
M     = 1000



propsd   = 0.002

t0   = Sys.time()
t1   = Sys.time()

TT   = dim(data)[1]
dm   = dim(data)[2]

R    = array(NA,c(dm, dm, TT))
pold   <- c(0.1,0.8)

Rbar   <- cor(data)
llold  <- rep(0,TT)
bi     = min(M,50000)
TIMING = rep(NA,M+bi)
LLH    <- rep(NA,M+bi)
R[,,1] <- Rbar
Rnew   = R
Pbar   = Reduce('+',Sig)/TT
Rtilde = (1-pold[2])*Rbar-pold[1]*Pbar
IPD    = rep(TRUE,TT)
is.positive.definite(Rtilde)

res <- matrix(NA,ncol=2,nrow=M+bi)
accdcc   = rep(0,bi+M)
Rpred    = vector(mode = "list", length = M)


for(t in 2:TT){
  R[,,t]   <- Rtilde+pold[1]*Sig[[t-1]]+pold[2]*R[,,t-1]
  llold[t] <- mvnfast::dmvn(data[t,], rep(0,dm), R[,,t],log=TRUE)
}

for(m in 1:(M+bi)){
  t2 = Sys.time()
  
  repeat{
    pnew  = rnorm(2,pold,sd=propsd)
    if(all(pnew>0)) break
  }
  
  llnew = rep(0,TT)
  for(t in 2:TT){
    Rnew[,,t]   <- (1-pnew[2])*Rbar-pnew[1]*Pbar+pnew[1]*Sig[[t-1]]+pnew[2]*Rnew[,,t-1]
    IPD[t] = prod(eigen(Rnew[,,t],symmetric = TRUE,only.values = TRUE)$values>0)==1 
  }
  
  if (all(IPD==TRUE)){
    for(t in 2:TT){
      llnew[t] = mvnfast::dmvn(data[t,], rep(0,dm), Rnew[,,t],log=TRUE)
    }
  } else {llnew = -Inf}
  
  if((sum(llnew)-sum(llold)+
      dbeta(pnew[1],3,10,log=TRUE)-dbeta(pold[1],3,10,log=TRUE)+
      dbeta(pnew[2],10,3,log=TRUE)-dbeta(pold[2],10,3,log=TRUE))>log(runif(1)))
  {
    llold  = llnew
    pold   = pnew
    R      = Rnew
    accdcc[m] = 1
  }
  
  ##-----
  ## Collect results and prediction
  ##-----
  
  res[m,] <- pold
  LLH[m]  <- sum(llold)
  
  if(m>bi){
    Rpred[[m-bi]] <- (1-pold[2])*Rbar-pold[1]*Pbar+pold[1]*Sig[[TT]]+pold[2]*R[,,TT]
  }
  
  if(!m%%100){
    print(paste(round(m/(M+bi)*100,2),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
    print(round(c(mean(accdcc[1:m])),2))
    t1   = Sys.time()
  }
  TIMING[m] = Sys.time()-t2
}

mean(accdcc[(bi+1):(bi+M)])

par(mfrow=c(2,3))
plot(TIMING,type='l')
plot(tail(LLH,M),type='l')
plot(tail(res[,1],M),type='l')
plot(tail(res[,2],M),type='l')

library(corrplot)
par(mfrow=c(1,1))
corrplot(cor(res)) 

apply(tail(res,M),2,median)


post.size = 5000
ind       = round(seq(1,M,length=post.size))
r         = res[(bi+1):(bi+M),]

res = list(Rpred[ind],r[ind,],
           accdcc[(bi+1):(bi+M)][ind],LLH[(bi+1):(bi+M)][ind])
names(res) = c('Rpred','r','accdcc','LLH')

save(res,file='empirical/temp/results_heavy_normal.Rdata')

