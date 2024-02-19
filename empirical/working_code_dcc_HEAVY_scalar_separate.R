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
M     = 50000



propsd = 0.0001


propsdnu = 0.1

t0   = Sys.time()
t1   = Sys.time()

TT   = dim(data)[1]
dm   = dim(data)[2]

R    = array(NA,c(dm, dm, TT))
aold   <- 0.2
bold   <- 0.5 
nuold  <- 16
tdata  <- qt(udata,nuold)
Rbar   <- cor(tdata)
llold  <- rep(0,TT)
bi     = min(M,25000)
TIMING = rep(NA,M+bi)
LLH    <- rep(NA,M+bi)
R[,,1] <- cor(tdata)
Pbar   = Reduce('+',Sig)/T0
Rtilde = (1-bold)*Rbar-aold*Pbar
is.positive.definite(Rtilde)

llold    <- rep(0,TT)
restdcch <- matrix(NA,ncol=3,nrow=M+bi)
accdcc = rep(0,bi+M)
accnu  <- rep(0,bi+M)
Rpred = vector(mode = "list", length = M)

for(t in 2:TT){
  R[,,t]   <- Rtilde+aold*Sig[[t-1]]+bold*R[,,t-1]
  inlik    <- sum(dt(tdata[t,],df=nuold,log=TRUE))
  llold[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nuold, log=TRUE)-
    inlik
}

for(m in 1:(M+bi)){
  t2 = Sys.time()
 
  ##-----
  ## bs
  ##-----

  repeat{
    bn = rnorm(2,c(aold,bold),sd=propsd)
  
    anew = bn[1]
    bnew = bn[2]
    Rtilde = (1-bnew)*Rbar-anew*Pbar
    cond1 = anew>0
    cond2 = bnew>0
    cond3 = (prod(eigen(Rtilde,symmetric = TRUE,only.values = TRUE)$values>0)==1)
    if(cond1 && cond2 && cond3) break
  }
  
  llnew <- rep(0,TT)
  
  for(t in 2:TT){
    R[,,t]   <- Rtilde+anew*Sig[[t-1]]+bnew*R[,,t-1]
    inlik    <- sum(dt(tdata[t,],df=nuold,log=TRUE))
    llnew[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nuold, log=TRUE)-
      inlik
  }
  
  if((sum(llnew)-sum(llold)+
      dbeta(anew,3,10,log=TRUE)-dbeta(aold,3,10,log=TRUE)+
      dbeta(bnew,10,3,log=TRUE)-dbeta(bold,10,3,log=TRUE))>log(runif(1)))
  {
    llold  = llnew
    aold   = anew
    bold   = bnew
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
  Rbar   <- cor(tdata)
  
  llnew  = rep(0,TT)
  R[,,1] = Rbar
  Rtilde = (1-bold)*Rbar-aold*Pbar
  
  for(t in 2:TT){
    R[,,t]   <- Rtilde+aold*Sig[[t-1]]+bold*R[,,t-1]
    inlik    <- sum(dt(tdata[t,],df=nunew,log=TRUE))
    llnew[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nunew, log=TRUE)-
      inlik
  }
  
  if((sum(llnew)-sum(llold)+
      dexp(nunew,0.1,log=TRUE)-dexp(nuold,0.1,log=TRUE))>log(runif(1)))
  {
    llold  = llnew
    nuold  = nunew
    accnu[m] = 1
  }
  
  
  tdata  = qt(udata,nuold)
  Rbar   <- cor(tdata)
  Rtilde = (1-bold)*Rbar-aold*Pbar
  
  ##-----
  ## Collect results and prediction
  ##-----
  
  restdcch[m,] <- c(nuold,aold,bold) 
  LLH[m]       <- sum(llold)
  
  if(m>bi){
    Rpred[[m-bi]] <- Rtilde+aold*Sig[[TT]]+bold*R[,,TT]
  }
  
  
  if(!m%%100){
    print(paste(round(m/(M+bi)*100,2),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
    print(round(c(mean(accnu[1:m]),mean(accdcc[1:m])),2))
    t1   = Sys.time()
  }
  TIMING[m] = Sys.time()-t2
}  


mean(accnu[(bi+1):(bi+M)])
mean(accdcc[(bi+1):(bi+M)])


nu = restdcch[,1]
b1 = restdcch[,2]
b2 = restdcch[,3]

par(mfrow=c(2,3))
plot(TIMING,type='l')
plot(tail(LLH,M),type='l')
plot(restdcch[,1],type='l')
plot(restdcch[,2],type='l')
plot(restdcch[,3],type='l')


library(corrplot)
par(mfrow=c(1,1))
corrplot(cor(restdcch)) 


post.size = 5000
ind       = round(seq(1,M,length=post.size))
r         = restdcch[(bi+1):(bi+M),]

res = list(Rpred[ind],r[ind,],accnu[(bi+1):(bi+M)][ind],
           accdcc[(bi+1):(bi+M)][ind],LLH[ind])
names(res) = c('Rpred','r','accnu','accdcc','LLH')

save(res,file='empirical/temp/results_heavy_scalar_separate.Rdata')




