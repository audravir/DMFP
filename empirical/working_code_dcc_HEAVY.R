rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')

library(Rfast)
library(mvnfast)


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
M     = 30000



propsd = 0.03


propsdnu = 1



TIMING = rep(NA,M)
t0   = Sys.time()
t1   = Sys.time()

TT   = dim(data)[1]
dm   = dim(data)[2]
bi   = min(M,10^4)

udata = pnorm(data)*TT/(TT+1) 

R    = array(NA,c(dm, dm, TT))
aold   <- rep(0.05,dm)
bold   <- rep(0.90,dm)
nuold  <- 30
tdata  <- qt(udata,nuold)
Rbar   <- cor(tdata)
llold  <- rep(0,TT)
LLH    <- rep(NA,M)
R[,,1] <- cor(tdata)
Pbar   = Reduce('+',Sig)/T0
# A      = Outer(aold,aold)
# B      = Outer(bold,bold)
A = diag(aold)
B = diag(bold)

llold    <- rep(0,TT)
restdcch <- matrix(NA,ncol=dm*2+1,nrow=M)
accdcc   <- rep(0,bi+M)
accnu    <- rep(0,bi+M)
Rpred = vector(mode = "list", length = M)



for(t in 2:TT){
  R[,,t]   <- Rbar+A*(Sig[[t-1]]-Pbar)+B*(R[,,t-1]-Rbar)
  inlik    <- sum(dt(tdata[t,],df=nuold,log=TRUE))
  llold[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nuold, log=TRUE)-
    inlik
}



for(m in 1:(M+bi)){
  
  t2 = Sys.time()
  ##-----
  ## bs
  ##-----
  # 10% of the time sample from large variance 
  fac = sample(c(1,sqrt(10)),size=2,replace=TRUE,prob = c(0.9,0.1))
  
  # repeat{
    bn   = rnorm(dm*2,c(aold,bold),sd=propsd*fac[1])
    anew = bn[1:dm]
    bnew = bn[(dm+1):(2*dm)]
    # A    = Outer(anew,anew)
    # B    = Outer(bnew,bnew)
    A = diag(anew)
    B = diag(bnew)
  # }
  
  llnew <- rep(0,TT)

  for(t in 2:TT){
    R[,,t]   <- Rbar+A*(Sig[[t-1]]-Pbar)+B*(R[,,t-1]-Rbar)
    inlik    <- sum(dt(tdata[t,],df=nuold,log=TRUE))
    llnew[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nuold, log=TRUE)-
      inlik
  }

  
  if((sum(llnew)-sum(llold)+
      sum(dbeta(anew,2,10,log=TRUE))-sum(dbeta(aold,2,10,log=TRUE))+
      sum(dbeta(bnew,10,2,log=TRUE))-sum(dbeta(bold,10,2,log=TRUE)))>log(runif(1)))
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
    nunew = rnorm(1,nuold,propsdnu*fac[2])
    if(nunew>dm) break
  }
  
  tdata  = qt(udata,nunew)
  Rbar   <- cor(tdata)

  llnew  = rep(0,TT)
  R[,,1] = Rbar
  # A     = Outer(aold,aold)
  # B     = Outer(bold,bold)
  A = diag(aold)
  B = diag(bold)
  
  
  for(t in 2:TT){
    R[,,t]   <- Rbar+A*(Sig[[t-1]]-Pbar)+B*(R[,,t-1]-Rbar)
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
  
  
  
  
  if(m>bi){
    restdcch[m-bi,] <- c(nuold,aold,bold) 
    LLH[m-bi]     <- sum(llold)
    Rpred[[m-bi]] <- Rbar+A*(Sig[[TT]]-Pbar)+B*(R[,,TT]-Rbar)
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
b1 = restdcch[,2:(dm+1)]
b2 = restdcch[,(dm+2):(2*dm+1)]

par(mfrow=c(2,2))
plot(TIMING,type='l')
plot(LLH,type='l')
plot(nu,type='l')

par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b1[,i],type='l',ylim=c(0,1))}

par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b2[,i],type='l',ylim=c(0,1))}



res = list(restdcch,accnu[(bi+1):(bi+M)],accdcc[(bi+1):(bi+M)],Rpred)
names(res) = c('restdcch','accnu','accdcc','Rpred')

save(res,file='empirical/temp/results_heavy.Rdata')


