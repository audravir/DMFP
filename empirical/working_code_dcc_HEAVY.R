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
M     = 5000



propsd = 0.001


propsdnu = 1



TIMING = rep(NA,M)
t0   = Sys.time()
t1   = Sys.time()

TT   = dim(data)[1]
dm   = dim(data)[2]
bi   = min(M,10^4)


R    = array(NA,c(dm, dm, TT))
aold   <- rep(0.1,dm)
bold   <- rep(0.8,dm)
nuold  <- 30
tdata  <- qt(udata,nuold)
Rbar   <- cor(tdata)
llold  <- rep(0,TT)
LLH    <- rep(NA,M)
R[,,1] <- cor(tdata)
Pbar   = Reduce('+',Sig)/T0
A      = Outer(aold,aold)
B      = Outer(bold,bold)
iota   = rep(1,dm)
Oiota  = Outer(iota,iota)
Rtilde = (Oiota-A-B)*Rbar
is.positive.definite(Rtilde)

llold    <- rep(0,TT)
restdcch <- matrix(NA,ncol=dm*2+1,nrow=M)
accdcc   <- rep(0,bi+M)
accnu    <- rep(0,bi+M)
Rpred = vector(mode = "list", length = M)


for(t in 2:TT){
  R[,,t]   <- Rtilde+A*Sig[[t-1]]+B*R[,,t-1]
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
  
  counter = 0
  repeat{
    counter = counter+1
    bn   = rnorm(dm*2,c(aold,bold),sd=propsd*fac[1])
    anew = bn[1:dm]
    bnew = bn[(dm+1):(2*dm)]
    A    = Outer(anew,anew)
    B    = Outer(bnew,bnew)
    Rtilde = (Oiota-A-B)*Rbar
    cond1 = anew[1]>0
    cond2 = bnew[1]>0
    cond3 = (prod(eigen(Rtilde,symmetric = TRUE,only.values = TRUE)$values>0)==1)
    cond4 = (sum(abs(A+B)<1)==dm^2)
    if(cond1 && cond2 && cond3 && cond4) break
    if (counter>5){
      print(c(cond1,cond2,cond3,cond4,m))
      fac[1] = fac[1]/2
    }
  }
  
  llnew <- rep(0,TT)

  for(t in 2:TT){
    R[,,t]   <- Rtilde+A*Sig[[t-1]]+B*R[,,t-1]
    inlik    <- sum(dt(tdata[t,],df=nuold,log=TRUE))
    llnew[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nuold, log=TRUE)-
      inlik
  }

  
  if((sum(llnew)-sum(llold)+
      sum(dnorm(anew,0,sqrt(10),log=TRUE))-sum(dnorm(aold,0,sqrt(10),log=TRUE))+
      sum(dnorm(bnew,0,sqrt(10),log=TRUE))-sum(dnorm(bold,0,sqrt(10),log=TRUE)))>log(runif(1)))
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
  A     = Outer(aold,aold)
  B     = Outer(bold,bold)
  Rtilde = (Oiota-A-B)*Rbar
  
  
  for(t in 2:TT){
    R[,,t]   <- Rtilde+A*Sig[[t-1]]+B*R[,,t-1]
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
  Rtilde = (Oiota-A-B)*Rbar
  
  
  if(m>bi){
    restdcch[m-bi,] <- c(nuold,aold,bold) 
    LLH[m-bi]     <- sum(llold)
    Rpred[[m-bi]] <- Rtilde+A*Sig[[TT]]+B*R[,,TT]
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
for(i in 1:dm) {plot(b1[,i],type='l')}

par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b2[,i],type='l')}



res = list(restdcch,accnu[(bi+1):(bi+M)],accdcc[(bi+1):(bi+M)],Rpred)
names(res) = c('restdcch','accnu','accdcc','Rpred')

save(res,file='empirical/temp/results_heavy.Rdata')


