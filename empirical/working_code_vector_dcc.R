rm(list=ls(all=TRUE))
load('data/FXdata.Rdata')

library(Rfast)
library(mvnfast)
library(profvis)

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

data = stand[1:T0,]
M    = 10000

# 0.0001 gives 25%
propsd = 0.0001

t0   = t1 = Sys.time()
TT   = dim(data)[1]
dm   = dim(data)[2]
bi   = min(M,10^4)
TIMING = rep(NA,M+bi)


Qold = array(NA,c(dm, dm, TT))
aold <- rep(0.12,dm) # ok starting value
bold <- rep(0.99,dm) # ok starting value
llold  <- rep(0,TT)


LLH <- rep(NA,M+bi)
Qold[,,1] = cor(data)

resdcc <- matrix(NA,ncol=dm*2,nrow=M+bi)
iota   = rep(1,dm)
Oiota  = Outer(iota,iota)
Sbar   = cov(data)
A      = Outer(aold,aold)
B      = Outer(bold,bold)
B0     = (Oiota-A-B)*Sbar  

accdcc1 = accdcc2= rep(0,bi+M)

Vpred  = Qpred = vector(mode='list',length = M)
  
for(t in 2:TT){
  Qold[,,t] <- B0+A*Outer(data[t-1,],data[t-1,])+B*Qold[,,(t-1)]
  t.ma  = Qold[,,t]
  t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
  t.R   = Outer(t.dv,t.dv)*t.ma
  llold[t] <- mvnfast::dmvn(data[t,], rep(0,dm), t.R, log=TRUE)
}
  
for(m in 1:(M+bi)){
  t2 <- Sys.time()  
  fac1=fac2=1 # in case the sampler is stuck, we can reduce the variance
  
  ##-----
  ## bs split randomly
  ##-----
  
  block1 <- sample(c(TRUE,FALSE),size=dm*2,replace = TRUE)
  block2 <- (!block1)
  
  ##-----
  ## bs
  ##-----
  
  # block 1
  counter = 0
  repeat{
    counter = counter +1
    b.prop = rnorm(dm*2,c(aold,bold),sd=propsd*fac1)
    bn     = b.prop*block1+c(aold,bold)*block2

    anew = bn[1:dm]
    bnew = bn[(dm+1):(2*dm)]
    A    = Outer(anew,anew)
    B    = Outer(bnew,bnew)
    B0   = (Oiota-A-B)*Sbar
    
    cond1 = anew[1]>0
    cond2 = bnew[1]>0
    cond3 = prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1 
    cond4 = (sum(abs(A+B)<1)==dm^2)
    if(cond1 && cond2 && cond3 && cond4) {
      break
    }
    if(counter >= 5){
      print(paste('BL1',cond1,cond2,cond3,cond4,'iter=',m,'fac=',fac1,sep=','))
      fac1=fac1/sqrt(10)
    }
  }
  
  llnew <- rep(0,TT)
  Qnew  = Qold
  
  for(t in 2:TT){
    Qnew[,,t] <- B0+A*Outer(data[t-1,],data[t-1,])+B*Qnew[,,(t-1)]
    t.ma  = Qnew[,,t]
    t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
    t.R   = Outer(t.dv,t.dv)*t.ma
    llnew[t]  <- mvnfast::dmvn(data[t,], rep(0,dm), t.R, log=TRUE)
  }
  
  if((sum(llnew)-sum(llold)+
      sum(dnorm(anew,0,sqrt(10),log=T))-sum(dnorm(aold,0,sqrt(10),log=T))+
      sum(dnorm(bnew,0,sqrt(10),log=T))-sum(dnorm(bold,0,sqrt(10),log=T)))>log(runif(1)))
  {
    llold  = llnew
    aold   = anew
    bold   = bnew
    Qold   = Qnew
    accdcc1[m] = 1
  }
  
  # block 2
  counter = 0
  repeat{
    counter = counter+1
    b.prop = rnorm(dm*2,c(aold,bold),sd=propsd*fac2)
    bn     = b.prop*block2+c(aold,bold)*block1
    
    anew = bn[1:dm]
    bnew = bn[(dm+1):(2*dm)]
    A    = Outer(anew,anew)
    B    = Outer(bnew,bnew)
    B0   = (Oiota-A-B)*Sbar
    
    cond1 = anew[1]>0
    cond2 = bnew[1]>0
    cond3 = prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1 
    cond4 = sum(abs(A+B)<1)==dm^2
    if(cond1 && cond2 && cond3 && cond4) {
      break
    }
    if(counter >= 5){
      print(paste('BL2',cond1,cond2,cond3,cond4,'iter=',m,'fac=',fac2,sep=','))
      fac2=fac2/sqrt(10)
    }
  }
  
  llnew <- rep(0,TT)
  Qnew  = Qold
  
  for(t in 2:TT){
    Qnew[,,t] <- B0+A*Outer(data[t-1,],data[t-1,])+B*Qnew[,,(t-1)]
    t.ma  = Qnew[,,t]
    t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
    t.R   = Outer(t.dv,t.dv)*t.ma
    llnew[t]  <- mvnfast::dmvn(data[t,], rep(0,dm), t.R, log=TRUE)
  }
  
  if((sum(llnew)-sum(llold)+
      sum(dnorm(anew,0,sqrt(10),log=T))-sum(dnorm(aold,0,sqrt(10),log=T))+
      sum(dnorm(bnew,0,sqrt(10),log=T))-sum(dnorm(bold,0,sqrt(10),log=T)))>log(runif(1)))
  {
    llold  = llnew
    aold   = anew
    bold   = bnew
    Qold   = Qnew
    accdcc2[m] = 1
  }
  
  ##-----
  ## Collect results and prediction
  ##-----
  
  resdcc[m,] <- c(aold,bold) 
  LLH[m]     <- sum(llold)
  
  if(m>bi){
    A     = Outer(aold,aold)
    B     = Outer(bold,bold)
    B0    = (Oiota-A-B)*Sbar
    Qpred[[m-bi]] <- B0+A*Outer(data[TT,],data[TT,])+B*Qold[,,TT]
    Vpred[[m-bi]] <- diag(diag(Qpred[[m-bi]])^{-1/2})%*%Qpred[[m-bi]]%*%diag(diag(Qpred[[m-bi]])^{-1/2})
  }
  
  
  if(!m%%100){
    print(paste(round(m/(M+bi)*100,2),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
    t1 = Sys.time()
    print(round(c(mean(accdcc1[1:m]),mean(accdcc2[1:m])),2))
  }
  TIMING[m] = Sys.time()-t2
}  

mean(accdcc1[(bi+1):(bi+M)])
mean(accdcc2[(bi+1):(bi+M)])

par(mfrow=c(2,1))
plot(LLH,type='l')
plot(TIMING,type='l')

b1=resdcc[,1:(dm)]
b2=resdcc[,(dm+1):(dm*2)]

par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b1[,i],type='l')}

par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b2[,i],type='l')}

res = list(Vpred,Qpred,resdcc,accdcc1[(bi+1):(bi+M)],accdcc2[(bi+1):(bi+M)],LLH)
names(res) = c('Vpred','Qpred','resdcc','accdcc1','accdcc2','LLH')
save(res,file=paste('empirical/temp/results_vectordcc.Rdata',sep=''))



  
