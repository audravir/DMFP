rm(list=ls(all=TRUE))
load('data/FXdata.Rdata')

library(Rfast)
library(mvnfast)
library(profvis)

# profvis({
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

data = stand[1:T0,-15]
M    = 5000



# 0.0005 gave accp 0.001
# 0.00005 gave 0.4684
propsd = 0.0001




t0   = Sys.time()
TT   = dim(data)[1]
dm   = dim(data)[2]
bi   = min(M,10^4)



Qold = array(NA,c(dm, dm, TT))
aold <- rep(0.1,dm)
bold <- rep(0.95,dm)
llold  <- rep(0,TT)


LLH <- rep(NA,M)
Qold[,,1] = cor(data)

resdcc <- matrix(NA,ncol=dm*2,nrow=M)
iota   = rep(1,dm)
Oiota  = Outer(iota,iota)
Sbar   = cov(data)
A      = Outer(aold,aold)
B      = Outer(bold,bold)
B0     = (Oiota-A-B)*Sbar  

accdcc <- rep(0,bi+M)

Vpred  = Qpred = vector(mode='list',length = M)
  
for(t in 2:TT){
  Qold[,,t] <- B0+A*Outer(data[t-1,],data[t-1,])+B*Qold[,,(t-1)]
  t.ma  = Qold[,,t]
  t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
  t.R   = Outer(t.dv,t.dv)*t.ma
  llold[t] <- mvnfast::dmvn(data[t,], rep(0,dm), t.R, log=TRUE)
}
  

for(m in 1:(M+bi)){
  
  t1 <- Sys.time()  
  ##-----
  ## bs
  ##-----
  # 10% of the time sample from large variance 
  fac = sample(c(1,sqrt(10)),1,prob = c(0.9,0.1))
  
  repeat{
    bn   = rnorm(dm*2,c(aold,bold),sd=propsd*fac)
    anew = bn[1:dm]
    bnew = bn[(dm+1):(2*dm)]
    A    = Outer(anew,anew)
    B    = Outer(bnew,bnew)
    B0   = (Oiota-A-B)*Sbar
    if(anew[1]>0 && bnew[1]>0 && (prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1 ) && (sum(abs(A+B)<1)==dm^2)) break
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
    accdcc[m] = 1
  }
  
  A     = Outer(aold,aold)
  B     = Outer(bold,bold)
  B0    = (Oiota-A-B)*Sbar
  
  if(m>bi){
    resdcc[m-bi,] <- c(aold,bold) 
    LLH[m-bi]     <- sum(llold)
    Qpred[[m-bi]] <- B0+A*Outer(data[TT,],data[TT,])+B*Qold[,,TT]
    Vpred[[m-bi]] <- diag(diag(Qpred[[m-bi]])^{-1/2})%*%Qpred[[m-bi]]%*%diag(diag(Qpred[[m-bi]])^{-1/2})
  }
  
  
  if(!m%%100){
    print(paste(round(m/(M+bi)*100),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
    }
}  
# })


mean(accdcc[(bi+1):(bi+M)])

par(mfrow=c(1,1))
plot(LLH,type='l')

b1=resdcc[,1:(dm)]
b2=resdcc[,(dm+1):(dm*2)]

par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b1[,i],type='l')}

par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b2[,i],type='l')}

res = list(Vpred,Qpred,resdcc,accdcc[(bi+1):(bi+M)],LLH)
names(res) = c('Vpred','Qpred','resdcc','accdcc','LLH')
save(res,file=paste('empirical/temp/results_vectordcc.Rdata',sep=''))



  
