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

data = stand[1:T0,]
M=30000

# 0.001 gave accp 0.002
propsd=0.0001

# 0.1 gave accp 0.681
propsdnu = 0.5

###

t0   = Sys.time()
TT   = dim(data)[1]
dm   = dim(data)[2]
bi   = min(M,10^4)

udata = pnorm(data)*TT/(TT+1) 

Qold = array(NA,c(dm, dm, TT))
aold   <- rep(0.1,dm)
bold   <- rep(0.99,dm)
nuold  <- 20
tdata  <- qt(udata,nuold)
llold  <- rep(0,TT)

Qold[,,1] = cor(tdata)

resdcc <- matrix(NA,ncol=dm*2+1,nrow=M)
iota   = rep(1,dm)
Oiota  = Outer(iota,iota)
Sbar   = cov(tdata)
A      = Outer(aold,aold)
B      = Outer(bold,bold)
B0     = (Oiota-A-B)*Sbar

accdcc <- rep(0,bi+M)
accnu  <- rep(0,bi+M)
Vpred  = Qpred = vector(mode = "list", length = M)

for(t in 2:TT){
  Qold[,,t] <- B0+A*Outer(tdata[t-1,],tdata[t-1,])+B*Qold[,,(t-1)]
  t.ma  = Qold[,,t]
  t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
  t.R   = Outer(t.dv,t.dv)*t.ma
  inlik = sum(dt(tdata[t,],df=nuold,log=TRUE))
  llold[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm),t.R, nuold, log=TRUE)-inlik
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
    B0  = (Oiota-Outer(anew,anew)-Outer(bnew,bnew))*Sbar
    if(anew[1]>0 && bnew[1]>0 && prod(eigen.sym(B0,dm-1,vectors = FALSE)$values>0)==1 && (sum(abs(Outer(anew,anew)+Outer(bnew,bnew))<1)==dm^2)) break
  }
  
  llnew <- rep(0,TT)
  Qnew  = Qold
  A     = Outer(anew,anew)
  B     = Outer(bnew,bnew)
  
  for(t in 2:TT){
    Qnew[,,t] <- B0+A*Outer(tdata[t-1,],tdata[t-1,])+B*Qnew[,,(t-1)]
    t.ma  = Qnew[,,t]
    t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
    t.R   = Outer(t.dv,t.dv)*t.ma
    inlik = sum(dt(tdata[t,],df=nuold,log=TRUE))
    llnew[t]  <- mvnfast::dmvt(tdata[t,], rep(0,dm), t.R, nuold, log=TRUE)-inlik
  }
  
  if((sum(llnew)-sum(llold)+
      sum(dnorm(anew,0,sqrt(10),log=TRUE))-sum(dnorm(aold,0,sqrt(10),log=TRUE))+
      sum(dnorm(bnew,0,sqrt(10),log=TRUE))-sum(dnorm(bold,0,sqrt(10),log=TRUE)))>log(runif(1)))
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
  
  llnew <- rep(0,TT)
  Qnew  = Qold
  A     = Outer(aold,aold)
  B     = Outer(bold,bold)
  B0    = (Oiota-A-B)*Sbar
  
  
  for(t in 2:TT){
    Qnew[,,t] <- B0+A*Outer(tdata[t-1,],tdata[t-1,])+B*Qnew[,,(t-1)]
    
    t.ma  = Qnew[,,t]
    t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
    t.R   = Outer(t.dv,t.dv)*t.ma
    inlik = sum(dt(tdata[t,],df=nunew,log=TRUE))
    llnew[t]  <- mvnfast::dmvt(tdata[t,], rep(0,dm), t.R, nunew, log=TRUE)-inlik
      
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
  A      = Outer(aold,aold)
  B      = Outer(bold,bold)
  B0     = (Oiota-A-B)*Sbar
  
  if(m>bi){
    resdcc[m-bi,] <- c(nuold,aold,bold) 
    Qpred[[m-bi]] <- B0+A*Outer(tdata[TT,],tdata[TT,])+B*Qold[,,TT]
    Vpred[[m-bi]] <- diag(diag(Qpred[[m-bi]])^{-1/2})%*%Qpred[[m-bi]]%*%diag(diag(Qpred[[m-bi]])^{-1/2})
  }

  
  if(!m%%100){
    print(paste(round(m/(M+bi)*100),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
  }
}  
# })


mean(accnu[(bi+1):(bi+M)])
mean(accdcc[(bi+1):(bi+M)])


res = list(Vpred,Qpred,resdcc,accdcc[(bi+1):(bi+M)],accnu[(bi+1):(bi+M)])
names(res) = c('Vpred','Qpred','resdcc','accdcc','accnu')
save(res,file=paste('empirical/temp/results_vectordcc_tcop.Rdata',sep=''))