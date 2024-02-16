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



propsd = 0.0001
# propsd = 0.0001, 0.01
# 0.001 gave accp 0.85

# propsdnu = 0.1, 0.01
propsdnu = 0.01

t0   = Sys.time()
t1   = Sys.time()

TT   = dim(data[1:T0,])[1]
dm   = dim(data[1:T0,])[2]

R    = array(NA,c(dm, dm, TT))
aold   <- 0.20
bold   <- 0.50
nuold  <- 16
tdata  <- qt(udata,nuold)
Rbar   <- cor(tdata)
llold  <- rep(0,TT)
bi     <- min(M,25000)
TIMING = rep(NA,M+bi)
LLH    <- rep(NA,M+bi)
R[,,1] <- cor(tdata)
Pbar = Reduce('+',Sig)/T0









restdcch <- matrix(NA,ncol=3,nrow=(bi+M))
acctdcch <- rep(0,bi+M)
Rpred    = vector(mode = "list", length = M)

for(t in 2:TT){
  R[,,t]   <- Rbar+aold*(Sig[[t-1]]-Pbar)+bold*(R[,,t-1]-Rbar)
  llold[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nuold, log=TRUE)-
    sum(dt(tdata[t,],df=nuold,log=TRUE))
}




for(m in 1:(M+bi)){
  t2 = Sys.time()
  
  parnew = rnorm(2,c(aold,bold),propsd)
  nunew  = truncnorm::rtruncnorm(1,a = dm+1,mean = nuold,sd = propsdnu) # 2 for 3variate
  anew   = parnew[1]
  bnew   = parnew[2]
  llnew  = rep(0,TT)
  tdata  = qt(udata,nunew)
  Rbar   = cor(tdata)
  R[,,1] = Rbar
  
  IPD    = rep(1,TT)
  
  for(t in 2:TT){
    R[,,t] <- Rbar+anew*(Sig[[t-1]]-Pbar)+bnew*(R[,,t-1]-Rbar)
    IPD[t] =  prod(eigen(R[,,t],symmetric = TRUE,only.values = TRUE)$values>0)==1  
  }
  
  if(sum(IPD)==TT){
    for(t in 2:TT){
      llnew[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nuold, log=TRUE)-
        sum(dt(tdata[t,],df=nunew,log=TRUE))
    } 
  } else {llnew=-Inf;print(paste('m=',m,',notPD'))}
  
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
  
  tdata  = qt(udata,nuold)
  Rbar   = cor(tdata)
  
  restdcch[m,] <- c(aold,bold,nuold) 
  LLH[m]       <- sum(llold)
  
  if(m>bi){
    Rpred[[m-bi]] <- Rbar+aold*(Sig[[TT]]-Pbar)+bold*(R[,,TT]-Rbar)
  }
  
  if(!m%%100){
    print(paste(round(m/(M+bi)*100,2),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
    print(round(c(mean(acctdcch[1:m])),2))
    t1   = Sys.time()
  }
  TIMING[m] = Sys.time()-t2
}  

mean(acctdcch[(bi+1):(bi+M)])

par(mfrow=c(2,3))
plot(TIMING,type='l')
plot(LLH,type='l')
plot(restdcch[,1],type='l')
plot(restdcch[,2],type='l')
plot(restdcch[,3],type='l')


post.size = 1000
ind       = round(seq(1,M,length=post.size))
r         = restdcch[(bi+1):(bi+M),]

res = list(Rpred[ind],r[ind,],acctdcch[(bi+1):(bi+M)][ind],LLH[ind])
names(res) = c('Rpred','r','acc','LLH')

save(res,file='empirical/temp/results_heavy_scalar_previous.Rdata')

