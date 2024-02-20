rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')

library(LaplacesDemon)
library(Rfast)
library(profvis)
library(future.apply)
plan(multisession, workers = 4)

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

data = Sigma[1:T0]
rm(Sigma)

# function arguments
M =25000

propsdb  = 0.0002 
propsdnu = 0.1

t0   = Sys.time()
t1   = Sys.time()
# dwish  = function(Sig,nu,S){dwishart(Sig, nu, S/nu, log=TRUE)}
# this density function is the same as in GOLOSNOY et al (2012), Eq (3)
Sig    = data
dm     = dim(Sig[[1]])[1]
TT     = length(Sig)
bi     = min(M,25000)
TIMING = rep(NA,M+bi)
resc   = matrix(NA,nrow=M+bi,ncol=dm*2+1)
LLH    = rep(NA,M+bi)
Vpred  = vector(mode = "list", length = M)
nu     = 18 # good starting value
b1     = rep(0.90,dm) # 0.9 good starting value
b2     = rep(0.40,dm) # 0.4 good starting value
Sbar   = Reduce('+',Sig)/TT
iota   = rep(1,dm)
Oiota  = Outer(iota,iota)
B1     = Outer(b1,b1) 
B2     = Outer(b2,b2)
B0     = (Oiota-B1-B2)*Sbar
is.positive.definite(B0)
llo    = lln = rep(0,TT)
accB1  = accB2 = accnu = rep(0,bi+M)
V      = Vn = list()
V[[1]] = Vn[[1]] = Sbar


for(t in 2:TT){
  V[[t]]   = B0+B1*V[[t-1]]+B2*Sig[[t-1]]
}
  
dwish.t   <- function(x,y){LaplacesDemon::dwishart(x, nu, y/nu, log=TRUE)}
# this is the best for home PC 
llo <- future_mapply(dwish.t,Sig,V)

for(m in 1:(bi+M)){
  t2=Sys.time()
  
  ##-----
  ## bs split randomly
  ##-----
  
  block1 <- sample(c(TRUE,FALSE),size=dm*2,replace = TRUE)
  block2 <- (!block1)
  
  # block 1
  repeat{
    b.prop = rnorm(dm*2,c(b1,b2),sd=propsdb)
    bn     = b.prop*block1+c(b1,b2)*block2
    
    b1n = bn[1:dm]
    b2n = bn[(dm+1):(2*dm)]
    B1  = Outer(b1n,b1n)
    B2  = Outer(b2n,b2n)
    B0  = (Oiota-B1-B2)*Sbar
    
    cond1 = b1n[1]>0
    cond2 = b2n[1]>0
    # cond3 = prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1 
    cond4 = sum(abs(B1+B2)<1)==dm^2
    if(cond1 && cond2 && cond4) {
      break
    }
  }
  
  IPD   <- rep(1,TT)
  
  for(t in 2:TT){
    Vn[[t]]   = B0+B1*Vn[[t-1]]+B2*Sig[[t-1]]
    IPD[t] =  prod(eigen(Vn[[t]],symmetric = TRUE,only.values = TRUE)$values>0)==1
  }
  
  if(sum(IPD)==TT){
    dwish.t <- function(x,y){LaplacesDemon::dwishart(x, nu, y/nu, log=TRUE)}
    lln     <- future_mapply(dwish.t,Sig,Vn)
  } else {lln=-Inf}
  
  if((sum(lln)-sum(llo)+
      sum(dnorm(b1n,0,sqrt(10),log=TRUE))-sum(dnorm(b1,0,sqrt(10),log=TRUE))+
      sum(dnorm(b2n,0,sqrt(10),log=TRUE))-sum(dnorm(b2,0,sqrt(10),log=TRUE)))>log(runif(1))){
    b1   = b1n
    b2   = b2n
    accB1[m] = 1
    llo  = lln
    V    = Vn
  }
 
  # block 2
  repeat{
    b.prop = rnorm(dm*2,c(b1,b2),sd=propsdb)
    bn     = b.prop*block2+c(b1,b2)*block1
    
    b1n = bn[1:dm]
    b2n = bn[(dm+1):(2*dm)]
    B1  = Outer(b1n,b1n)
    B2  = Outer(b2n,b2n)
    B0  = (Oiota-B1-B2)*Sbar

    cond1 = b1n[1]>0
    cond2 = b2n[1]>0
    # cond3 = prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1 
    cond4 = sum(abs(B1+B2)<1)==dm^2
    if(cond1 && cond2 && cond4) {
      break
    }
  }
  
  IPD   <- rep(1,TT)
  
  for(t in 2:TT){
    Vn[[t]]   = B0+B1*Vn[[t-1]]+B2*Sig[[t-1]]
    IPD[t] =  prod(eigen(Vn[[t]],symmetric = TRUE,only.values = TRUE)$values>0)==1
  }
  
  if(sum(IPD)==TT){
    dwish.t <- function(x,y){LaplacesDemon::dwishart(x, nu, y/nu, log=TRUE)}
    lln     <- future_mapply(dwish.t,Sig,Vn)
  } else {lln=-Inf}
  
  if((sum(lln)-sum(llo)+
      sum(dnorm(b1n,0,sqrt(10),log=TRUE))-sum(dnorm(b1,0,sqrt(10),log=TRUE))+
      sum(dnorm(b2n,0,sqrt(10),log=TRUE))-sum(dnorm(b2,0,sqrt(10),log=TRUE)))>log(runif(1))){
    b1   = b1n
    b2   = b2n
    accB2[m] = 1
    llo  = lln
    V    = Vn
  }
  B1 = Outer(b1,b1) 
  B2 = Outer(b2,b2)
  B0 = (Oiota-B1-B2)*Sbar
  
  #-----
  # nu
  #-----

  nun  = truncnorm::rtruncnorm(1,a = dm+1,mean = nu,sd = propsdnu)

  dwish.t   <- function(x,y){LaplacesDemon::dwishart(x, nun, y/nun, log=TRUE)}
  lln <- future_mapply(dwish.t,Sig,V)

  if((sum(lln)-sum(llo)+
      dexp(nun,0.01,log=TRUE)-dexp(nu,0.01,log=TRUE)+
      log(truncnorm::dtruncnorm(nun,a = dm+1,b=Inf,mean = nu,sd = propsdnu))-
      log(truncnorm::dtruncnorm(nu,a = dm+1,b=Inf,mean = nun,sd = propsdnu)))>log(runif(1))){
    llo   = lln
    accnu[m] = 1
    nu    = nun
  }

  ##-----
  ## Collect results and prediction
  ##-----
  resc[m,] = c(nu,b1,b2)
  LLH[m]   = sum(llo)
  
  if(m>bi){
    Vpred[[m-bi]] = B0+B1*V[[TT]]+B2*Sig[[TT]]
  }
  
  if(!m%%100){
    print(paste(round(m/(M+bi)*100,2),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
    print(round(c(mean(accnu[1:m]),mean(accB1[1:m]),mean(accB2[1:m])),2))
    t1   = Sys.time()
  }
  TIMING[m] = Sys.time()-t2
}


mean(accnu[(bi+1):(bi+M)])  
mean(accB1[(bi+1):(bi+M)])
mean(accB2[(bi+1):(bi+M)])


par(mfrow=c(2,2))
plot(resc[,1],type='l')
plot(LLH,type='l')
plot(TIMING,type='l')

b1=resc[,2:(dm+1)]
b2=resc[,(dm+2):(dm*2+1)]

par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b1[,i],type='l')}

par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b2[,i],type='l')}

library(corrplot)
par(mfrow=c(1,1))
corrplot(cor(resc[,-1])) 

post.size = 1000
ind       = round(seq(1,M,length=post.size))
r         = resc[(bi+1):(bi+M),]

res = list(Vpred[ind],r[ind,],
           accnu[(bi+1):(bi+M)][ind],accB1[(bi+1):(bi+M)][ind],accB2[(bi+1):(bi+M)][ind],
           LLH[ind])
  names(res) = c('Vpred','r','accnu','accB1','accB1','LLH')
  save(res,file='empirical/temp/results_caw.Rdata')


