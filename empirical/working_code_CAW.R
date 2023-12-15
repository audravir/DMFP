rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')

library(LaplacesDemon)
library(Rfast)
library(profvis)
library(future.apply)
plan(multisession, workers = 4)

# profvis({
nn       = length(date)
end.date = which(zoo::as.yearmon(date)=="ene 2021")[1]-1
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

# function arguments
M = 10000

# 0.001 give accp of 0.004
# 0.0001 0.04
propsdb = 0.00005 

# 0.001 gave accp of 0.5126
propsdnu = 0.001


t0   = Sys.time()
# dwish  = function(Sig,nu,S){dwishart(Sig, nu, S/nu, log=TRUE)}
# this density function is the same as in GOLOSNOY et al (2012), Eq (3)
Sig    = data
dm     = dim(Sig[[1]])[1]
TT     = length(Sig)
bi     = min(M,10^4)
resc   = matrix(NA,nrow=M,ncol=dm*2+1)
Vpred  = vector(mode = "list", length = M)
nu     = 17
b1     = rep(0.95,dm)
b2     = rep(0.3,dm)
Sbar   = Reduce('+',Sig)/TT
iota   = rep(1,dm)
Oiota  = Outer(iota,iota)
B1     = Outer(b1,b1) 
B2     = Outer(b2,b2)
B0     = (Oiota-B1-B2)*Sbar
llo    = lln = rep(0,TT)
accB   = accnu = rep(0,bi+M)
V      = Vn = list()
V[[1]] = Vn[[1]] = Sbar


for(t in 2:TT){
  V[[t]]   = B0+B1*V[[t-1]]+B2*Sig[[t-1]]
  # llo[t]   = dwish(Sig[[t]],nu,V[[t]])
}
  
dwish.t   <- function(x,y){LaplacesDemon::dwishart(x, nu, y/nu, log=TRUE)}

# this is the best for home PC 
llo <- future_mapply(dwish.t,Sig,V)

for(m in 1:(bi+M)){
  t1=Sys.time()
  ##-----
  ## bs
  ##-----
  repeat{
    bn  = rnorm(dm*2,c(b1,b2),sd=propsdb)
    b1n = bn[1:dm]
    b2n = bn[(dm+1):(2*dm)]
    B1  = Outer(b1n,b1n)
    B2  = Outer(b2n,b2n)
    B0  = (Oiota-B1-B2)*Sbar
    if(b1n[1]>0 && b2n[1]>0 && (prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1 ) && (sum(abs(B1+B2)<1)==dm^2)) break
  }
  for(t in 2:TT){
    Vn[[t]]   = B0+B1*Vn[[t-1]]+B2*Sig[[t-1]]
    # lln[t]   = dwish(Sig[[t]],nu,Vn[[t]])
  }
  dwish.t <- function(x,y){LaplacesDemon::dwishart(x, nu, y/nu, log=TRUE)}
  lln     <- future_mapply(dwish.t,Sig,Vn)
  
  
  if((sum(lln)-sum(llo)+
      sum(dnorm(b1n,0,sqrt(10),log=TRUE))-sum(dnorm(b1,0,sqrt(10),log=TRUE))+
      sum(dnorm(b2n,0,sqrt(10),log=TRUE))-sum(dnorm(b2,0,sqrt(10),log=TRUE)))>log(runif(1))){
    b1   = b1n
    b2   = b2n
    accB[m] = 1
    llo  = lln
    V    = Vn
  }
  B1 = Outer(b1,b1) 
  B2 = Outer(b2,b2)
  B0 = (Oiota-B1-B2)*Sbar
    
  ##-----
  ## nu
  ##-----
  repeat{
    nun = rnorm(1,nu,sd=propsdnu)
    if(nun>(dm)) break
  }
    
  dwish.t   <- function(x,y){LaplacesDemon::dwishart(x, nun, y/nun, log=TRUE)}
  lln <- future_mapply(dwish.t,Sig,V)

  if((sum(lln)-sum(llo)+
      dexp(nun,1/10,log=TRUE)-dexp(nu,1/10,log=TRUE))>log(runif(1))){
    llo   = lln
    accnu[m] = 1
    nu    = nun
  }
    
  if(m>bi){
    ##-----
    ## Collect results
    ##-----
    resc[m-bi,] = c(nu,b1,b2)
    
    ##-----
    ## Prediction
    ##-----
    Vpred[[m-bi]]    = B0+B1*V[[TT]]+B2*Sig[[TT]]
  }
  
  
  if(!m%%100){
    print(paste(round(m/(M+bi)*100,2),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
  }
}
 
# })

mean(accnu[(bi+1):(bi+M)])  
mean(accB[(bi+1):(bi+M)])



nu=resc[,1]
par(mfrow=c(1,1))
plot(nu,type='l')

b1=resc[,2:(dm+1)]
b2=resc[,(dm+2):(dm*2+1)]

par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b1[,i],type='l')}

par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b2[,i],type='l')}



res = list(Vpred,resc,accnu[(bi+1):(bi+M)],accB[(bi+1):(bi+M)])
  names(res) = c('Vpred','resc','accnu','accB')
  save(res,file='empirical/temp/results_caw.Rdata')


