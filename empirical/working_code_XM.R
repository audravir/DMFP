rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')

library(LaplacesDemon)
library(Rfast)
library(profvis)
library(future.apply)
parallel::detectCores()
plan(multisession, workers = 4)

profvis({

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
M = 1000

# 0.001 gives accp 0.516
propsdb = 0.005
propsdnu = 0.001



t0   = Sys.time()

diwish = function(Sig,nu,S){dinvwishart(Sig, nu, S, log=TRUE)}

Sig    = data
dm     = dim(Sig[[1]])[1]
TT     = length(Sig)
bi     = M
resc   = matrix(NA,nrow=M,ncol=dm*2+2)
Vpred  = vector(mode = "list", length = M)
nu     = 17
lag    = 15
b1     = b2 = rep(0.2,dm)
Sbar   = Reduce('+',Sig)/TT
iota   = rep(1,dm)
Oiota  = Outer(iota,iota)
B1     = Outer(b1,b1) 
B2     = Outer(b2,b2)
B0     = (Oiota-B1-B2)*Sbar
llo    = lln = rep(0,TT)
accB   = accnu = accl = rep(0,bi+M)
V      = Vn = list()
G1     = G2 = G2n = c(list(matrix(0,nrow=dm,ncol=dm)),Sig[-TT])
for(t in 2:TT) G2[[t]]     = Reduce('+',Sig[max(1,t-lag):(t-1)])/min(c(t-1,lag))

for(t in 1:TT){
  V[[t]]   = (B0+(b1%*%t(b1))*G1[[t]]+(b2%*%t(b2))*G2[[t]])
  # llo[t]   = diwish(Sig[[t]],nu,(nu-dm-1)*V[[t]])
}

diwish.t = function(x,y){LaplacesDemon::dinvwishart(x, nu, (nu-dm-1)*y, log=TRUE)}
llo <- future_mapply(diwish.t,Sig,V)

for(m in 1:(bi+M)){
  t1=Sys.time()
  ##-----
  ## l
  ##-----
  repeat{
    pos  = rbinom(1,1,0.5)
    lagnew = lag+sample(c(1,2,3,4),1,prob = c(6/12,3/12,2/12,1/12))*((-1)^(1-pos))
    if(lagnew>1) break
  }
  
  for(t in 2:TT) G2n[[t]]     = Reduce('+',Sig[max(1,t-lagnew):(t-1)])/min(c(t-1,lagnew))
  for(t in 1:TT){
    V[[t]]   = (B0+B1*G1[[t]]+B2*G2n[[t]])
    # lln[t]   = diwish(Sig[[t]],nu,(nu-dm-1)*V[[t]])
  }
  
  diwish.t <-  function(x,y){LaplacesDemon::dinvwishart(x, nu, (nu-dm-1)*y, log=TRUE)}
  lln      <- future_mapply(diwish.t,Sig,V)
  
  if(m%%(rbinom(1,200,0.5))!=0){
    if((sum(lln)-sum(llo))>log(runif(1))){
      llo     = lln
      lag     = lagnew
      G2      = G2n
      accl[m] = 1
    }
  } else {
    llo     = lln
    lag     = lagnew
    G2      = G2n
    accl[m] = 1
  }

  ##-----
  ## bs
  ##-----
  # 10% of the time sample from large variance to 
  fac = sample(c(1,10),1,prob = c(0.9,0.1))
  
  repeat{
    bn  = rnorm(dm*2,c(b1,b2),sd=propsdb*fac)
    b1n = bn[1:dm]
    b2n = bn[(dm+1):(2*dm)]
    B1  = Outer(b1n,b1n)
    B2  = Outer(b2n,b2n)
    B0  = (Oiota-B1-B2)*Sbar
    if(b1n[1]>0 && b2n[1]>0 && (prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1 )&& (sum(abs(B1+B2)<1)==dm^2)) break
  }
  for(t in 1:TT){
    Vn[[t]]  = (B0+B1*G1[[t]]+B2*G2[[t]])
    # lln[t]   = diwish(Sig[[t]],nu,(nu-dm-1)*Vn[[t]])
  }
  
  lln <- future_mapply(diwish.t,Sig,Vn)
  
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
    nun = rnorm(1,nu,sd=propsdnu*fac)
    if(nun>(dm+1)) break
  }
  
  diwish.t <- function(x,y){LaplacesDemon::dinvwishart(x, nun, (nun-dm-1)*y, log=TRUE)}
  lln      <- future_mapply(diwish.t,Sig,V)
  
  if((sum(lln)-sum(llo)+
      dexp(nun,1/10,log=TRUE)-dexp(nu,1/10,log=TRUE))>log(runif(1))){
    llo   = lln
    accnu[m] = 1
    nu    = nun
  }
  
  if (m>bi){
    ##-----
    ## Collect results
    ##-----
    resc[m-bi,] = c(lag,nu,b1,b2)
    
    ##-----
    ## Prediction
    ##-----
    Vpred[[m-bi]]    = B0+B1*Sig[[TT]]+B2*Reduce('+',Sig[(TT+1-lag):TT])/lag
  }

  
  if(!m%%100){
    print(paste(round(m/(M+bi)*100,2),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
    }
  
}
})


mean(accl[(bi+1):(bi+M)])
mean(accnu[(bi+1):(bi+M)])
mean(accB[(bi+1):(bi+M)])

par(mfrow=c(2,1))
plot(resc[(bi+1):(bi+M),1],type='l')
plot(resc[(bi+1):(bi+M),2],type='l')


res = list(Vpred,resc,
           accl[(bi+1):(bi+M)],
           accnu[(bi+1):(bi+M)],
           accB[(bi+1):(bi+M)])
names(res) = c('Vpred','resc','accl','accnu','accB')
save(res,file='empirical/temp/results_xm.Rdata')

