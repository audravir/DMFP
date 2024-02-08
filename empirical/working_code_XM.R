rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')

library(LaplacesDemon)
library(Rfast)
library(profvis)
library(future.apply)
parallel::detectCores()
plan(multisession, workers = 6)

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

data = Sigma[1:T0]
rm(Sigma)
# function arguments
M = 1000

# 0.001 gives accp 0.516
propsdb  = 0.001
propsdnu = 0.01

t0   = Sys.time()
t1   = Sys.time()

diwish = function(Sig,nu,S){dinvwishart(Sig, nu, S, log=TRUE)}

Sig    = data
dm     = dim(Sig[[1]])[1]
TT     = length(Sig)
bi     = min(M,25000)
TIMING = rep(NA,M+bi)
resc   = matrix(NA,nrow=M+bi,ncol=dm*2+2)
Vpred  = vector(mode = "list", length = M)
nu     = 25 #20 too low
lag    = 15
b1     = rep(0.43,dm) #0.43 ok starting value
b2     = rep(0.76,dm) #0.76 ok starting value
Sbar   = Reduce('+',Sig)/TT
iota   = rep(1,dm)
Oiota  = Outer(iota,iota)
B1     = Outer(b1,b1) 
B2     = Outer(b2,b2)
B0     = (Oiota-B1-B2)*Sbar
llo    = lln = rep(0,TT)
LLH    = rep(NA,M+bi)
accB1  = accnu = accl = accB2 = rep(0,bi+M)
V      = Vn = list()
G1     = G2 = G2n = c(list(matrix(0,nrow=dm,ncol=dm)),Sig[-TT])
for(t in 2:TT) {G2[[t]] = Reduce('+',Sig[max(1,t-lag):(t-1)])/min(c(t-1,lag))}

for(t in 1:TT){
  V[[t]]   = (B0+(b1%*%t(b1))*G1[[t]]+(b2%*%t(b2))*G2[[t]])
  # llo[t]   = diwish(Sig[[t]],nu,(nu-dm-1)*V[[t]])
}

diwish.t = function(x,y){LaplacesDemon::dinvwishart(x, nu, (nu-dm-1)*y, log=TRUE)}
llo <- future_mapply(diwish.t,Sig,V)



for(m in 1:(bi+M)){
  t2   = Sys.time()
  fac1 = fac2 = 1
  
  ##-----
  ## l
  ##-----
  repeat{
    pos  = rbinom(1,1,0.5)
    lagnew = lag+(-1)^(pos)
    if(lagnew>1) break
  }
  
  for(t in 2:TT) G2n[[t]] = Reduce('+',Sig[max(1,t-lagnew):(t-1)])/min(c(t-1,lagnew))
  for(t in 1:TT){
    V[[t]]   = (B0+B1*G1[[t]]+B2*G2n[[t]])
    # lln[t]   = diwish(Sig[[t]],nu,(nu-dm-1)*V[[t]])
  }
  
  diwish.t <- function(x,y){LaplacesDemon::dinvwishart(x, nu, (nu-dm-1)*y, log=TRUE)}
  lln      <- future_mapply(diwish.t,Sig,V)
  
  # if(m%%(rbinom(1,200,0.5))!=0){
  if((sum(lln)-sum(llo))>log(runif(1))){
    llo     = lln
    lag     = lagnew
    G2      = G2n
    accl[m] = 1
  }
  # } else {
  #   llo     = lln
  #   lag     = lagnew
  #   G2      = G2n
  #   accl[m] = 1
  # }
  
  ##-----
  ## bs split randomly
  ##-----
  block1 <- sample(c(TRUE,FALSE),size=dm*2,replace = TRUE)
  block2 <- (!block1)
  
  # block 1
  counter = 0
  repeat{
    counter <- counter + 1
    b.prop = rnorm(dm*2,c(b1,b2),sd=propsdb*fac1)
    bn     = b.prop*block1+c(b1,b2)*block2
    
    b1n = bn[1:dm]
    b2n = bn[(dm+1):(2*dm)]
    B1  = Outer(b1n,b1n)
    B2  = Outer(b2n,b2n)
    B0  = (Oiota-B1-B2)*Sbar
    
    cond1 = b1n[1]>0
    cond2 = b2n[1]>0
    cond3 = (prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1 )
    cond4 = (sum(abs(B1+B2)<1)==dm^2)
    if(cond1 && cond2 && cond3 && cond4) {
      break
    }
    if(counter >= 5){
      print(paste('BL1',cond1,cond2,cond3,cond4,'iter=',m,'fac=',fac1,sep=','))
      fac1=fac1/2
    }
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
    accB1[m] = 1
    llo  = lln
    V    = Vn
  }

  # block 2
  counter = 0
  repeat{
    counter <- counter + 1
    b.prop = rnorm(dm*2,c(b1,b2),sd=propsdb*fac2)
    bn     = b.prop*block2+c(b1,b2)*block1
    
    b1n = bn[1:dm]
    b2n = bn[(dm+1):(2*dm)]
    B1  = Outer(b1n,b1n)
    B2  = Outer(b2n,b2n)
    B0  = (Oiota-B1-B2)*Sbar
    
    cond1 = b1n[1]>0
    cond2 = b2n[1]>0
    cond3 = (prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1 )
    cond4 = (sum(abs(B1+B2)<1)==dm^2)
    if(cond1 && cond2 && cond3 && cond4) {
      break
    }
    if(counter >= 5){
      print(paste('BL2',cond1,cond2,cond3,cond4,'iter=',m,'fac=',fac2,sep=','))
      fac2=fac2/2
    }
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
    accB2[m] = 1
    llo  = lln
    V    = Vn
  }

  B1  = Outer(b1n,b1n)
  B2  = Outer(b2n,b2n)
  B0  = (Oiota-B1-B2)*Sbar
  
  ##-----
  ## nu
  ##-----
  repeat{
    nun = rnorm(1,nu,sd=propsdnu)
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
  
  ##-----
  ## Collect results and predict
  ##-----
  
  resc[m,] = c(lag,nu,b1,b2)
  LLH[m]   = sum(llo)
  
  if (m>bi){
    Vpred[[m-bi]]    = B0+B1*Sig[[TT]]+B2*Reduce('+',Sig[(TT+1-lag):TT])/lag
  }
  
  if(!m%%100){
    print(paste(round(m/(M+bi)*100,2),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
    print(c(mean(accl[1:m]),mean(accnu[1:m]),mean(accB1[1:m]),mean(accB2[1:m])))
    t1   = Sys.time()
  }
  TIMING[m] = Sys.time()-t2
}

mean(accl[(bi+1):(bi+M)])
mean(accnu[(bi+1):(bi+M)])
mean(accB1[(bi+1):(bi+M)])
mean(accB2[(bi+1):(bi+M)])


par(mfrow=c(2,2))
plot(TIMING,type='l')
plot(LLH,type='l')
plot(resc[,1],type='l')
plot(resc[,2],type='l')

b1=resc[,3:(dm+2)]
b2=resc[,(dm+3):(dm*2+2)]

par(mfrow=c(4,5)) 
for(i in 1:dm) {plot(b1[,i],type='l')}

par(mfrow=c(4,5)) 
for(i in 1:dm) {plot(b2[,i],type='l')}


library(corrplot)
par(mfrow=c(1,1))
corrplot(cor(resc[,-c(1,2)])) 


res = list(Vpred,resc,
           accl[(bi+1):(bi+M)],
           accnu[(bi+1):(bi+M)],
           accB1[(bi+1):(bi+M)],
           accB2[(bi+1):(bi+M)],LLH)
names(res) = c('Vpred','resc','accl','accnu','accB1','accB2','LLH')
save(res,file='empirical/temp/results_xm.Rdata')

