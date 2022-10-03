#' @export
xm1 = function(Sigma,M){
  library(matrixcalc)
  library(mixAK)
  library(LaplacesDemon)
  #library(countreg)
  diwish = function(Sig,nu,S){dinvwishart(Sig, nu, S, log=TRUE)}
  riwish = function(nu,S) rinvwishart(nu, S)
  
  Sig    = Sigma[1:T0]
  dm     = dim(Sig[[1]])[1]
  bi     = M
  resc   = matrix(NA,nrow=bi+M,ncol=dm*2+2)
  Vpred  = list()
  nu     = 5
  lag    = 5
  b1     = b2 = rep(0.2,dm)
  Sbar   = Reduce('+',Sig)/T0
  iota   = rep(1,dm)
  B0     = (iota%*%t(iota)-b1%*%t(b1)-b2%*%t(b2))*Sbar
  llo    = lln = rep(0,T0)
  accB   = accnu = accl = rep(0,bi+M)
  V      = Vn = list()
  G1     = G2 = G2n = c(list(matrix(0,nrow=dm,ncol=dm)),Sig[-T0])
  for(t in 2:T0) G2[[t]]     = Reduce('+',Sig[max(1,t-lag):(t-1)])/min(c(t-1,lag))
  
  for(t in 1:T0){
    V[[t]]   = (B0+(b1%*%t(b1))*G1[[t]]+(b2%*%t(b2))*G2[[t]])
    llo[t]   = diwish(Sig[[t]],nu,(nu-dm-1)*V[[t]])
  }
  
  for(m in 1:(bi+M)){
    ##-----
    ## l
    ##-----
    repeat{
      pos  = rbinom(1,1,0.5)
      lagnew = lag+sample(c(1,2,3,4,5),1,prob = c(5/15,4/15,3/15,2/15,1/15))*(-1)^(1-pos)
      if(lagnew>1) break
    }
    
    for(t in 2:T0) G2n[[t]]     = Reduce('+',Sig[max(1,t-lagnew):(t-1)])/min(c(t-1,lagnew))
    for(t in 1:T0){
      V[[t]]   = (B0+(b1%*%t(b1))*G1[[t]]+(b2%*%t(b2))*G2n[[t]])
      lln[t]   = diwish(Sig[[t]],nu,(nu-dm-1)*V[[t]])
    }
    
    if((sum(lln)-sum(llo))>log(runif(1))){
      llo     = lln
      lag     = lagnew
      G2      = G2n
      accl[m] = 1
    }
    
    #if(rbinom(1,1,0.05)==1){lag = lag+rpois(1,1)*(-1)^(1-pos)}
    
    ##-----
    ## bs
    ##-----
    repeat{
      bn  = rnorm(dm*2,c(b1,b2),sd=0.007)#0.015
      b1n = bn[1:dm]
      b2n = bn[(dm+1):(2*dm)]
      B1  = b1n%*%t(b1n)
      B2  = b2n%*%t(b2n)
      B0  = (iota%*%t(iota)-B1-B2)*Sbar
      if(b1n[1]>0 && b2n[1]>0 && is.positive.definite(B0) && (sum(abs(B1+B2)<1)==dm^2)) break
    }
    for(t in 1:T0){
      Vn[[t]]  = (B0+(b1n%*%t(b1n))*G1[[t]]+(b2n%*%t(b2n))*G2[[t]])
      lln[t]   = diwish(Sig[[t]],nu,(nu-dm-1)*Vn[[t]])
    }
    if((sum(lln)-sum(llo)+
        sum(dnorm(b1n,0,sqrt(10),log=T0))-sum(dnorm(b1,0,sqrt(10),log=T0))+
        sum(dnorm(b2n,0,sqrt(10),log=T0))-sum(dnorm(b2,0,sqrt(10),log=T0)))>log(runif(1))){
      b1   = b1n
      b2   = b2n
      accB[m] = 1
      llo  = lln
      V    = Vn
    }
    B0  = (iota%*%t(iota)-(b1%*%t(b1))-(b2%*%t(b2)))*Sbar
    
    ##-----
    ## nu
    ##-----
    repeat{
      nun = rnorm(1,nu,sd=0.1)#0.5
      # 0.1 for 3-variate too large
      if(nun>(dm+1)) break
    }
    
    di  = function(Sig,S) diwish(Sig,nun,(nun-dm-1)*S)
    lln = mapply(di,Sig,V)
    
    if((sum(lln)-sum(llo)+
        dexp(nun,1/10,log=T)-dexp(nu,1/10,log=T))>log(runif(1))){
      llo   = lln
      accnu[m] = 1
      nu    = nun
    }
    
    ##-----
    ## Collect results
    ##-----
    resc[m,] = c(lag,nu,b1,b2)
    
    ##-----
    ## Prediction
    ##-----
    Vpred[[m]]    = B0+(b1%*%t(b1))*Sig[[T0]]+(b2%*%t(b2))*Reduce('+',Sig[(T0+1-lag):T0])/lag
  
    }
  res = list(Vpred[(bi+1):(bi+M)],resc[(bi+1):(bi+M),],
             accl[(bi+1):(bi+M)],
             accnu[(bi+1):(bi+M)],
             accB[(bi+1):(bi+M)])
  names(res) = c('Vpred','resc','accl','accnu','accB')
  save(res,file='temp/results_xm1_EX.Rdata')

}


