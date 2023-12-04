#' @export
caw = function(data,M,propsdb,propsdnu){
  t0   = Sys.time()
  dwish  = function(Sig,nu,S){dwishart(Sig, nu, S/nu, log=TRUE)}
  # this density function is the same as in GOLOSNOY et al (2012), Eq (3)
  Sig    = data
  dm     = dim(Sig[[1]])[1]
  TT     = length(Sig)
  bi     = min(M,10^4)
  resc   = matrix(NA,nrow=bi+M,ncol=dm*2+1)
  Vpred  = list()
  nu     = 20
  b1     = rep(0.95,dm)
  b2     = rep(0.25,dm)
  Sbar   = Reduce('+',Sig)/TT
  iota   = rep(1,dm)
  B0     = (iota%*%t(iota)-b1%*%t(b1)-b2%*%t(b2))*Sbar
  llo    = lln = rep(0,TT)
  accB   = accnu = rep(0,bi+M)
  V      = Vn = list()
  V[[1]] = Vn[[1]] = Sbar
  
  for(t in 2:TT){
    V[[t]]   = B0+(b1%*%t(b1))*V[[t-1]]+(b2%*%t(b2))*Sig[[t-1]]
    llo[t]   = dwish(Sig[[t]],nu,V[[t]])
  }
  
  for(m in 1:(bi+M)){
    t1=Sys.time()
    ##-----
    ## bs
    ##-----
    repeat{
      bn  = rnorm(dm*2,c(b1,b2),sd=propsdb)
      b1n = bn[1:dm]
      b2n = bn[(dm+1):(2*dm)]
      B1  = b1n%*%t(b1n)
      B2  = b2n%*%t(b2n)
      B0  = (iota%*%t(iota)-B1-B2)*Sbar
      if(b1n[1]>0 && b2n[1]>0 && is.positive.definite(B0) && (sum(abs(B1+B2)<1)==dm^2)) break
    }
    for(t in 2:TT){
      Vn[[t]]   = B0+(b1n%*%t(b1n))*Vn[[t-1]]+(b2n%*%t(b2n))*Sig[[t-1]]
      lln[t]   = dwish(Sig[[t]],nu,Vn[[t]])
    }
    
    
    if((sum(lln)-sum(llo)+
        sum(dnorm(b1n,0,sqrt(10),log=TRUE))-sum(dnorm(b1,0,sqrt(10),log=TRUE))+
        sum(dnorm(b2n,0,sqrt(10),log=TRUE))-sum(dnorm(b2,0,sqrt(10),log=TRUE)))>log(runif(1))){
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
      nun = rnorm(1,nu,sd=propsdnu)
      if(nun>(dm)) break
    }
    
    di  = function(Sig,S) dwish(Sig,nun,S)
    lln = mapply(di,Sig,V)
    
    if((sum(lln)-sum(llo)+
        dexp(nun,1/10,log=TRUE)-dexp(nu,1/10,log=TRUE))>log(runif(1))){
      llo   = lln
      accnu[m] = 1
      nu    = nun
    }
    
    ##-----
    ## Collect results
    ##-----
    resc[m,] = c(nu,b1,b2)
    
    ##-----
    ## Prediction
    ##-----
    Vpred[[m]]    = B0+(b1%*%t(b1))*V[[TT]]+(b2%*%t(b2))*Sig[[TT]]
    
    if(!m%%100){
      print(paste(round(m/(M+bi)*100,2),"%",sep=""))
      print(Sys.time()-t1)
      print(Sys.time()-t0)
    }
  }
  res = list(Vpred[(bi+1):(bi+M)],resc[(bi+1):(bi+M),],
               accnu[(bi+1):(bi+M)],
               accB[(bi+1):(bi+M)])
    names(res) = c('Vpred','resc','accnu','accB')
    save(res,file='empirical/temp/results_caw.Rdata')
}