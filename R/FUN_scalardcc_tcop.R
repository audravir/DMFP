#' @export
scalartdcc = function(standrets,M){
  
  data  = standrets[1:T0,]
  udata = pnorm(data)*T0/(T0+1) 
  
  TT   = dim(data[1:T0,])[1]
  dm   = dim(data[1:T0,])[2]
  bi   = M
 
  Qold = array(NA,c(dm, dm, TT))
  R    = array(NA,c(dm, dm, TT))
  aold   <- 0.01
  bold   <- 0.98
  nuold  <- 10
  tdata  <- qt(udata,nuold)
  
  Qold[,,1] = cor(tdata)
  R[,,1] <- diag(diag(Qold[,,1])^{-1/2})%*%Qold[,,1]%*%diag(diag(Qold[,,1])^{-1/2})
  
  llold   <- rep(0,TT)
  restdcc <- matrix(NA,ncol=3,nrow=(bi+M))
  acctdcc <- rep(0,bi+M)
  S       = cov(tdata)
  
  for(t in 2:TT){
    Qold[,,t] <- S*(1-aold-bold)+aold*(tdata[t-1,]%*%t(tdata[t-1,]))+bold*Qold[,,(t-1)]
    R[,,t]   <- diag(diag(Qold[,,t])^{-1/2})%*%Qold[,,t]%*%diag(diag(Qold[,,t])^{-1/2})
    llold[t] <- mvtnorm::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nuold, log=T)-
      sum(dt(tdata[t,],df=nuold,log=TRUE))
  }
  
  for(m in 1:(M+bi)){
    parnew = rnorm(2,c(aold,bold),0.0075) #0.0005 too small for 10-variate
    # 0.00075 too small for 3 variate
    nunew  = truncnorm::rtruncnorm(1,a = 4,mean = nuold,sd = 2) # 2 for 3variate
    anew   = parnew[1]
    bnew   = parnew[2]
    llnew  = rep(0,TT)
    Qnew   = Qold
    tdata  = qt(udata,nunew)
    S      = cov(tdata)
    
    for(t in 2:TT){
      Qnew[,,t]<- S*(1-anew-bnew)+anew*(tdata[t-1,]%*%t(tdata[t-1,]))+bnew*Qnew[,,(t-1)]
      R[,,t]   <- diag(diag(Qnew[,,t])^{-1/2})%*%Qnew[,,t]%*%diag(diag(Qnew[,,t])^{-1/2})
      llnew[t] <- mvtnorm::dmvt(tdata[t,], rep(0,dm), R[,,t], df = nunew, log=TRUE)-
        sum(dt(tdata[t,],df=nunew,log=TRUE))
    }
    
    if((sum(llnew)-sum(llold)+
        dbeta(anew,3,10,log=TRUE)-dbeta(aold,3,10,log=TRUE)+
        dbeta(bnew,10,3,log=TRUE)-dbeta(bold,10,3,log=TRUE)+
        dexp(nunew,0.1,log=TRUE)-dexp(nuold,0.1,log=TRUE)+
        log(truncnorm::dtruncnorm(nunew,a = 4,mean = nuold,sd = 2))-
        log(truncnorm::dtruncnorm(nuold,a = 4,mean = nunew,sd = 2)))>log(runif(1))
       &&
       (sum(parnew)<1))
    {
      llold  = llnew
      aold   = anew
      bold   = bnew
      nuold  = nunew
      Qold   = Qnew
      acctdcc[m] = 1
    }
    restdcc[m,] <- c(aold,bold,nuold) 
}  
  res = list(restdcc[(bi+1):(bi+M),],acctdcc[(bi+1):(bi+M)])
  names(res) = c('restdcc','acctdcc')
  save(res,file='temp/results_scalar_dcc_t_EX.Rdata')
}