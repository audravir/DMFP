#' @export
vectordcc_tcop = function(data,M,propsd,propsdnu){
  t0   = Sys.time()
  TT   = dim(data)[1]
  dm   = dim(data)[2]
  bi   = min(M,10^4)
  
  udata = pnorm(data)*TT/(TT+1) 
  
  Qold = array(NA,c(dm, dm, TT))
  R    = array(NA,c(dm, dm, TT))
  aold   <- rep(0.1,dm)
  bold   <- rep(0.99,dm)
  nuold  <- 20
  tdata  <- qt(udata,nuold)
  llold  <- rep(0,TT)
 
  Qold[,,1] = cor(tdata)
  R[,,1] <- diag(diag(Qold[,,1])^{-1/2})%*%Qold[,,1]%*%diag(diag(Qold[,,1])^{-1/2})
  Sbar = cov(tdata)
  
  resdcc <- matrix(NA,ncol=dm*2+1,nrow=(bi+M))
  iota   = rep(1,dm)
  B0     = (iota%*%t(iota)-aold%*%t(aold)-bold%*%t(bold))*Sbar
  
  accdcc <- rep(0,bi+M)
  accnu  <- rep(0,bi+M)
  Vpred  = Qpred = list()
  
  for(t in 2:TT){
    Qold[,,t]   <- B0+(aold%*%t(aold))*(data[t-1,]%*%t(data[t-1,]))+(bold%*%t(bold))*Qold[,,(t-1)]
    R[,,t]   <- diag(diag(Qold[,,t])^{-1/2})%*%Qold[,,t]%*%diag(diag(Qold[,,t])^{-1/2})
    llold[t] <- mvtnorm::dmvt(tdata[t,], rep(0,dm), R[,,t], df=nuold, log=TRUE)-
      sum(dt(tdata[t,],df=nuold,log=TRUE))
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
      B0  = (iota%*%t(iota)-anew%*%t(anew)-bnew%*%t(bnew))*Sbar
      if(anew[1]>0 && bnew[1]>0 && is.positive.definite(B0) && (sum(abs(anew%*%t(anew)+bnew%*%t(bnew))<1)==dm^2)) break
    }
    
    llnew <- rep(0,TT)
    Qnew  = Qold
    
    for(t in 2:TT){
      Qnew[,,t]   <- B0+(anew%*%t(anew))*(tdata[t-1,]%*%t(tdata[t-1,]))+(bnew%*%t(bnew))*Qnew[,,(t-1)]
      R[,,t]   <- diag(diag(Qnew[,,t])^{-1/2})%*%Qnew[,,t]%*%diag(diag(Qnew[,,t])^{-1/2})
      llnew[t] <- mvtnorm::dmvt(tdata[t,], rep(0,dm), R[,,t], df=nuold, log=TRUE)-
        sum(dt(tdata[t,],df=nuold,log=TRUE))
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
    
    
    ##-----
    ## nu
    ##-----
    repeat{
      nunew = rnorm(1,nuold,propsdnu)
      if(nunew>dm) break
    }
    
    tdata  = qt(udata,nunew)
    Sbar   = cov(tdata)
    B0     = (iota%*%t(iota)-aold%*%t(aold)-bold%*%t(bold))*Sbar
    
    llnew <- rep(0,TT)
    Qnew  = Qold
    
    for(t in 2:TT){
      Qnew[,,t]   <- B0+(aold%*%t(aold))*(tdata[t-1,]%*%t(tdata[t-1,]))+(bold%*%t(bold))*Qnew[,,(t-1)]
      R[,,t]   <- diag(diag(Qnew[,,t])^{-1/2})%*%Qnew[,,t]%*%diag(diag(Qnew[,,t])^{-1/2})
      llnew[t] <- mvtnorm::dmvt(tdata[t,], rep(0,dm), R[,,t], df=nunew, log=TRUE)-
        sum(dt(tdata[t,],df=nunew,log=TRUE))
    }
    
    if((sum(llnew)-sum(llold)+
        dexp(nunew,0.1,log=TRUE)-dexp(nuold,0.1,log=TRUE))>log(runif(1)))
    {
      llold  = llnew
      nuold  = nunew
      Qold   = Qnew
      accdnu[m] = 1
    }
    
    tdata  = qt(udata,nuold)
    Sbar   = cov(tdata)
    B0     = (iota%*%t(iota)-aold%*%t(aold)-bold%*%t(bold))*Sbar
    
    
    resdcc[m,] <- c(aold,bold,nuold) 
    Qpred[[m]] <- B0+aold*(tdata[TT,]%*%t(tdata[TT,]))+bold*Qold[,,TT]
    Vpred[[m]] <- diag(diag(Qpred[[m]])^{-1/2})%*%Qpred[[m]]%*%diag(diag(Qpred[[m]])^{-1/2})
    
    if(!m%%100){
      print(paste(round(m/(M+bi)*100),"%",sep=""))
      print(Sys.time()-t1)
      print(Sys.time()-t0)
    }
  }  
  res = list(Vpred[(bi+1):(bi+M)],Qpred[(bi+1):(bi+M)],resdcc[(bi+1):(bi+M),],
             accdcc[(bi+1):(bi+M)],accnu[(bi+1):(bi+M)])
  names(res) = c('Vpred','Qpred','resdcc','accdcc','accnu')
  save(res,file=paste('empirical/temp/results_vectordcc_tcop.Rdata',sep=''))
}

