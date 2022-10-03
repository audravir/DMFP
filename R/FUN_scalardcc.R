#' @export
scalardcc = function(standrets,M){
  data = standrets[1:T0,]
  TT   = dim(data)[1]
  dm   = dim(data)[2]
  bi   = M
  Qold = array(NA,c(dm, dm, TT))
  R       = array(NA,c(dm, dm, TT))
  Qold[,,1] = cor(data)
  R[,,1] <- diag(diag(Qold[,,1])^{-1/2})%*%Qold[,,1]%*%diag(diag(Qold[,,1])^{-1/2})
  aold   <- 0.01
  bold   <- 0.95
  llold  <- rep(0,TT)
  resdcc <- matrix(NA,ncol=2,nrow=(bi+M))
  accdcc <- rep(0,bi+M)
  Vpred  = Qpred = list()
  S      = cov(data)
  
  for(t in 2:TT){
    Qold[,,t]   <- S*(1-aold-bold)+aold*(data[t-1,]%*%t(data[t-1,]))+bold*Qold[,,(t-1)]
    R[,,t]   <- diag(diag(Qold[,,t])^{-1/2})%*%Qold[,,t]%*%diag(diag(Qold[,,t])^{-1/2})
    llold[t] <- mvtnorm::dmvnorm(data[t,], rep(0,dm), R[,,t], log=T)
  }
  
  for(m in 1:(M+bi)){
    
    parnew = rnorm(2,c(aold,bold),0.005) # changed from 0.0005
    
    anew  <- parnew[1]
    bnew  <- parnew[2]
    llnew <- rep(0,TT)
    Qnew  = Qold
    
    for(t in 2:TT){
      Qnew[,,t]   <- S*(1-anew-bnew)+anew*(data[t-1,]%*%t(data[t-1,]))+bnew*Qnew[,,(t-1)]
      R[,,t]   <- diag(diag(Qnew[,,t])^{-1/2})%*%Qnew[,,t]%*%diag(diag(Qnew[,,t])^{-1/2})
      llnew[t] <- mvtnorm::dmvnorm(data[t,], rep(0,dm), R[,,t], log=T)
    }
    
    if((sum(llnew)-sum(llold)+
        dbeta(anew,3,10,log=TRUE)-dbeta(aold,3,10,log=TRUE)+dbeta(bnew,10,3,log=TRUE)-dbeta(bold,10,3,log=TRUE))>log(runif(1))
       &&
       (sum(parnew)<1))
    {
      llold  = llnew
      aold   = anew
      bold   = bnew
      Qold   = Qnew
      accdcc[m] = 1
    }
    
    resdcc[m,] <- c(aold,bold) 
    Qpred[[m]] <- S*(1-aold-bold)+aold*(data[TT,]%*%t(data[TT,]))+bold*Qold[,,TT]
    Vpred[[m]] <- diag(diag(Qpred[[m]])^{-1/2})%*%Qpred[[m]]%*%diag(diag(Qpred[[m]])^{-1/2})
    
    # if(m>bi){save(Qold,file=paste('temp/dcc_',m-bi,'.Rdata',sep=''))}
  }  
  
  res = list(Vpred[(bi+1):(bi+M)],Qpred[(bi+1):(bi+M)],resdcc[(bi+1):(bi+M),],accdcc[(bi+1):(bi+M)])
  names(res) = c('Vpred','Qpred','resdcc','accdcc')
  save(res,file=paste('temp/results_scalar_dcc_EX.Rdata',sep=''))
}

