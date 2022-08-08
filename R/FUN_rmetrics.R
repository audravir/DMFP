#' @export
rmetrics = function(standrets,M){
  data = standrets[1:T0,]
  
  # Sigma is a list containing realized covariances
  # The file also contains the data (standard normally distibuted)
  dm   = dim(data)[2]
  T    = dim(data)[1]
  bi   = M

  lamo    = 0.99
  llo     = rep(0,T)
  llikRMe = matrix(0,ncol=T,nrow=bi+M)
  resRMe  = rep(NA,bi+M)
  accRme  = rep(0,bi+M)
  sdprop  = 0.001 # 0.0001 too small 0.0005 too small
  
  Vpred   = list()
  Qold    = array(NA,c(dm, dm, T))
  R       = array(NA,c(dm, dm, T))
  Qold[,,1]  = cor(data)
  R[,,1]  = diag(diag(Qold[,,1])^{-1/2})%*%Qold[,,1]%*%diag(diag(Qold[,,1])^{-1/2})
  
  for(t in 2:T){
    Qold[,,t]   = (1-lamo)*(data[t-1,]%*%t(data[t-1,]))+lamo*Qold[,,(t-1)]
    # R[,,t]   = nearcor(round(cov2cor(Q[,,t]),6))$cor
    R[,,t]   = diag(diag(Qold[,,t])^{-1/2})%*%Qold[,,t]%*%diag(diag(Qold[,,t])^{-1/2})
    llo[t]   = mvtnorm::dmvnorm(data[t,],rep(0,dm),R[,,t],log=T)
  }
  
  for(m in 1:(M+bi)){
    lamn = truncnorm::rtruncnorm(1,0,1,lamo,sdprop) #asymmetric TN proposal
    Qnew = Qold
    lln  = rep(0,T)
    for(t in 2:T){
      Qnew[,,t] = (1-lamn)*(data[t-1,]%*%t(data[t-1,]))+lamn*Qnew[,,(t-1)]
      R[,,t] = diag(diag(Qnew[,,t])^{-1/2})%*%Qnew[,,t]%*%diag(diag(Qnew[,,t])^{-1/2})
      lln[t] = mvtnorm::dmvnorm(data[t,],rep(0,dm),R[,,t],log=T)
    }
    if((sum(lln)-sum(llo)+dbeta(lamn,10,3,log=T)-dbeta(lamo,10,3,log=T)
        +log(dtruncnorm(lamo,0,1,lamn,sdprop))-
        log(dtruncnorm(lamn,0,1,lamo,sdprop)))>log(runif(1))){
      llo  = lln
      lamo = lamn
      Qold = Qnew
      accRme[m]  = 1
    }
    
    resRMe[m]   = lamo
    llikRMe[m,] = llo
    Qpred      <- (1-lamo)*(data[T,]%*%t(data[T,]))+lamo*Qold[,,T]
    Vpred[[m]] <- diag(diag(Qpred)^{-1/2})%*%Qpred%*%diag(diag(Qpred)^{-1/2})
  }
  
  res = list(Vpred[(bi+1):(bi+M)],resRMe[(bi+1):(bi+M)],accRme[(bi+1):(bi+M)])
  names(res) = c('Vpred','resRMe','accRMe')
  save(res,file=paste('temp/results_RMe.Rdata',sep=''))
}

