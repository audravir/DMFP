#' @export
PMCMC_delNegro = function(first_lik,second_lik,Npart,prior){
  y1 = apply(first_lik,2,median)
  y2 = apply(second_lik,2,median)
  K = length(y1)
  priorm = prior[1]
  priorsd = prior[2]
  
  M = bi = dim(first_lik)[1]
  
  ##-------
  ## Particle MCMC
  ##-------
  
  betas = rep(NA,M+bi)
  xssp  = wssp = matrix(NA,ncol=M+bi,nrow=K)
  acc   = rep(0,M+bi)
  lams  = pnorm(rnorm(K,0,1))
  llo   = sum(log(lams*y1+(1-lams)*y2))
  beta  = truncnorm::rtruncnorm(1,a=-1,b=1,0.1,sqrt(0.1))
  tau2  = 1-beta^2
  xso = xsn = wsso = rep(0,K)
  
  for(m in 1:(M+bi)){
    betan  = truncnorm::rtruncnorm(1,a=-1,b=1,beta,0.2)
    tau2n  = 1-betan^2
    
    xn    = rnorm(Npart,0,1)
    wssn   = rep(0,K)
    
    for(t in 1:K){
      # blind propagate
      xs  = rnorm(Npart,betan*xn,sqrt(tau2n))
      lams = pnorm(xs)
      # weights
      ws  = lams*y1[t]+(1-lams)*y2[t]
      wssn[t] = mean(ws)
      #if(sum(ws)==0){ws=1/Npart}
      ind = sample(1:Npart,Npart,replace = T,prob=ws/sum(ws))
      xn = xs[ind]
      xsn[t] = mean(xn)
    }
    
    lln = sum(log(wssn))
    if((lln-llo+log(truncnorm::dtruncnorm(betan,0,1,priorm,priorsd))-
        log(truncnorm::dtruncnorm(beta,0,1,priorm,priorsd)))>log(runif(1))){
      llo  = lln
      beta = betan
      tau2 = tau2n
      acc[m]  = 1
      xso = xsn
      wsso = wssn
    }
    
    betas[m]  = beta
    xssp[,m]  = xso
    wssp[,m]  = wsso
  }
  resdn = list(betas[(bi+1):(bi+M)],xssp[,(bi+1):(bi+M)],
               wssp[,(bi+1):(bi+M)],acc[(bi+1):(bi+M)])
  names(resdn) = c('beta','weights_xs','likelihood','acc')
  return(resdn)
}