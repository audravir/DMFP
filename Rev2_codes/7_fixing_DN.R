  
  first_lik = lL_xm1
  second_lik = lL_tdcc
  Npart = 1000
  prior = c(0,1)
  
    y1 = apply(first_lik,2,median)
    y2 = apply(second_lik,2,median)
    
    y1 = apply(lL_xm1,2,median)
    y2 = lL_static[1,]
    
    plot(y1,type='l')
    lines(y2,col=2)
    
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
    beta  = truncnorm::rtruncnorm(1,a=0,b=1,0.9,sqrt(0.1))
    tau2  = 1-beta^2
    xso = xsn = wsso = rep(0,K)
    
    for(m in 1:(M+bi)){
      betan  = truncnorm::rtruncnorm(1,a=-1,b=1,beta,1)
      tau2n  = 1-betan^2
      
      xn     = rnorm(Npart,0,1)
      wssn   = rep(0,K)
      
      for(t in 1:K){
        # blind propagate
        xs  = rnorm(Npart,betan*xn,sqrt(tau2n))
        lams = pnorm(xs)
        # weights
        ws  = lams*y1[t]+(1-lams)*y2[t]
        # re-sample
        ind = sample(1:Npart,Npart,replace = T,prob=ws/sum(ws))
        wssn[t] = mean(ws[ind])
        xn = xs[ind]
        xsn[t] = mean(xn)
      }
      
      lln = sum(log(wssn))
      if((lln-llo+log(truncnorm::dtruncnorm(betan,-1,1,priorm,priorsd))-
          log(truncnorm::dtruncnorm(beta,-1,1,priorm,priorsd)))>log(runif(1))){
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
  
    mean(resdn$acc)
    par(mfrow=c(2,2))
    plot(resdn$beta,type='l')
    plot(apply(pnorm(resdn$weights_xs),1,median),type='l',lwd=2,ylim=c(0,1))  
    lines(apply(ws_jore1,2,median),col=2)
hist(resdn$beta,freq=FALSE)
grid=seq(-1,1,length=1000)    
lines(grid,dtruncnorm(grid,-1,1,priorm,priorsd),lwd=2)        
    
    