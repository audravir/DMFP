rm(list=ls(all=TRUE))
library(truncnorm)
load("empirical/temp/individual_llh.RData")

#---- will pool:
# 1. DCC-t + DCC
# 2. DCC-t + Static
# 3. CAW+AIW
# 4. CAW+Static

#---- weights:
# Geweke
# Jore1
# Equally weighted





resDN1= PMCMC_delNegro(lL_tdcc,lL_dcc,1000,c(0,2),10000,propsd = 0.5)
resDN2= PMCMC_delNegro(lL_tdcc,mat.st,1000,c(0,2),10000,propsd = 0.5)
resDN3= PMCMC_delNegro(lL_caw,lL_xm,1000,c(0,2),10000,propsd = 0.5)
resDN4= PMCMC_delNegro(lL_caw,mat.st,1000,c(0,2),10000,propsd = 0.5)

ws_DN1 = t(pnorm(resDN1$weights_xs))
if (dim(ws_DN1)[1] > post.sample){
  dnind = round(seq(1,dim(ws_DN1)[1],length=post.sample))
  ws_DN1 = ws_DN1[dnind,]
}

ws_DN2 = t(pnorm(resDN2$weights_xs))
if (dim(ws_DN2)[1] > post.sample){
  dnind = round(seq(1,dim(ws_DN2)[1],length=post.sample))
  ws_DN2 = ws_DN2[dnind,]
}

ws_DN3 = t(pnorm(resDN3$weights_xs))
if (dim(ws_DN3)[1] > post.sample){
  dnind = round(seq(1,dim(ws_DN3)[1],length=post.sample))
  ws_DN3 = ws_DN3[dnind,]
}

ws_DN4 = t(pnorm(resDN4$weights_xs))
if (dim(ws_DN4)[1] > post.sample){
  dnind = round(seq(1,dim(ws_DN4)[1],length=post.sample))
  ws_DN4 = ws_DN4[dnind,]
}

mypool = function(ll1,ll2,ws_DN){
  ws_gew    = lL_gew = matrix(NA,ncol=K,nrow=post.sample)
  ws_jore1  = lL_jore1 = matrix(NA,ncol=K,nrow=post.sample)
  lL_equal  = lL_DN =  matrix(NA,ncol=K,nrow=post.sample)
  
  lhf = exp(ll1)
  llf = exp(ll2)
  
  for(m in 1:post.sample){
    for (t in 1:K){
      # geweke
      a1p = lhf[m,1:t]
      a2p = llf[m,1:t]
      opw = function(x){
        -sum(log(x*a1p+(1-x)*a2p))
      }
      ws_gew[m,t] = optim(0.5, opw, gr = NULL,
                          method = c("L-BFGS-B"),
                          lower = 0, upper = 1, hessian = FALSE)$par
      lL_gew[m,t] = log(ws_gew[m,t]*lhf[m,t]+(1-ws_gew[m,t])*llf[m,t])
      
      # jore1
      ws_jore1[m,t] = a1p[t]/(a1p[t]+a2p[t])
      lL_jore1[m,t] = log(ws_jore1[m,t]*lhf[m,t]+(1-ws_jore1[m,t])*llf[m,t])
      
      # equally-weighted
      lL_equal[m,t] = log(0.5*lhf[m,t]+0.5*llf[m,t])
      
      # DN
      lL_DN[m,t] = log(ws_DN[m,t]*lhf[m,t]+(1-ws_DN[m,t])*llf[m,t])
      
    }
  }
  res=list(ws_gew,lL_gew,ws_jore1,lL_jore1,lL_equal,lL_DN)
  names(res) = c('ws_gew','lL_gew','ws_jore1','lL_jore1','lL_equal','lL_DN')
  return(res)
}

# best.model = mypool(ll1=lL_caw,ll2=lL_tdcc)
mat.st <- matrix(rep(lL_static,post.sample), nrow = post.sample, byrow = TRUE)
mat.rm <- matrix(rep(lL_rmf,post.sample), nrow = post.sample, byrow = TRUE)


# 1. DCC-t + DCC
res1 = mypool(ll1=lL_tdcc,ll2=lL_dcc,ws_DN1)
# 2. DCC-t + Static
res2 = mypool(ll1=lL_tdcc,ll2= mat.st,ws_DN2)
# 3. CAW+AIW
res3 = mypool(ll1=lL_caw,ll2=lL_xm,ws_DN3)
# 4. CAW+Static
res4 = mypool(ll1=lL_caw,ll2=mat.st,ws_DN4)

# # 5. DCC-t + RMf
# res5 = mypool(ll1=lL_tdcc,ll2= mat.rm)
# # 6. CAW+RMf
# res6 = mypool(ll1=lL_caw,ll2=mat.rm)



# save.image('empirical/temp/other_pools.RData')


move.axis = 100
mkvol     = apply(tail(apply(RVs^2,2,scale),K),1,mean)
atx       <- seq(date[(nn+1)], date[(nn+K)], by=25)

pdf('tables_and_figures/other_pools_FX.pdf',height=6,width=12)
par(mfrow=c(2,2), mar=c(4, 3, 1, 1) + 0.1)
plot(tail(date,K),mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]),main='DCC-t + DCC')
par(new = TRUE)
plot(tail(date,K),rep(0,K),type='l',ylim=c(-600,10),ylab='',xlab='',xaxt="n",
     xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]))
lines(tail(date,K),cumsum(apply(res1$lL_gew,2,median))-cumsum(apply(lL_caw,2,median)),
      col=1,lwd=2)
lines(tail(date,K),cumsum(apply(res1$lL_jore1,2,median))-cumsum(apply(lL_caw,2,median)),
      col='pink',lwd=2)
lines(tail(date,K),cumsum(apply(res1$lL_DN,2,median))-cumsum(apply(lL_caw,2,median)),
      col='coral',lwd=2,lty=2)
axis(1, at=atx, labels=format(atx, "%Y/%m"),las=2,cex.axis=0.75)
legend(x=date[(nn+1)]-move.axis-35,y=0,col=c(1,'pink','coral'),
       lty=c(1,1,2),lwd=c(2,2,2),
       legend=c('Geweke','Jore1','DelNegro'))

plot(tail(date,K),mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]),main='DCC-t + Static')
par(new = TRUE)
plot(tail(date,K),rep(0,K),type='l',ylim=c(-750,10),ylab='',xlab='',xaxt="n",
     xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]))
lines(tail(date,K),cumsum(apply(res2$lL_gew,2,median))-cumsum(apply(lL_caw,2,median)),
      col=1,lwd=2)
lines(tail(date,K),cumsum(apply(res2$lL_jore1,2,median))-cumsum(apply(lL_caw,2,median)),
      col='pink',lwd=2)
lines(tail(date,K),cumsum(apply(res2$lL_DN,2,median))-cumsum(apply(lL_caw,2,median)),
      col='coral',lwd=2,lty=2)
axis(1, at=atx, labels=format(atx, "%Y/%m"),las=2,cex.axis=0.75)
legend(x=date[(nn+1)]-move.axis-35,y=0,col=c(1,'pink','coral'),lty=c(1,1,2),lwd=c(2,2,2),
       legend=c('Geweke','Jore1','DelNegro'))

plot(tail(date,K),mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]),main='CAW + AIW')
par(new = TRUE)
plot(tail(date,K),rep(0,K),type='l',ylim=c(-10,100),ylab='',xlab='',xaxt="n",
     xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]))
lines(tail(date,K),cumsum(apply(res3$lL_gew,2,median))-cumsum(apply(lL_caw,2,median)),
      col=1,lwd=2)
lines(tail(date,K),cumsum(apply(res3$lL_jore1,2,median))-cumsum(apply(lL_caw,2,median)),
      col='pink',lwd=2)
lines(tail(date,K),cumsum(apply(res3$lL_DN,2,median))-cumsum(apply(lL_caw,2,median)),
      col='coral',lwd=2,lty=2)
axis(1, at=atx, labels=format(atx, "%Y/%m"),las=2,cex.axis=0.75)
legend(x=date[(nn+1)]-move.axis-35,y=100,col=c(1,'pink','coral'),lty=c(1,1,2),lwd=c(2,2,2),
       legend=c('Geweke','Jore1','DelNegro'))

plot(tail(date,K),mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]),main='CAW + Static')
par(new = TRUE)
plot(tail(date,K),rep(0,K),type='l',ylim=c(-10,100),ylab='',xlab='',xaxt="n",
     xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]))
lines(tail(date,K),cumsum(apply(res4$lL_gew,2,median))-cumsum(apply(lL_caw,2,median)),
      col=1,lwd=2)
lines(tail(date,K),cumsum(apply(res4$lL_jore1,2,median))-cumsum(apply(lL_caw,2,median)),
      col='pink',lwd=2)
lines(tail(date,K),cumsum(apply(res4$lL_DN,2,median))-cumsum(apply(lL_caw,2,median)),
      col='coral',lwd=2,lty=2)
axis(1, at=atx, labels=format(atx, "%Y/%m"),las=2,cex.axis=0.75)
legend(x=date[(nn+1)]-move.axis-35,y=100,col=c(1,'pink','coral'),
       lty=c(1,1,2),lwd=c(2,2,2),
       legend=c('Geweke','Jore1','DelNegro'))
dev.off()




