rm(list=ls(all=TRUE))
load('data/FXdata.Rdata')

dm = dim(rets)[2]

M = bi = 5000
sc.fac = 0.5
nn = dim(rets)[1]
llo = lln = rep(NA,nn)
acc=res = rep(0,M+bi)

for(t in 1:nn){
  llo[t] = mvnfast::dmvn(rets[t,],rep(0,dm),sc.fac*RCov[,,t],log=TRUE)
}

for(m in 1:(M+bi)){
  sc.fac.n = rnorm(1,sc.fac,sd=0.01)
  
  for(t in 1:nn){
    lln[t] = mvnfast::dmvn(rets[t,],rep(0,dm),sc.fac.n*RCov[,,t],log=TRUE)
  }
  
  if((sum(lln)-sum(llo)+
      sum(dnorm(sc.fac.n,1,sqrt(10),log=T))-
      sum(dnorm(sc.fac,1,sqrt(10),log=T)))>log(runif(1))){
    llo  = lln
    sc.fac = sc.fac.n
    acc[m] = 1
  }
  res[m] = sc.fac
}

mean(acc)
plot(res,type='l')

mean(tail(res,500))


