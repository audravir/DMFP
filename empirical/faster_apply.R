rm(list=ls(all=TRUE))

library(parallel)
detectCores()

rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')

nn       = length(date)
end.date = which(zoo::as.yearmon(date)=="ene 2021")[1]-1
if(is.na(date[end.date])){end.date = which(zoo::as.yearmon(date)=="jan 2020")[1]-1}
  
date[end.date]
  
T0  = end.date  #
T0 # for estimation
K = nn-end.date
  K # for oos evaluation
  
  T0+K
  nn
  # the same
  
  data = Sigma[1:T0]
  
  # function arguments
  M = 1000
  
  # 0.001 give accp of 0.004
  propsdb = 0.0001 
  
  # 0.001 gave accp of
  propsdnu = 0.001
  
  
  dwish  = function(Sig,nu,S){dwishart(Sig, nu, S/nu, log=TRUE)}
  # this density function is the same as in GOLOSNOY et al (2012), Eq (3)
  Sig    = data
  dm     = dim(Sig[[1]])[1]
  TT     = length(Sig)
  bi     = min(M,10^4)
  resc   = matrix(NA,nrow=M,ncol=dm*2+1)
  Vpred  = vector(mode = "list", length = M)
  nu     = 20
  b1     = rep(0.95,dm)
  b2     = rep(0.3,dm)
  Sbar   = Reduce('+',Sig)/TT
  iota   = rep(1,dm)
  Oiota  = Outer(iota,iota)
  B1     = Outer(b1,b1) 
  B2     = Outer(b2,b2)
  B0     = (Oiota-B1-B2)*Sbar
  llo    = lln = ll3=rep(0,TT)
  accB   = accnu = rep(0,bi+M)
  V      = Vn = list()
  V[[1]] = Vn[[1]] = Sbar
  
  
  for(t in 2:TT){
    V[[t]]   = B0+B1*V[[t-1]]+B2*Sig[[t-1]]
  }


system.time(  for(t in 2:TT){
  llo[t]   = dwish(Sig[[t]],nu,V[[t]])
})

sum(llo)

system.time(ll2<-mapply(dwish,Sig,nu,V))

sum(ll2[-1])

dwish.t  = function(x,y){dwishart(x, nu, y/nu, log=TRUE)}

system.time(  ll4<-parallel::mcmapply(function(x,y){
  return(dwish.t(x,y))
},x=Sig,y=V,mc.cores=1)
)

sum(ll4[-1])


