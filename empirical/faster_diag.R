M=100000
dm=15

Qold  = rWishart(1,dm,diag(dm))[,,1]
data=rnorm(dm)
l1=l2=l3=rep(NA,M)
library(profvis)
library(Rfast)
library(mvnfast)


x <-rnorm(10) 
y <-rnorm(10)
res<-Outer(x,y)
outer(x,y)

# library("microbenchmark")
# 
# microbenchmark(mvnfast::dmvt(data,rep(0,dm), R1 , dm+10, log=TRUE), 
#                Rfast::dmvt(data,rep(0,dm), R1 , dm+10, logged=TRUE),times=1000000)  
# 

# AUDRA: compare Rfast with mvtnfast for dmvt
# Rfast is faster
# use Rfast: dmvt and Outer

fm = diag(dm)

profvis({
  
  for(t in 1:M){
    Ma = Qold 
    dv <- Ma[ col(Ma)==row(Ma) ]^{-1/2} 
    R3  <- Outer(dv,dv)*Ma
    l3[t] = dmvt(data,rep(0,dm), R3 , dm+10, log=TRUE)
  }
  
for(t in 1:M){
  R1 <- diag(diag(Qold )^{-1/2})%*%Qold %*%diag(diag(Qold )^{-1/2})
  l1[t] = dmvt(data,rep(0,dm), R1 , dm+10, log=TRUE)
}
  
for(t in 1:M){
  Ma = Qold 
  fm[col(fm)==row(fm)] <- Ma[ col(Ma)==row(Ma) ]^{-1/2} 
  R2  <- fm%*%Ma%*%fm
  l2[t] = dmvt(data, rep(0,dm), R2 ,dm+10, log=TRUE)
}
  

  
})

sum(l1)
sum(l2)
sum(l3)
