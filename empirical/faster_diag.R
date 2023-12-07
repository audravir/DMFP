rm(list=ls(all=TRUE))
library(rbenchmark)
library(Rfast)
library(mvnfast)
library(matrixcalc)
library(RSpectra)


dm=15
data=rnorm(dm)
means=rep(0,dm)
df=dm+10
Sigma=rWishart(1,dm+10,diag(dm))[,,1]


benchmark(Rfast::dmvt(data,means,Sigma,df,logged=TRUE), 
          mvnfast::dmvt(data,means,Sigma,df,log=TRUE), 
          order="relative", replications=50000)


b=rnorm(dm)
benchmark(Rfast::Outer(b,b), 
          base::outer(b,b), 
          order="relative", replications=50000)


benchmark(is.positive.definite(Sigma), 
          eigen.sym(Sigma,k = dm-1,vectors = FALSE), 
          order="relative", replications=5000)


benchmark(eigs_sym(Sigma,k = dm-1, opts = list(retvec = FALSE)), 
          eigen.sym(Sigma,k = dm-1,vectors = FALSE), 
          order="relative", replications=5000)

Sigma




# # faster Q matrices
# for(t in 1:M){
#   Ma = Qold 
#   dv <- Ma[ col(Ma)==row(Ma) ]^{-1/2} 
#   R3  <- Outer(dv,dv)*Ma
#   l3[t] = dmvt(data,rep(0,dm), R3 , dm+10, log=TRUE)
# }
