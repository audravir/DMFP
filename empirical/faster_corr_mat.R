rm(list=ls(all=TRUE))

load('empirical/temp/cm_xm.Rdata')


dm = 5

CM = CM[1:dm,1:dm]
CMchol = chol(CM)

M   = 5000
x   = matrix(rnorm(M*dm,sd=0.01),ncol=dm)
x.c = matrix(NA,ncol=dm,nrow=M)

for(i in 1:M){
  x.c[i,] = x[i,] %*% CMchol
}

cor(x.c)
CM

apply(x.c,2,sd)


library(Rfast)
x[i,] %*% CMchol

mm=matrix(x[i,],ncol=dm)

mat.mult(mm,CMchol)

library(microbenchmark)

microbenchmark(mat.mult(mm,CMchol),mm%*%CMchol,unit = "relative",times = M)


