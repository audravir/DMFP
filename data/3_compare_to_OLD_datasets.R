rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')
fx15 = list(rets,RVs,stand,RCor,RCov)
names(fx15) = c('rets','RVs','stand','RCor','RCov')
rm(rets,RVs,stand,RCor,RCov)

load('data/old_5data_FX.Rdata')
fx5 = list(rets,RVs,stand,RCor,RCov)
names(fx5) = c('rets','RVs','stand','RCor','RCov')
rm(rets,RVs,stand,RCor,RCov)

load('data/old_10data.Rdata')
var10 = list(rets,RVs,stand,RCor,RCov)
names(var10) = c('rets','RVs','stand','RCor','RCov')

rm(list=setdiff(ls(), c("fx15" , "fx5","var10")))

dm1 = dim(fx15$RVs)[2]
dm2 = dim(fx5$RVs)[2]
dm3 = dim(var10$RVs)[2]

##-----------------
## compare the smoothness of RVs
##-----------------

par(mfrow=c(3,1))

plot(fx15$RVs[,1],type='l')
for(i in 1:dm1) lines(fx15$RVs[,i])

plot(fx5$RVs[,1],type='l')
for(i in 1:dm2) lines(fx5$RVs[,i])

plot(var10$RVs[,1],type='l')
for(i in 1:dm3) lines(var10$RVs[,i])

smt <- function(x)sd(diff(x)) / mean(abs(diff(x)))

apply(fx15$RVs,2,smt)
apply(fx5$RVs,2,smt)
apply(var10$RVs,2,smt)

##-----------------
## look at the correlations
##-----------------

C = combn(dm1,2)
corsd = rep(NA,ncol(C))
for (i in 1:ncol(C)) {
   corsd[i] = sd(fx15$RCor[C[1,i],C[2,i],])
}
hist(corsd)

C = combn(dm2,2)
corsd = rep(NA,ncol(C))
for (i in 1:ncol(C)) {
  corsd[i] = sd(fx5$RCor[C[1,i],C[2,i],])
}
hist(corsd)

C = combn(dm3,2)
corsd = rep(NA,ncol(C))
for (i in 1:ncol(C)) {
  corsd[i] = sd(var10$RCor[C[1,i],C[2,i],])
}
hist(corsd)


library(corrplot)
corrplot(cor(fx15$rets))
corrplot(cor(fx5$rets))
corrplot(cor(var10$rets))
