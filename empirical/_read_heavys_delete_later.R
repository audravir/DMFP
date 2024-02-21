load('empirical/temp/results_heavy_scalar_joint.Rdata')
dim(res$r)
par(mfrow=c(2,2))
for(i in 1:3) plot(res$r[,i],type='l')
plot(res$LLH,type='l')
apply(res$r,2,median)
mean(res$LLH)
mean(res$acc)

load('empirical/temp/results_heavy_scalar_separate.Rdata')
dim(res$r)
par(mfrow=c(2,2))
for(i in 1:3) plot(res$r[,i],type='l')
plot(res$LLH,type='l')
apply(res$r,2,median)
mean(res$LLH)
mean(res$accnu)
mean(res$accdcc)

load('empirical/temp/results_heavy.Rdata')
dim(res$r)
mean(res$accnu)
mean(res$accdcc1)
mean(res$accdcc2)

dm=14
nu = res$r[,1]
b1 = res$r[,2:(dm+1)]
b2 = res$r[,(dm+2):(2*dm+1)]

par(mfrow=c(2,1))
plot(res$LLH,type='l')
plot(nu,type='l')


par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b1[,i],type='l')}

par(mfrow=c(3,5)) 
for(i in 1:dm) {plot(b2[,i],type='l')}
