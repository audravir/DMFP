rm(list=ls(all=TRUE))
library(copula)
library(zoo)
library(QTLRel)

load('10data.Rdata')

dm = dim(rets)[2]

C     = combn(dm,2)
names = c('BAC','JPM','IBM')

cop_model <- tCopula(dim = 2)


pdf(file="../JAE_v1/Figures/motiv.pdf",width=9,height=3)

TT=dim(rets)[1]
par(mfrow=c(1,3))
i=1
plot(date,RCor[C[1,i],C[2,i],],type='l',ylim=c(-0.1,1),col='gray80',
     main=paste(names[C[1,i]],'/',names[C[2,i]]),ylab='Rcor',xlab='')

stcor  =c(rep(cor(rets[1:(TT/2),1:2])[2,1],TT/2),
          rep(cor(rets[(TT/2+1):TT,1:2])[2,1],TT/2))
lines(date,stcor,lwd=2)

rcor=rollapply(rets[,1:2], 50, function(x) cor(x[,1],x[,2]), by.column=FALSE)
length(rcor)
lines(date,c(rep(NA,50-1),rcor))
abline(v=date[TT/2],lwd=2,lty=4,col='gray60')
legend(x=date[1],y=0.3,col=c('gray80',1,1),lwd=c(2,2,2),
       lty=c(1,1,1),
       legend=c('Rcor','Static','Rolling'))

dt = cbind(udata[1:(TT/2),C[1,i]],udata[1:(TT/2),C[2,i]])
fit <- fitCopula(cop_model, dt, method = 'ml')
stu <- tCopula(dim = 2, param =coef(fit)[[1]],
               df = coef(fit)[[2]])
contour(stu, dCopula, xlim = c(0, 1), ylim=c(0, 1),lwd=1,col='grey80',
        main=paste('rho=',round(coef(fit)[[1]],2),
                   ',nu=',round(coef(fit)[[2]]),sep='',
                   ylab=NULL,xlab=NULL))
points(dt,pch=20)

dt = cbind(udata[(TT/2+1):TT,C[1,i]],udata[(TT/2+1):TT,C[2,i]])
fit <- fitCopula(cop_model, dt, method = 'ml')
stu <- tCopula(dim = 2, param =coef(fit)[[1]],
               df = coef(fit)[[2]])
contour(stu, dCopula, xlim = c(0, 1), ylim=c(0, 1),lwd=1,col='grey80',
        main=paste('rho=',round(coef(fit)[[1]],2),
                   ',nu=',round(coef(fit)[[2]]),sep='',
                   ylab=NULL,xlab=NULL))
points(dt,pch=20)
dev.off()

date[c(1,(TT/2))]
date[c((TT/2+1),TT)]


### desc1 plot

pdf(file="../JAE_v1/Figures/desc1.pdf",width=9,height=5)
par(mfrow=c(2,3))
plot(date,rets[,1],type='l',main='BAC: rt and RVt',col='gray80',ylab='',
     xlab='')
lines(date,RVs[,1])
qqPlot(stand[,1],x="norm",main='BAC: QQ-plot',
       ylab='',xlab='',pch=20,cex=1.2,
       plot.it=TRUE,confidence=.95)
hist(udata[,1],freq=FALSE,ylab='',xlab='',xlim=c(0,1),
     main='BAC: ut')
abline(h=1,lwd=2)

plot(date,rets[,2],type='l',main='JPM: rt and RVt',col='gray80',ylab='',
     xlab='')
lines(date,RVs[,2])
qqPlot(stand[,2],x="norm",main='JPM: QQ-plot',
       ylab='',xlab='',pch=20,cex=1.2,
       plot.it=TRUE,confidence=.95)
hist(udata[,2],freq=FALSE,ylab='',xlab='',xlim=c(0,1),
     main='JPM: ut')
abline(h=1,lwd=2)
dev.off()




