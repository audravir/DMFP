rm(list=ls(all=TRUE))
library(xtable)

load('empirical/temp/FX_portfolio_at_x_10k.Rdata')

par(mfrow=c(3,1))

#second asset for example
plot(ws_gmvr[1,,2,1],type='l',ylim=c(0,0.3),main='GMVr')
for(i in 1:postM) lines(ws_gmvr[1,,2,i])

#second asset for example
plot(ws_CVAR10r[1,,2,1],type='l',ylim=c(0,0.3),main='CVAR10r')
for(i in 1:postM) lines(ws_CVAR10r[1,,2,i])

#second asset for example
plot(ws_CVAR05r[1,,2,1],type='l',ylim=c(0,0.3),main='CVAR05r')
for(i in 1:postM) lines(ws_CVAR05r[1,,2,i])

############################################################

esfun=function(x,p){
  es=mean(x[which(x<quantile(x,p))])
  return(es)
}

p05 = function(x) quantile(x,0.05)
p95 = function(x) quantile(x,0.95)

p10 = function(x) quantile(x,0.10)
p90 = function(x) quantile(x,0.90)


load('data/rf.RData')

##

all.ws = list(ws_gmvr,ws_CVAR05r,ws_CVAR10r)

names(all.ws) = c('gmvr','CVAR05r','CVAR10r')

all.prets = all.cos = all.tos = list()
pret      = cos = tos = array(NA,dim=c(length(models),K,postM))

for(j in 1:length(all.ws)){
  for(m in 1:postM){
    for(i in 1:length(models)){
      for(t in 1:K){
        pret[i,t,m] = sum(all.ws[[j]][i,t,,m]*rets[nn+t,])
        cos[i,t,m]  = (sum((all.ws[[j]][i,t,,m])^2))^(1/2)
        if(t<K)  tos[i,t,m]   = sum(abs(all.ws[[j]][i,t+1,,m]-all.ws[[j]][i,t,,m]*((1+rets[nn+t,])/(1+pret[i,t,m]))))
      }
    }
  }
  all.prets[[j]] = pret
  all.cos[[j]]   = cos
  all.tos[[j]] = tos
}

names(all.prets) = names(all.cos) = names(all.tos)=c('gmvr','CVAR05r','CVAR10r')



# avrg_return sd SR CVAR05 CVAR10 CO TO

all.avrg = all.sds = all.SR = all.cvar05 = all.cvar10 = list()
avrg = sds     = SR = cvar05 = cvar10 = array(NA,dim=c(length(models),postM))


for(j in 1:length(all.ws)){
  for(m in 1:postM){
    for(i in 1:length(models)){
        avrg[i,m] = mean(all.prets[[j]][i,,m])*252
        sds[i,m] = sd(all.prets[[j]][i,,m])*sqrt(252)
        SR[i,m]  = mean(all.prets[[j]][i,,m]*252-rf)/(sd(all.prets[[j]][i,,m])*sqrt(252))
        cvar05[i,m] = esfun(all.prets[[j]][i,,m]*252,0.05)
        cvar10[i,m] = esfun(all.prets[[j]][i,,m]*252,0.10)
      }
  }
  all.avrg[[j]] = avrg
  all.sds[[j]] = sds
  all.SR[[j]]  = SR
  all.cvar05[[j]] = cvar05
  all.cvar10[[j]] = cvar10
}

# # 1. jore's1
# # 2. geweke's
# # 3. equally w.
# # 4. caw
# # 5. dcct
# # 6. dcc-heavy-t

all.colors  = c('pink',1,'coral','royalblue','coral3','violet')
all.ltys    = c(1,1,2,1,1,1)
all.legends = c('Jore1','Geweke','Equal','CAW','DCC-t','DCC-HEAVY-t')


##-------------
# GVM 
##-------------

pdf('tables_and_figures/portfolio_cis_GMV.pdf',height=8,width=12)
par(mfrow=c(4,2),mar= c(2,2,2,2))
# plot 1: avrg return
plot(density(all.avrg[[1]][5,]),col=0,ylim=c(0,30),xlim=c(6.4,7.2),
     main='Average return',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.avrg[[1]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.avrg[[1]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=6.4,y=30,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 2: stdev
plot(density(all.sds[[1]][5,]),col=0,ylim=c(0,350),xlim=c(1.31,1.39),
     main='Standard deviation',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.sds[[1]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.sds[[1]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=1.31,y=350,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 3: SR
plot(density(all.SR[[1]][5,]),col=0,ylim=c(0,40),xlim=c(3.35,3.95),
     main='Adjusted Sharpe ratio',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.SR[[1]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.SR[[1]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=3.35,y=40,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 4: CVaR05
plot(density(all.cvar05[[1]][5,]),col=0,ylim=c(0,5),xlim=c(-44,-39),
     main='CVaR05',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.cvar05[[1]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.cvar05[[1]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=-44,y=5,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 5: CVaR10
plot(density(all.cvar10[[1]][5,]),col=0,ylim=c(0,8),xlim=c(-33,-30.5),
     main='CVaR10',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.cvar10[[1]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.cvar10[[1]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=-33,y=8,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 6: CO
plot(density(apply(all.cos[[1]], c(1,3), mean)[5,]),col=0,ylim=c(0,14000),xlim=c(0.46,0.47),
     main='Average concentration',xlab='',ylab='')
for(i in 1:length(models)) {
  zzz=apply(all.cos[[1]], c(1,3), mean)[i,]
  lines(density(zzz),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(apply(all.cos[[1]], c(1,3), mean)[i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=0.46,y=14000,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 7: TO
plot(density(apply(all.tos[[1]], c(1,3), mean,na.rm=TRUE)[5,]),col=0,ylim=c(0,2500),xlim=c(0.405,0.43),
     main='Average turnover',xlab='',ylab='')
for(i in 1:length(models)) {
  zzz=apply(all.tos[[1]], c(1,3), mean,na.rm=TRUE)[i,]
  lines(density(zzz),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(apply(all.tos[[1]], c(1,3), mean,na.rm=TRUE)[i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=0.405,y=2500,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)
dev.off()


##-------------
# CVAR05 
##-------------

pdf('tables_and_figures/portfolio_cis_cvar05.pdf',height=8,width=12)
par(mfrow=c(4,2),mar= c(2,2,2,2))
# plot 1: avrg return
plot(density(all.avrg[[2]][5,]),col=0,ylim=c(0,9),xlim=c(6.3,7.2),
     main='Average return',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.avrg[[2]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.avrg[[2]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=6.3,y=9,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 2: stdev
plot(density(all.sds[[2]][5,]),col=0,ylim=c(0,175),xlim=c(1.32,1.4),
     main='Standard deviation',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.sds[[2]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.sds[[2]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=1.32,y=175,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 3: SR
plot(density(all.SR[[2]][5,]),col=0,ylim=c(0,15),xlim=c(3.35,3.95),
     main='Adjusted Sharpe ratio',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.SR[[2]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.SR[[2]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=3.35,y=15,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 4: CVaR05
plot(density(all.cvar05[[2]][5,]),col=0,ylim=c(0,3.5),xlim=c(-44,-40),
     main='CVaR05',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.cvar05[[2]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.cvar05[[2]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=-44,y=3.5,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 5: CVaR10
plot(density(all.cvar10[[2]][5,]),col=0,ylim=c(0,4.2),xlim=c(-33,-30.5),
     main='CVaR10',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.cvar10[[2]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.cvar10[[2]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=-33,y=4.2,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 6: CO
plot(density(apply(all.cos[[2]], c(1,3), mean)[5,]),col=0,ylim=c(0,11000),xlim=c(0.462,0.471),
     main='Average concentration',xlab='',ylab='')
for(i in 1:length(models)) {
  zzz=apply(all.cos[[2]], c(1,3), mean)[i,]
  lines(density(zzz),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(apply(all.cos[[2]], c(1,3), mean)[i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=0.462,y=11000,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 7: TO
plot(density(apply(all.tos[[2]], c(1,3), mean,na.rm=TRUE)[5,]),col=0,ylim=c(0,1200),xlim=c(0.405,0.425),
     main='Average turnover',xlab='',ylab='')
for(i in 1:length(models)) {
  zzz=apply(all.tos[[2]], c(1,3), mean,na.rm=TRUE)[i,]
  lines(density(zzz),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(apply(all.tos[[2]], c(1,3), mean,na.rm=TRUE)[i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=0.405,y=1200,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)
dev.off()



##-------------
# CVAR10 
##-------------

pdf('tables_and_figures/portfolio_cis_cvar10.pdf',height=8,width=12)
par(mfrow=c(4,2),mar= c(2,2,2,2))
# plot 1: avrg return
plot(density(all.avrg[[3]][5,]),col=0,ylim=c(0,13),xlim=c(6.3,7.2),
     main='Average return',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.avrg[[3]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.avrg[[3]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=6.3,y=13,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 2: stdev
plot(density(all.sds[[3]][5,]),col=0,ylim=c(0,230),xlim=c(1.32,1.39),
     main='Standard deviation',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.sds[[3]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.sds[[3]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=1.32,y=230,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 3: SR
plot(density(all.SR[[3]][5,]),col=0,ylim=c(0,30),xlim=c(3.35,3.95),
     main='Adjusted Sharpe ratio',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.SR[[3]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.SR[[3]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=3.35,y=30,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 4: CVaR05
plot(density(all.cvar05[[3]][5,]),col=0,ylim=c(0,3),xlim=c(-44,-39.5),
     main='CVaR05',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.cvar05[[3]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.cvar05[[3]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=-44,y=3,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 5: CVaR10
plot(density(all.cvar10[[3]][5,]),col=0,ylim=c(0,5.2),xlim=c(-33,-30.5),
     main='CVaR10',xlab='',ylab='')
for(i in 1:length(models)) {
  lines(density(all.cvar10[[3]][i,]),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(all.cvar10[[3]][i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=-33,y=5.2,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 6: CO
plot(density(apply(all.cos[[3]], c(1,3), mean)[5,]),col=0,ylim=c(0,8000),xlim=c(0.462,0.471),
     main='Average concentration',xlab='',ylab='')
for(i in 1:length(models)) {
  zzz=apply(all.cos[[3]], c(1,3), mean)[i,]
  lines(density(zzz),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(apply(all.cos[[3]], c(1,3), mean)[i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=0.462,y=8000,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)

# plot 7: TO
plot(density(apply(all.tos[[3]], c(1,3), mean,na.rm=TRUE)[5,]),col=0,ylim=c(0,1400),xlim=c(0.405,0.425),
     main='Average turnover',xlab='',ylab='')
for(i in 1:length(models)) {
  zzz=apply(all.tos[[3]], c(1,3), mean,na.rm=TRUE)[i,]
  lines(density(zzz),col=all.colors[i],lty=all.ltys[i],lwd=2)
  points(x=mean(apply(all.tos[[3]], c(1,3), mean,na.rm=TRUE)[i,]),y=0,col=all.colors[i],pch=16)
}
legend(x=0.405,y=1400,col=all.colors,lty=all.ltys,lwd=2,legend=all.legends)
dev.off()

##----------
## Table
##----------

ALL.res = matrix(NA,nrow=21,ncol=length(models))

# 1-7 GMV portfolio
# 8-14 CVAR05
# 15-21 CVAR10

# 1 avrg_return 2 sd 3 SR 4 CVAR05 5 CVAR10 6 TO 7 CO 

for(i in 1:3){
  ALL.res[i*7-6,] = apply(all.avrg[[i]], 1, mean)
  ALL.res[i*7-5,] = apply(all.sds[[i]], 1, mean)
  ALL.res[i*7-4,] = apply(all.SR[[i]], 1, mean)
  ALL.res[i*7-3,] = apply(all.cvar05[[i]], 1, mean)
  ALL.res[i*7-2,] = apply(all.cvar10[[i]], 1, mean)
  ALL.res[i*7-1,] = apply(apply(all.cos[[3]], c(1,3), mean), 1, mean)
  ALL.res[i*7,]   = apply(apply(all.tos[[3]], c(1,3), mean,na.rm=TRUE), 1, mean)
}

colnames(ALL.res) = models
rownames(ALL.res) = rep(c('mean','sdev','A.Sh.','CVaR05','CVaR10','TO','CO'),3)


tableLines <- print(xtable(ALL.res,digits = 3,caption="Portfolio allocation results based on 1-step-ahead
              predictions for 2020/01/02 to 2023/01/31 out-of-sample period ($K$ = 797 observations).
              The three portfolios are:
             Global Minimum Variance (GMV), and minimum Conditional Value at Risk
             for lower 5 and 10 percentiles (CVaR05 and CVaR10), all with short-sale constraints.
             The table reports average portfolio return,
             overall standard deviation (in \\%),
             adjusted Sharpe ratio, Conditional Value at Risk
             for lower 5 and 10 percentiles,
             turnover, and concentration (all quantities annualized) for the
             pooled models (Geweke's, Jore's and
             equally weighted),
              two best individual models (CAW  and DCC-t) and a competitor model (DCC-HEAVY-t).",
                           align = "lccc|ccc",label='table:gmvfull_FX_new'),
                    scalebox=0.8,sanitize.text.function=function(x){x})
writeLines (tableLines, con = "tables_and_figures/gmvfull_FX_new.tex")


##-----------
## Which are best

# GMV
i=1
order(ALL.res[i*7-6,],decreasing = TRUE)
order(ALL.res[i*7-5,])
order(ALL.res[i*7-4,],decreasing = TRUE)
order(ALL.res[i*7-3,],decreasing = TRUE)
order(ALL.res[i*7-2,],decreasing = TRUE)
order(ALL.res[i*7-1,])
order(ALL.res[i*7,])

# CVAR05
i=2
order(ALL.res[i*7-6,],decreasing = TRUE)
order(ALL.res[i*7-5,])
order(ALL.res[i*7-4,],decreasing = TRUE)
order(ALL.res[i*7-3,],decreasing = TRUE)
order(ALL.res[i*7-2,],decreasing = TRUE)
order(ALL.res[i*7-1,])
order(ALL.res[i*7,])

# CVAR10
i=3
order(ALL.res[i*7-6,],decreasing = TRUE)
order(ALL.res[i*7-5,])
order(ALL.res[i*7-4,],decreasing = TRUE)
order(ALL.res[i*7-3,],decreasing = TRUE)
order(ALL.res[i*7-2,],decreasing = TRUE)
order(ALL.res[i*7-1,])
order(ALL.res[i*7,])
