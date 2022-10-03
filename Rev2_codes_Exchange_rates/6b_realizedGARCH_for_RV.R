rm(list=ls(all=TRUE))

load('data/10data.Rdata')
oos = 252
nn  = dim(RVs)[1]
dm  = dim(RVs)[2]

#-----------------------------------
# 
# RVs are standard deviations!
#-----------------------------------

load('temp/RV_forc2.Rdata')
par(mfrow=c(2,5))
for(i in 1:dm){
  plot(tail(RVs[,i],oos),type='l',lwd=2)
  lines(RV_forc$`1sa`[,i],col=2)
}


###########
# Forecast via realized GARCH
###########
# https://www.r-bloggers.com/2014/01/the-realized-garch-model/

library(rugarch)
require(xts)

Forcs = NULL

spec = ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
      variance.model = list(model = 'realGARCH', garchOrder = c(1, 1)))


for(i in 1:dm){
  rvdep = xts(RVs[,i]^2, order.by = date)
  retdep = xts(rets[,i], order.by = date)
  
  fit = ugarchfit(spec, retdep, out.sample = 252, solver = 'hybrid',
                  realizedVol =rvdep)
  
  forc1  = ugarchforecast(fit, n.ahead = 1, n.roll = 252)
  forcs1 = xts(sigma(forc1)[1, ], move(as.Date(names(sigma(forc1)[1, ])), by = 1))
  Forcs = cbind(Forcs,forcs1)
}

Forcs=as.matrix(Forcs[-oos,])


par(mfrow=c(2,5))
for(i in 1:dm){
  plot(tail((RVs[,i]),oos),type='l')
  lines(RV_forc$`1sa`[,i],lwd=2,col=2)
  lines(Forcs[,i],col=3,lwd=2)
}



em1ONEsa = em2ONEsa =  matrix(NA,nrow=oos,ncol=dm)
ONEsa_mae  = matrix(NA,ncol=dm,nrow=2)
ONEsa_rmse = matrix(NA,ncol=dm,nrow=2)

for(i in 1:dm){
  em1ONEsa[,i] = RV_forc$`1sa`[,i] - tail(RVs[,i],oos)
  em2ONEsa[,i] = Forcs[,i] - tail(RVs[,i],oos)

  ONEsa_rmse[1,i] = sqrt(mean((em1ONEsa[,i])^2))
  ONEsa_rmse[2,i] = sqrt(mean((em2ONEsa[,i])^2))

  ONEsa_mae[1,i] = mean(abs(em1ONEsa[,i]))
  ONEsa_mae[2,i] = mean(abs(em2ONEsa[,i]))
}

ONEsa = rbind(ONEsa_rmse/matrix(apply(ONEsa_rmse, 2, min),ncol=dm,nrow=2,byrow=TRUE)
              ,ONEsa_mae/matrix(apply(ONEsa_mae, 2, min),ncol=dm,nrow=2,byrow=TRUE)
)

ONEsa




RV_forc = list(tail(date,oos),Forcs)
names(RV_forc) = c('Date','1sa')

save(RV_forc,file='temp/RV_forc4.Rdata')

