rm(list=ls(all=TRUE))

library(forecast)

load('data/10data.Rdata')
oos = 252
nn  = dim(RVs)[1]
dm  = dim(RVs)[2]

#-----------------------------------
# Produced by HAR(1,5,35) model
# RVs are standard deviations!
#-----------------------------------

load('temp/RV_forc2.Rdata')
par(mfrow=c(2,5))
for(i in 1:dm){
  plot(tail(RVs[,i],oos),type='l',lwd=2)
  lines(RV_forc$`1sa`[,i],col=2)
}


###########
# Forecast via ARFIMA
###########


library(fracdiff)

Forcs = matrix(NA,nrow=oos,ncol=dm)

for(i in 1:dm){
  for(t in 1:oos){
    fit = fracdiff::fracdiff(log(RVs[1:(nn-oos+t-1),i]^2),nar=1)
    Forcs[t,i] = exp(forecast(fit,h = 1)[[2]][1]/2)
  }
}

par(mfrow=c(2,5))
for(i in 1:dm){
  plot(tail(RVs[,i],oos),type='l',lwd=2)
  lines(RV_forc$`1sa`[,i],col=2)
  lines(Forcs[,i],col=3)
}

##

RVsforc1_m1 = RV_forc$`1sa`
RVsforc1_m2 = Forcs

RV_forc = list(tail(date,oos),
               Forcs)
names(RV_forc) = c('Date','1sa')
save(RV_forc,file='temp/RV_forc3.Rdata')



ONEsa_mae  = matrix(NA,ncol=dm,nrow=2)
ONEsa_rmse = matrix(NA,ncol=dm,nrow=2)
#FIVEsa  = matrix(NA,ncol=dm,nrow=4)

em1ONEsa = em2ONEsa = em3ONEsa = em4ONEsa = matrix(NA,nrow=oos,ncol=dm)
#em1FIVEsa = em2FIVEsa = matrix(NA,nrow=oos,ncol=dm)


for(i in 1:dm){
  em1ONEsa[,i] = RVsforc1_m1[,i] - tail(RVs[,i],oos)
  em2ONEsa[,i] = RVsforc1_m2[,i] - tail(RVs[,i],oos)
  
  
  ONEsa_rmse[1,i] = sqrt(mean((em1ONEsa[,i])^2))
  ONEsa_rmse[2,i] = sqrt(mean((em2ONEsa[,i])^2))
  
  
  ONEsa_mae[1,i] = mean(abs(em1ONEsa[,i]))
  ONEsa_mae[2,i] = mean(abs(em2ONEsa[,i]))

}


ONEsa = rbind(ONEsa_rmse/matrix(apply(ONEsa_rmse, 2, min),ncol=dm,nrow=2,byrow=TRUE)
              ,ONEsa_mae/matrix(apply(ONEsa_mae, 2, min),ncol=dm,nrow=2,byrow=TRUE)
)

ONEsa

# library(MCS)
# # https://www.rdocumentation.org/packages/MCS/versions/0.1.3/topics/LossVol
# 
# lossf= "AE1"
# 
# mean(LossVol(tail(RVs[,i],oos), RVsforc1_m1[,i], which = lossf))
# mean(LossVol(tail(RVs[,i],oos), RVsforc1_m2[,i], which = lossf))

