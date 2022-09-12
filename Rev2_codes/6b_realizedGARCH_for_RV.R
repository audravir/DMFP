rm(list=ls(all=TRUE))

load('data/10data.Rdata')
oos = 252
nn  = dim(RVs)[1]
dm  = dim(RVs)[2]

#-----------------------------------
# Produced by HAR(1,5,22) model
# RVs are standard deviations!
#-----------------------------------

load('temp/RV_forc.Rdata')
par(mfrow=c(2,5))
for(i in 1:dm){
  plot(tail(RVs[,i],oos),type='l',lwd=2)
  lines(RV_forc$`1sa`[,i],col=2)
}


###########
# Forecast via realized GARCH
###########
# https://www.r-bloggers.com/2014/01/the-realized-garch-model/
