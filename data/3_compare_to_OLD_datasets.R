rm(list=ls(all=TRUE))

load('FXdata.Rdata')
fx15 = list(rets,RVs,stand,RCor,RCov)
rm(rets,RVs,stand,RCor,RCov)

load('old_5data_FX.Rdata')
fx5 = list(rets,RVs,stand,RCor,RCov)
rm(rets,RVs,stand,RCor,RCov)

load('old_10data.Rdata')
var10 = list(rets,RVs,stand,RCor,RCov)

rm(list=setdiff(ls(), c("fx15" , "fx5","var10")))



