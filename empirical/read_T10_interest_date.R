rm(list=ls(all=TRUE))
load('empirical/temp/res_FX.Rdata')
load('data/FXdata.Rdata')
library(zoo)

K
tail(date,K)
class(date)

DGS10$date = as.Date(DGS10$DATE, "%Y-%m-%d")
DGS10$rf   = as.numeric(DGS10$DGS10)
DGS10$rf   = na.locf(DGS10$rf)

ind = DGS10$date %in% tail(date,K)
rf  = DGS10$rf[ind]
mean(rf)
plot(tail(date,K),rf,type='l')

save(rf,file='data/rf.RData')
