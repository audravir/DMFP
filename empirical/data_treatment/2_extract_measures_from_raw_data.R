rm(list=ls(all=TRUE))
library(xts)
library(highfrequency)
library(tidyquant)
library(moments)
library(sfsmisc)
library(lubridate)
library(plyr)

load('data/raw1min_xts.Rdata')
assets  = names(all.xts)
dm      = ncol(all.xts)

# remove saturdays and sundays

wkdays  = weekdays(as.Date(index(all.xts)))
table(wkdays)

# > table(wkdays)
# wkdays
# Friday    Monday  Saturday    Sunday  Thursday   Tuesday Wednesday 
# 1246804   1247113        56       958   1251791   1249617   1251259 

sat.ind = which(wkdays=="Saturday")
sun.ind = which(wkdays=="Sunday")

all.xts = all.xts[-c(sat.ind,sun.ind),]

# # remove X-mas and January 1st
# FX trading is de-centralized and not affected by national holidays
# but Dec 25th and Jan 1s is not open anywhere

Year    = seq(2007,2022)
xMonth  = rep('12',length(Year))
xDay    = rep('25',length(Year))
nyMonth = rep('01',length(Year))
nyDay   = rep('01',length(Year))

xdates  = as.Date(paste(Year, xMonth, xDay, sep = '-'))
nydates = as.Date(paste(Year, nyMonth, nyDay, sep = '-'))
xdates
nydates

dayindex    = c(xdates,nydates)
indexasdate = as.Date(index(all.xts))
holi.days   = NULL

for(i in 1:length(dayindex)){
  ind  = which(dayindex[i]==indexasdate)
  holi.days = c(holi.days,ind)
}

unique(as.Date(index(all.xts[holi.days,])))

all.xts = all.xts[-holi.days,]

# start at '2008-01-02'

which(date(index(all.xts))=='2008-01-02')[1]
dim(all.xts)[1]
dim(all.xts)[1]-which(date(index(all.xts))=='2008-01-02')[1]

all.xts = all.xts[(which(date(index(all.xts))=='2008-01-02')[1]):(dim(all.xts)[1]),]
dim(all.xts)


# delete such days where at least one asset has no data
ep         = endpoints(all.xts,on="days")
all.na     = function(x){sum(is.na(x))==length(x)}
empty.days = NULL

for(j in 1:dm){
  empty.days = cbind(empty.days,period.apply(all.xts[,j],INDEX=ep,FUN=all.na))
}

apply(empty.days,2,sum)
# how many empty days per asset
# > apply(empty.days,2,sum)
# EURUSD EURGBP EURJPY EURAUD USDCAD USDJPY GBPJPY GBPUSD AUDUSD 
# 0      0      0      5      0      0      2      0      0 
# NZDJPY EURCAD AUDCAD CADJPY GBPAUD AUDNZD 
# 0      0      2      0      0      0 

empty.index = which(apply(empty.days,1,sum)!=0)

# > names(empty.index)
# [1] "2008-06-16 23:59:00" "2008-06-17 23:59:00"
# [3] "2008-06-18 23:59:00" "2008-06-19 23:59:00"
# [5] "2008-06-20 23:59:00" "2009-02-19 23:59:00"
# [7] "2009-02-20 23:59:00" "2009-03-12 23:59:00"
# [9] "2009-03-13 23:59:00"


# check one as example

apply(all.xts[which(as.Date(index(all.xts))=="2008-06-16"),],2,all.na)

# remove the empty days
indexasdate = as.Date(index(all.xts))
ind         = as.Date(index(empty.days[empty.index,]))
remove      = NULL

for(i in 1:length(ind)){
  remove=c(remove,which(indexasdate==ind[i]))
}

all.xts = all.xts[-remove,]

# there are some missing hours, fill-in all NA's with the previous observation
all.xts.original = all.xts
all.xts          = na.locf(all.xts)

## Remove days where we have only a few hours of trading

dayindex = unique(as.Date(index(all.xts)))
hrs      = rep(NA,length(dayindex))
indexasdate = as.Date(index(all.xts))
fewhrs      = NULL

for(i in 1:length(dayindex)){
  ind  = which(dayindex[i]==indexasdate)
  temp = index(all.xts[ind,])
  hrs[i] = as.numeric(difftime(tail(temp,1),head(temp,1),units="hours"))
  if(hrs[i]<23){
    fewhrs = c(fewhrs,ind)
  } 
}

par(mfrow=c(1,1))
plot(hrs,pch=20);abline(h=23,col=2)
sum((hrs<23))
unique((as.Date(index(all.xts[fewhrs,]))))

all.xts = all.xts[-fewhrs,]


## calculate daily log returns
# open-to-close returns

hfrets  = makeReturns(all.xts)*100
ep      = endpoints(all.xts,on='days')
rets    = NULL
sum_f1  = function(x)sum(x[-1]) # exclude 1st obs to avoid overnight

for(i in 1:dm){
  rets = cbind(rets,period.apply(hfrets[,i],ep,sum_f1))
}
rets = xts(coredata(rets),order.by = as.Date(index(rets)))
dim(rets)  

rt = coredata(rets)
par(mfrow=c(1,1))
for(i in 1:dm){
  plot(rt[,i],main =i,type='l')
}




ass=8

rets[which.max(rt[,ass]),ass]



## extract daily realized covariances and correlations

rv = rKernelCov(rData = all.xts, alignBy = "minutes",
           alignPeriod = 1, makeReturns = TRUE)

length(rv)

RCov = array(unlist(rv)*100^2,c(dm,dm,length(rv)))
RCor = array(1,c(dm,dm,length(rv)))
RVs  = matrix(NA,ncol = dm,nrow=length(rv))

for(t in 1:length(rv)){
  RVs[t,]   = sqrt(diag(RCov[,,t]))
  RCor[,,t] = nearcor(round(cov2cor(RCov[,,t]),8))$cor
}

# Check here if everything is the same length, if all OK then proceed

dim(rets)
dim(RVs)
dim(RCor)
dim(RCov)

stand = stand.full = matrix(NA,ncol=dm,nrow=dim(rets)[1])

for(t in 1:dim(rets)[1]){
  stand[t,] = rets[t,]/sqrt(diag(RCov[,,t]))
}

stand_means = apply(stand, 2, mean)
stand_sds   = apply(stand, 2, sd)

marginals = list(stand_means,stand_sds)
names(marginals) = c('mean','sd')
save(marginals,file='temp/marginals.Rdata')

for(i in 1:dm){
  stand.full[,i] = (stand[,i]-stand_means[i])/stand_sds[i]
}

apply(stand.full,2,mean)
apply(stand.full,2,sd)


##

Sigma = alply(RCor,3) #
date  = as.Date(index(rets))
rets  = coredata(rets)
stand = stand.full
udata = pnorm(stand.full)


save(Sigma,rets,RVs,stand,udata,RCov,RCor,assets,date,file='data/FXdata.Rdata')

