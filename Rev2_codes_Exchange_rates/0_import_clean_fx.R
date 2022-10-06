rm(list=ls(all=TRUE))
# https://search.r-project.org/CRAN/refmans/highfrequency/html/rCov.html

library(xts)
library(highfrequency)
library(lubridate)
library(tseries)
library(moments)
library(nortsTest)

lret.daily = function(y){
  z=to.period(y,period = 'days')[,4]
  z2=as.vector(z)
  lret = c(NA,(log(z2[-1])-log(z[-length(z2)]))*100)
  res = xts(lret,order.by = as_date(index(z)))
  return(res)
}

rets = xts(order.by = seq(as.Date("2010/1/1"), as.Date("2023/1/1"), "days"))
min.seq = seq.POSIXt(as.POSIXct("2010-01-01 00:00:00",tz="UTC"), 
                     as.POSIXct("2021-12-31 23:59:59",tz="UTC"), 
                     by = "5 min")

years = seq(2010,2021)

##------------
## EUR-USD
##------------

eur.usd.full = NULL

for(j in 1:length(years)){
  x <- read.table(paste("data/DAT_ASCII_EURUSD_M1_",years[j],".csv",sep=''), header = FALSE, sep = ";", dec = ".")
  colnames(x) <- c("DateTime Stamp", "OPEN", "HIGH", "LOW", "CLOSE", "Volume")
  # format the first column as a date.time
  temp <- as.POSIXct(strptime(x[,1], "%Y%m%d %H%M%S"))
  x$temp <- temp
  index   = match(min.seq,temp)
  index   = na.omit(index)
  eur.usd.full = rbind(eur.usd.full,xts(x[index,5], order.by=x[index,7]))
}

eur.usd.lret = lret.daily(eur.usd.full)

plot(eur.usd.lret)

##------------
## USD-JPY
##------------

usd.jpy.full = NULL

for(j in 1:length(years)){
  x <- read.table(paste("data/DAT_ASCII_USDJPY_M1_",years[j],".csv",sep=''), header = FALSE, sep = ";", dec = ".")
  colnames(x) <- c("DateTime Stamp", "OPEN", "HIGH", "LOW", "CLOSE", "Volume")
  # format the first column as a date.time
  temp <- as.POSIXct(strptime(x[,1], "%Y%m%d %H%M%S"))
  x$temp <- temp
  index   = match(min.seq,temp)
  index   = na.omit(index)
  usd.jpy.full = rbind(usd.jpy.full,xts(x[index,5], order.by=x[index,7]))
}

usd.jpy.lret = lret.daily(usd.jpy.full)

plot(usd.jpy.lret)

##------------
## EUR/GBP
##------------

eur.gbp.full = NULL

for(j in 1:length(years)){
  x <- read.table(paste("data/DAT_ASCII_EURGBP_M1_",years[j],".csv",sep=''), header = FALSE, sep = ";", dec = ".")
  colnames(x) <- c("DateTime Stamp", "OPEN", "HIGH", "LOW", "CLOSE", "Volume")
  # format the first column as a date.time
  temp <- as.POSIXct(strptime(x[,1], "%Y%m%d %H%M%S"))
  x$temp <- temp
  index   = match(min.seq,temp)
  index   = na.omit(index)
  eur.gbp.full = rbind(eur.gbp.full,xts(x[index,5], order.by=x[index,7]))
}

eur.gbp.lret = lret.daily(eur.gbp.full)

plot(eur.gbp.lret)


##---------------
## Match the dates
##---------------

df = merge(eur.usd.full,usd.jpy.full,eur.gbp.full,all=TRUE,fill=NA)
df = na.omit(df)
dm = dim(df)[2]

#rm(eur.usd.full,usd.jpy.full,eur.gbp.full)

##---------------
## Extract RVs and RCovs
##---------------

# can change to 5, 10, 30 for different frequency returns

rv = rKernelCov(rData = df, alignBy = "minutes",
           alignPeriod = 10, makeReturns = TRUE,kernelType = "Bartlett")

RCov = array(unlist(rv)*100^2,c(dm,dm,length(rv)))
RCor = array(1,c(dm,dm,length(rv)))
RVs  = matrix(NA,ncol = dm,nrow=length(rv))

for(t in 1:length(rv)){
  RVs[t,]   = sqrt(diag(RCov[,,t]))
  RCor[,,t] = nearcor(round(cov2cor(RCov[,,t]),8))$cor
}

rets = merge(rets,eur.usd.lret,usd.jpy.lret,eur.gbp.lret,all=TRUE,fill=NA)
rets = na.omit(rets)
rets = rets[as.POSIXct(names(rv),tz="UTC")]
RCor = RCor[,,-1]
RCov = RCov[,,-1]
RVs  = RVs[-1,]

# Check here if everything is the same length, if all OK the proceed

dim(rets)
dim(RVs)
dim(RCor)
dim(RCov)

#---------------
# remove weekends, i.e. SUNDAY
#---------------

wkend = which(weekdays(index(rets))=="Sunday")
retsx = rets[-wkend,]
RVsx  = RVs[-wkend,]
RCorx = RCor[,,-wkend]
RCovx = RCov[,,-wkend] 

par(mfrow=c(3,1))
for(i in 1:dm) plot(RVsx[1:500,i],type='l')
 
dim(retsx)
dim(RVsx)
dim(RCorx)
dim(RCovx)

#---------------
# remove Holidays
#---------------

index(retsx)[1]
tail(index(retsx),1)

library(RQuantLib)
holid = getHolidayList("UnitedStates", index(retsx)[1], tail(index(retsx),1))

z=(index(retsx))
holidind = NULL
for(i in 1:length(holid)){
  if(length(which(z==holid[i]))==1){
    holidind = c(holidind,which(z==holid[i]))
  } 
}

retsx = retsx[-holidind,]
RVsx  = RVsx[-holidind,]
RCorx = RCorx[,,-holidind]
RCovx = RCovx[,,-holidind] 

par(mfrow=c(3,1))
for(i in 1:dm) plot(RVsx[1:500,i],type='l')

dim(retsx)
dim(RVsx)
dim(RCorx)
dim(RCovx)

#---------------
# Exclude nearly singular matrices; cannot be used in estimation
#---------------

determinants = rep(NA,dim(retsx)[1])

for(t in 1:dim(retsx)[1]){
  determinants[t] = det(RCorx[,,t])
}

par(mfrow=c(1,2))
hist(determinants,breaks = 100)
exclude = which(determinants<0.01)
hist(determinants[-exclude],breaks = 100)


rets = retsx[-exclude]
par(mfrow=c(3,1))
plot(rets[,1])
plot(rets[,2])
plot(rets[,3])

RCov = RCovx[,,-exclude]
RCor = RCorx[,,-exclude]
RVs  = RVsx[-exclude,]

save(rets,RCov,RCor,RVs,file='Rev2_codes_Exchange_rates/EX.Rdata')