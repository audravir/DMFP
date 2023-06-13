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

min.seq = seq.POSIXt(as.POSIXct("2010-01-01 00:00:00",tz="UTC"), 
                     as.POSIXct("2021-12-31 23:59:59",tz="UTC"), 
                     by = "5 min")

years = seq(2012,2021)

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

head(eur.usd.full)

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



##------------
## sp500
##------------

sp500.full = NULL

for(j in 1:length(years)){
  x <- read.table(paste("data/DAT_ASCII_SPXUSD_M1_",years[j],".csv",sep=''), header = FALSE, sep = ";", dec = ".")
  colnames(x) <- c("DateTime Stamp", "OPEN", "HIGH", "LOW", "CLOSE", "Volume")
  # format the first column as a date.time
  temp <- as.POSIXct(strptime(x[,1], "%Y%m%d %H%M%S"))
  x$temp <- temp
  index   = match(min.seq,temp)
  index   = na.omit(index)
  sp500.full = rbind(sp500.full,xts(x[index,5], order.by=x[index,7]))
}

##------------
## wti
##------------

wti.full = NULL

for(j in 1:length(years)){
  x <- read.table(paste("data/DAT_ASCII_WTIUSD_M1_",years[j],".csv",sep=''), header = FALSE, sep = ";", dec = ".")
  colnames(x) <- c("DateTime Stamp", "OPEN", "HIGH", "LOW", "CLOSE", "Volume")
  # format the first column as a date.time
  temp <- as.POSIXct(strptime(x[,1], "%Y%m%d %H%M%S"))
  x$temp <- temp
  index   = match(min.seq,temp)
  index   = na.omit(index)
  wti.full = rbind(wti.full,xts(x[index,5], order.by=x[index,7]))
}

##---------------
## Match the dates
##---------------

df = merge(eur.usd.full,usd.jpy.full,eur.gbp.full,
           sp500.full,wti.full,all=TRUE,fill=NA)
df = na.omit(df)
dm = dim(df)[2]
dim(df)

##---------------
## Extract RVs and RCovs
##---------------

rets = xts(order.by = seq(as.Date("2010/1/1"), as.Date("2023/1/1"), "days"))

for(i in 1:dm){
  tmp = lret.daily(df[,i])
  rets = merge(rets,tmp,all=TRUE,fill=NA)
}
rets = na.omit(rets)
head(rets)
tail(rets)

# can change to 5, 10, 30 for different frequency returns

rv = rKernelCov(rData = df, alignBy = "minutes",
           alignPeriod = 5, makeReturns = TRUE)

head(rv,3)
tail(rv,3)

## Use 5-minute returns! makes the estimation much more precise, almost no
# close to singular matrices

RCov = array(unlist(rv)*100^2,c(dm,dm,length(rv)))
RCor = array(1,c(dm,dm,length(rv)))
RVs  = matrix(NA,ncol = dm,nrow=length(rv))

for(t in 1:length(rv)){
  RVs[t,]   = sqrt(diag(RCov[,,t]))
  RCor[,,t] = nearcor(round(cov2cor(RCov[,,t]),8))$cor
}

# rets = xts(order.by = seq(as.Date("2010/1/1"), as.Date("2023/1/1"), "days"))
# rets = merge(rets,eur.usd.lret,usd.jpy.lret,eur.gbp.lret,
#              sp500.lret,wti.lret,all=TRUE,fill=NA)
# rets = na.omit(rets)
RCor = RCor[,,-1]
RCov = RCov[,,-1]
RVs  = RVs[-1,]

# Check here if everything is the same length, if all OK then proceed

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

par(mfrow=c(2,3))
for(i in 1:dm) plot(RVsx[1:500,i],type='l')
 
dim(retsx)
dim(RVsx)
dim(RCorx)
dim(RCovx)

#---------------
# remove Holidays
#---------------


library(RQuantLib)
holid = isHoliday("UnitedStates/NYSE", index(retsx))  

retsx = retsx[!holid,]
RVsx  = RVsx[!holid,]
RCorx = RCorx[,,!holid]
RCovx = RCovx[,,!holid] 

par(mfrow=c(3,1))
for(i in 1:dm) plot(RVsx[1:500,i],type='l')

dim(retsx)
dim(RVsx)
dim(RCorx)
dim(RCovx)

#---------------
# Exclude nearly singular matrices; cannot be used in estimation
#---------------

determinants.cov = rep(NA,dim(retsx)[1])

for(t in 1:dim(retsx)[1]){
  determinants.cov[t] = det(RCovx[,,t])
}

exclude = which(determinants.cov<1e-7)

rets = retsx[-exclude]
RCov = RCovx[,,-exclude]
RCor = RCorx[,,-exclude]
RVs  = RVsx[-exclude,]

stand =  matrix(NA,ncol=dm,nrow=dim(rets)[1])

for(t in 1:dim(rets)[1]){
  stand[t,] = rets[t,]/sqrt(diag(RCov[,,t]))
}

apply(stand,2,mean)
apply(stand,2,sd)

par(mfrow=c(dm,1))
for(i in 1:dm) plot(stand[,i],type='l')
for(i in 1:dm) acf(stand[,i]^2,ylim=c(-0.1,0.1))
for(i in 1:dm) pacf(stand[,i]^2)


save(rets,RCov,RCor,RVs,file='Rev2_codes_Exchange_rates/EX.Rdata')
