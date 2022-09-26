rm(list=ls(all=TRUE))
library(xts)
library(realized)

##------------
## EUR-USD
##------------

eur.usd.2020 <- read.table("data/DAT_ASCII_EURUSD_M1_2020.csv", header = FALSE, sep = ";", dec = ".")
colnames(eur.usd.2020) <- c("DateTime Stamp", "Bar OPEN Bid Quote", "Bar HIGH Bid Quote", "Bar LOW Bid Quote", "Bar CLOSE Bid Quote", "Volume")

# format the first column as a date.time
temp <- as.POSIXct(strptime(eur.usd.2020[,1], "%Y%m%d %H%M%S"))
eur.usd.2020$temp <- temp

min.seq = seq.POSIXt(temp[1], tail(temp,1), by = "10 min")
index   = match(min.seq,temp)
index   = na.omit(index)

reduced1 = xts(eur.usd.2020[index,5], order.by=eur.usd.2020[index,7])
rm(eur.usd.2020)

##------------
## EUR-USD
##------------

usd.jpy.2020 <- read.table("data/DAT_ASCII_USDJPY_M1_2020.csv", header = FALSE, sep = ";", dec = ".")
colnames(usd.jpy.2020) <- c("DateTime Stamp", "Bar OPEN Bid Quote", "Bar HIGH Bid Quote", "Bar LOW Bid Quote", "Bar CLOSE Bid Quote", "Volume")

# format the first column as a date.time
temp <- as.POSIXct(strptime(usd.jpy.2020[,1], "%Y%m%d %H%M%S"))
usd.jpy.2020$temp <- temp

min.seq = seq.POSIXt(temp[1], tail(temp,1), by = "10 min")
index   = match(min.seq,temp)
index   = na.omit(index)

reduced2 = xts(usd.jpy.2020[index,5], order.by=usd.jpy.2020[index,7])
rm(usd.jpy.2020)


length(reduced1)
length(reduced2)

plot(reduced2)





