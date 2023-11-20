rm(list=ls(all=TRUE))
library(xts)
library(highfrequency)
library(tidyquant)
library(moments)

sumst = function(x){
  vec = c(mean(x,na.rm = TRUE),sd(x,na.rm = TRUE),
          skewness(x,na.rm = TRUE),kurtosis(x,na.rm = TRUE))
  names(vec) = c('mean','sd','skew','kurt')
  return(round(vec,2))
}

min.seq = seq.POSIXt(as.POSIXct("2007-01-01 00:00:00",tz="UTC"), 
                     as.POSIXct("2023-12-31 23:59:59",tz="UTC"), 
                     by = "30 min")
all.dfs = list()

assets  = c('EURUSD','EURCHF','EURGBP','EURJPY','EURAUD',
            'USDCAD','USDCHF','USDJPY','GBPCHF','GBPJPY',
            'GBPUSD','AUDUSD','CHFJPY','NZDJPY','NZDUSD',
            'EURCAD','AUDCAD','CADJPY','GBPAUD','AUDNZD',
            'GBPCAD')

# EURCHF has a sharp drop on January 15th, 2015
# it is NOT an error in the data
# On Jan. 15, 2015, Switzerland announced that it was 
# going to scrap its currency peg of 1.20 to the euro. 
# The Swiss franc immediately skyrocketed 20%. 
# https://www.forex.in.rs/usd-chf-15-january-2015/

# select any asset

i=10

# for(i in 1:length(assets)){

x       = utils::read.table(paste('C:/Users/USER/Desktop/DMFP/raw_histdata_data/',assets[i],
                           ".csv",sep=''),header = FALSE, sep = ";",dec=".")
colnames(x) <- c("DateTime Stamp", "OPEN", "HIGH", "LOW", "CLOSE", "Volume")
# format the first column as a date.time
temp    <- as.POSIXct(strptime(x[,1], "%Y%m%d %H%M%S"),tz = 'UTC')
x$temp  <- temp+60*60*7 #change to GMT+2; now data is weekdays only
index   = match(min.seq,temp)
index   = stats::na.omit(index)
all.dfs[[i]] = xts::xts(x[index,5], order.by=x[index,7])
# }


hfrets   = makeReturns(all.dfs[[i]])
dailyret = aggregateTS(hfrets,FUN = "sum",alignBy = "days")

newts1 = xts(coredata(dailyret)*100,order.by = date(dailyret))

# compare with data from YahooFinance
data = getSymbols(paste(assets[i],'=X',sep=''), from = date(dailyret)[1],
                  to = tail(date(dailyret),1),warnings = FALSE,
                  auto.assign = FALSE,periodicity = "daily",
                  return.class = 'xts')

yf = diff(log(data[,6]))*100
newts2 = xts(coredata(yf),order.by = date(yf))

plot(newts1)
lines(newts2,col=2,lwd=2)

cbind(sumst((dailyret)*100),sumst(yf))

# check more closely some small subsample: 
plot(newts1[1:150])
lines(newts2[1:150],col=2,lwd=2)

# very similar.

