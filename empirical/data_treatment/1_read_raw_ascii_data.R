rm(list=ls(all=TRUE))
library(xts)
library(highfrequency)
library(tidyquant)
library(moments)

min.seq = seq.POSIXt(as.POSIXct("2007-01-01 00:00:00",tz="UTC"), 
                     as.POSIXct("2023-12-31 23:59:59",tz="UTC"), 
                     by = "1 min")
all.dfs = list()

assets  = c('EURUSD','EURGBP','EURJPY','EURAUD','USDCAD',
            'USDJPY','GBPJPY','GBPUSD','AUDUSD','NZDJPY',
            'EURCAD','AUDCAD','CADJPY','GBPAUD','AUDNZD')

for(i in 1:length(assets)){
  x       = utils::read.table(paste('C:/Users/USER/Desktop/raw_histdata_data/',assets[i],
                             ".csv",sep=''),header = FALSE, sep = ";",dec=".")
  colnames(x) <- c("DateTime Stamp", "OPEN", "HIGH", "LOW", "CLOSE", "Volume")
  # format the first column as a date.time
  temp    <- as.POSIXct(strptime(x[,1], "%Y%m%d %H%M%S"),tz = "UTC")
  x$temp  <- temp+60*60*7 #change to GMT+2; now data is weekdays only
  # does not affect the co-dependence. Also, in line with the data
  # from Yahoo fiance
  index   = match(min.seq,temp)
  index   = na.omit(index)
  all.dfs[[i]] = xts(x[index,5], order.by=x[index,7])
}

names(all.dfs) = assets
all.xts        = do.call("merge", all.dfs)

save(all.xts,file='data/raw1min_xts.Rdata')



