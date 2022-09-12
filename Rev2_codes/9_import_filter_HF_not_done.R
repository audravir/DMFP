rm(list=ls(all=TRUE))

eur.usd.2020 <- read.table("data/DAT_ASCII_EURUSD_M1_2020.csv", header = FALSE, sep = ";", dec = ".")
colnames(eur.usd.2020) <- c("DateTime Stamp", "Bar OPEN Bid Quote", "Bar HIGH Bid Quote", "Bar LOW Bid Quote", "Bar CLOSE Bid Quote", "Volume")

# format the first column as a date.time
temp <- as.POSIXct(strptime(eur.usd.2020[,1], "%Y%m%d %H%M%S"))
eur.usd.2020$temp <- temp

min.seq = seq.POSIXt(temp[1], tail(temp,1), by = "10 min")
index   = match(min.seq,temp)
index   = na.omit(index)

reduced = eur.usd.2020[index,c(5,7)]
rm(eur.usd.2020)


