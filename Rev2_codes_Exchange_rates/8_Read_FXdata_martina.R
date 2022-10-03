#setwd("C:/Users/zahar/OneDrive - Colegio Universitario de Estudios Financieros (CUNEF)/Documentos/R/DMFP/data")
rm(list=ls(all=TRUE))

setwd("~/R/DMFP")


# read data EUR/USD for 2020 and headers according web site
eur.usd.2020 <- read.table("data/DAT_ASCII_EURUSD_M1_2020.csv", header = FALSE, sep = ";", dec = ".")
colnames(eur.usd.2020) <- c("DateTime Stamp", "Bar OPEN Bid Quote", "Bar HIGH Bid Quote", "Bar LOW Bid Quote", "Bar CLOSE Bid Quote", "Volume")

# format the first column as a date.time
temp <- strptime(eur.usd.2020[,1], "%Y%m%d %H%M%S")
eur.usd.2020 <- cbind(temp,eur.usd.2020[,-1])

# change header to shorter title
colnames(eur.usd.2020) <- c("DateTime", "OPEN", "HIGH", "LOW", "CLOSE", "Volume")

# calculate returns
logp <- log(eur.usd.2020$CLOSE) 
ret <- log(eur.usd.2020$CLOSE)[-1] - log(eur.usd.2020$CLOSE)[-length(logp)]
plot(ret, type = 'l')

# the vector difft gives the time gap associated to every return in minutes
difft <- difftime(head(eur.usd.2020[,1],-1),tail(eur.usd.2020[,1],-1))
# you can table it to see the absolute frequency of every possible time gap (it´s supposed to be 1 minute returns, but it isn´t always)
table(-difft)