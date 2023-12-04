rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')

library(matrixcalc)
library(mixAK)
library(LaplacesDemon)
library(countreg)
library(DMFP)

nn       = length(date)
end.date = which(zoo::as.yearmon(date)=="ene 2021")[1]-1
if(is.na(date[end.date])){end.date = which(zoo::as.yearmon(date)=="jan 2021")[1]-1}

date[end.date]

T0  = end.date  #
T0 # for estimation
K = nn-end.date
K # for oos evaluation

T0+K
nn
# the same

data = stand[1:T0,]

# 0.0005 gave accp 0.0034
# 0.00005 gave accp 0.786
# 0.000275 gave accp 0.0006

c(25000,0.00005)

vectordcc(data,10000,0.00005)

