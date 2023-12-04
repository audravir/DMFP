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

data = Sigma[1:T0]
data = Sigma[1:1000]

c(2000,0.001,0.001)
xm(data,M=2000,propsdb=0.001,propsdnu=0.001)

