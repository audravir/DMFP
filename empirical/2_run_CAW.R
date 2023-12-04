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

c(5000,0.0001,0.01)
caw(data,M=5000,propsdb = 0.0001,propsdnu = 0.01)