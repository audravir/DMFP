rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')

library(matrixcalc)
library(mixAK)
library(LaplacesDemon)
library(countreg)

nn       = length(date)
#end.date = which(zoo::as.yearmon(date)=="jan 2021")[1]-1
end.date = which(zoo::as.yearmon(date)=="ene 2021")[1]-1
if(is.na(date[end.date])){stop("Change the name of the month of end.date to system language")}


T0  = end.date  #
T0 # for estimation
K = nn-end.date
K # for oos evaluation

T0+K
nn
# the same

##???????????
data = Sigma[1:T0]


xm1(Sigma,200)

