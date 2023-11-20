rm(list=ls(all=TRUE))
library(sfsmisc)
library(truncnorm)
library(DMFP)
library(zoo)

load('data/3data_EX.Rdata')

nn       = length(date)
end.date = which(zoo::as.yearmon(date)=='Jan 2021')[1]-1
date[end.date]

T0  = end.date  #
T0
K = nn-end.date
K

T0+K
nn
# the length of the estimation window, 1:T0, which is 

# 250 is MCMC size. The codes have been run for 50k and the results are in
# Rdata files in folder temp. The files are too large for GitHub
# they are stored locally 
# 25k for burn-in+25k 
# thinned every 25th, the resulting posterior samples of 1000k

scalartdcc(stand,15000)
rmetrics(stand,15000)
scalardcc(stand,15000)
xm1(Sigma,15000)

