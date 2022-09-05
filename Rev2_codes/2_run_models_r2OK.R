rm(list=ls(all=TRUE))
library(sfsmisc)
library(truncnorm)
library(DMFP)

load('data/10data.Rdata')
T0  = 1990 

# the length of the estimation window, 1:T0, which is 
# from "2001-02-01" till "2008-12-31"
# the rest will be used for model evaluation

# 250 is MCMC size. The codes have been run for 50k and the results are in
# Rdata files in folder temp. The files are too large for GitHub
# they are stored locally 
# 25k for burn-in+25k 
# thinned every 25th, the resulting posterior samples of 1000k

scalartdcc(stand,250)
rmetrics(stand,250)
scalardcc(stand,250)
xm1(Sigma,250)

