rm(list=ls(all=TRUE))
library(sfsmisc)
library(truncnorm)

load('data/10data.Rdata')
T0  = 1990 

# the length of the estimation window, 1:T0, which is 
# from "2001-02-01" till "2008-12-31"
# the rest will be used for model evaluation

scalartdcc(stand,250)
rmetrics(stand,250)
scalardcc(stand,250)
xm1(Sigma,250)

