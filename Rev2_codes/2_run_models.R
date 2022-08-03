rm(list=ls(all=TRUE))
library(sfsmisc)
library(truncnorm)

load('10data.Rdata')
T0  = 1990 

# the length of the estimation window, 1:T0, which is 
# from "2001-02-01" till "2008-12-31"
# the rest will be used for model evaluation

source('FUN_scalardcc_tcop.R')
scalartdcc(stand,25000)

source('FUN_rmetrics.R')
rmetrics(stand,25000)

source('FUN_scalardcc.R')
scalardcc(stand,25000)

source('FUN_xm.R')
xm1(Sigma,25000)

