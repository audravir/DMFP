rm(list=ls(all=TRUE))
load('data/FXdata.Rdata')

library(xtable)

post.sample = 1000
dm          = dim(rets)[2]

##--------------
## rank-1-DCC-t
##--------------

load('empirical/temp/results_vectordcc_tcop.Rdata')
M   = dim(res$r)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

dim(res$r)

pdf(file='tables_and_figures/parameters_dcct.pdf',width=12,height=8)
par(mfrow=c(6, 5), mar=c(2, 2, 1, 1))
for(i in 1:dim(res$r)[2]) plot(res$r[ind,i],type='l')
dev.off()

##--------------
## DCC-HEAVY-t
##--------------

load('empirical/temp/results_heavy_t.Rdata')
M   = dim(res$r)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

dim(res$r)

pdf(file='tables_and_figures/parameters_heavyt.pdf',width=10,height=2)
par(mfrow=c(1, 3), mar=c(2, 2, 1, 1))
for(i in 1:dim(res$r)[2]) plot(res$r[ind,i],type='l')
dev.off()

##--------------
## AIW
##--------------

load('empirical/temp/results_xm.Rdata')
M   = dim(res$r)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

dim(res$r)

pdf(file='tables_and_figures/parameters_xm.pdf',width=12,height=8)
par(mfrow=c(6, 5), mar=c(2, 2, 1, 1))
for(i in 1:dim(res$r)[2]) plot(res$r[ind,i],type='l')
dev.off()


##--------------
## CAIW
##--------------

load('empirical/temp/results_caw_iw.Rdata')
M   = dim(res$r)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

dim(res$r)

pdf(file='tables_and_figures/parameters_caiw.pdf',width=12,height=8)
par(mfrow=c(6, 5), mar=c(2, 2, 1, 1))
for(i in 1:dim(res$r)[2]) plot(res$r[ind,i],type='l')
dev.off()
