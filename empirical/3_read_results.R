rm(list=ls(all=TRUE))
load('data/FXdata.Rdata')

library(xtable)
library(coda)
p25 = function(x) quantile(x,0.25)
p75 = function(x) quantile(x,0.75)


post.sample = 1000
dm          = dim(rets)[2]
ess         = function(x) effectiveSize(mcmc(x))/length(x)

##--------------
## rank-1-DCC-t
##--------------

load('empirical/temp/results_vectordcc_tcop.Rdata')
M   = dim(res$r)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

dim(res$r)

par.table = rbind(apply(res$r,2,p25),
                  apply(res$r,2,median),
                  apply(res$r,2,p75),
                  c(mean(res$accnu),rep(mean(res$accdcc1),dm),rep(mean(res$accdcc2),dm)))

round(apply(res$r,2,ess),3)

rownames(par.table) = c('P25','Median','P75','acc.P.')
colnames(par.table) = c('$\\nu$','$b_{1,1}$','$b_{1,2}$','$b_{1,3}$','$b_{1,4}$','$b_{1,5}$','$b_{1,6}$','$b_{1,7}$','$b_{1,8}$','$b_{1,9}$','$b_{1,10}$','$b_{1,11}$','$b_{1,12}$','$b_{1,13}$','$b_{1,14}$',
                        '$b_{2,1}$','$b_{2,2}$','$b_{2,3}$','$b_{2,4}$','$b_{2,5}$','$b_{2,6}$','$b_{2,7}$','$b_{2,8}$','$b_{2,9}$','$b_{2,10}$','$b_{2,11}$','$b_{2,12}$','$b_{2,13}$','$b_{2,14}$')

print(xtable::xtable(par.table,
                     caption = "Posterior descriptive statistics for the parameters
                     of the rank-1 DCC-t model. P25 and P75 are quartiles 1 and 3, 
                     respectively, and acc.P is the acceptance probability.",
                     label = 'table:dcct_FX', digits = 2),
      file='tables_and_figures/pars_dcct_FX.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,rotate.colnames = FALSE,scalebox = 0.5,
      sanitize.text.function = function(x) {x})

# pdf(file='tables_and_figures/parameters_dcct.pdf',width=12,height=8)
par(mfrow=c(6, 5), mar=c(2, 2, 1, 1))
for(i in 1:dim(res$r)[2]) plot(res$r[ind,i],type='l')
# dev.off()

##--------------
## DCC-HEAVY-t
##--------------

load('empirical/temp/results_heavy_t.Rdata')
M   = dim(res$r)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

dim(res$r)

par.table = rbind(apply(res$r,2,p25),
                  apply(res$r,2,median),
                  apply(res$r,2,p75),
                  c(mean(res$accnu),rep(mean(res$accdcc1),dm),rep(mean(res$accdcc2),dm)))

round(apply(res$r,2,ess),3)

rownames(par.table) = c('P25','Median','P75','acc.P.')
colnames(par.table) = c('$\\nu$','$b_{1,1}$','$b_{1,2}$','$b_{1,3}$','$b_{1,4}$','$b_{1,5}$','$b_{1,6}$','$b_{1,7}$','$b_{1,8}$','$b_{1,9}$','$b_{1,10}$','$b_{1,11}$','$b_{1,12}$','$b_{1,13}$','$b_{1,14}$',
                        '$b_{2,1}$','$b_{2,2}$','$b_{2,3}$','$b_{2,4}$','$b_{2,5}$','$b_{2,6}$','$b_{2,7}$','$b_{2,8}$','$b_{2,9}$','$b_{2,10}$','$b_{2,11}$','$b_{2,12}$','$b_{2,13}$','$b_{2,14}$')

print(xtable::xtable(par.table,
                     caption = "Posterior descriptive statistics for the parameters
                     of the rank-1 DCC-t model. P25 and P75 are quartiles 1 and 3, 
                     respectively, and acc.P is the acceptance probability.",
                     label = 'table:dcct_FX', digits = 2),
      file='tables_and_figures/pars_dcct_FX.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,rotate.colnames = FALSE,scalebox = 0.5,
      sanitize.text.function = function(x) {x})

#pdf(file='tables_and_figures/parameters_heavyt.pdf',width=10,height=2)
par(mfrow=c(1, 3), mar=c(2, 2, 1, 1))
for(i in 1:dim(res$r)[2]) plot(res$r[ind,i],type='l')
#dev.off()

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
