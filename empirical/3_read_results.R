rm(list=ls(all=TRUE))
load('data/FXdata.Rdata')
library(xtable)
post.sample = 1000

pars_lf = NULL
dm      = dim(rets)[2]

##------
## vector dcc
##------

load('empirical/temp/results_vectordcc_tcop.Rdata')
M   = dim(res$resdcc)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

mean(res$accdcc)
mean(res$accnu)
dim(res$resdcc)

a      <- res$resdcc[ind,2:(dm+1)]
b      <- res$resdcc[ind,(dm+2):(dm*2+1)]
nu     <- res$resdcc[ind,1]

par(mfrow=c(1,1)) 
plot(nu,type='l')

par(mfrow=c(3,5)) 
for(i in 1:dm) plot(a[,i],type='l')

par(mfrow=c(3,5)) 
for(i in 1:dm) plot(b[,i],type='l')



# pars_lf = rbind(pars_lf,
#                 c(median(a),sqrt(var(a)),mean(res$accdcc)),
#                 c(median(b),sqrt(var(b)),mean(res$accdcc)))
# 
# pdf(file='tables_and_figures/dcc_EX.pdf',width=10,height=5)
# par(mfrow=c(2,3))
# plot(a,type='l',main='DCC: a',ylab='',xlab='')
# acf(a,lwd=2,main='ACF',xlab='')
# pacf(a,lwd=2,main='PACF',xlab='')
# plot(b,type='l',main='DCC: b',ylab='',xlab='')
# acf(b,lwd=2,main='ACF',xlab='')
# pacf(b,lwd=2,main='PACF',xlab='')
# dev.off()

##------
## scalar dcc t-copula
##------

load('temp/results_scalar_dcc_t_EX.Rdata')
M   = dim(res$restdcc)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

a      <- res$restdcc[ind,1]
b      <- res$restdcc[ind,2]
nu     <- res$restdcc[ind,3]

pdf(file='tables_and_figures/dcct_EX.pdf',width=12,height=8)
par(mfrow=c(3,3))
plot(a,type='l',main='DCC-t: a',ylab='',xlab='')
acf(a,lwd=2,main='ACF',xlab='')
pacf(a,lwd=2,main='PACF',xlab='')
plot(b,type='l',main='DCC-t: b',ylab='',xlab='')
acf(b,lwd=2,main='ACF',xlab='')
pacf(b,lwd=2,main='PACF',xlab='')
plot(nu,type='l',main='DCC-t: eta',ylab='',xlab='')
acf(nu,lwd=2,main='ACF',xlab='')
pacf(nu,lwd=2,main='PACF',xlab='')
dev.off()

pars_lf = rbind(pars_lf,
                c(median(a),sqrt(var(a)),mean(res$acctdcc)),
                c(median(b),sqrt(var(b)),mean(res$acctdcc)),
                c(median(nu),sqrt(var(nu)),mean(res$acctdcc)))

##------
## DCC-HEAVY-t
##------

load('temp/results_dcch_t_EX.Rdata')
M   = length(res$acctdcch)
ind = round(seq(1,M,length=post.sample)) #thin every xth

a      <- res$restdcch[ind,1]
b      <- res$restdcch[ind,2]
nu     <- res$restdcch[ind,3]

pdf(file='tables_and_figures/dcct_heavy_EX.pdf',width=12,height=8)
par(mfrow=c(3,3))
plot(a,type='l',main='DCC-HEAVY-t: a',ylab='',xlab='')
acf(a,lwd=2,main='ACF',xlab='')
pacf(a,lwd=2,main='PACF',xlab='')
plot(b,type='l',main='DCC-HEAVY-t: b',ylab='',xlab='')
acf(b,lwd=2,main='ACF',xlab='')
pacf(b,lwd=2,main='PACF',xlab='')
plot(nu,type='l',main='DCC-t: eta',ylab='',xlab='')
acf(nu,lwd=2,main='ACF',xlab='')
pacf(nu,lwd=2,main='PACF',xlab='')
dev.off()


pars_lf = rbind(pars_lf,
                c(median(a),sqrt(var(a)),mean(res$acctdcc)),
                c(median(b),sqrt(var(b)),mean(res$acctdcc)),
                c(median(nu),sqrt(var(nu)),mean(res$acctdcc)))

df <- data.frame(c('RMe','DCC','','DCC-t','','','DCC-HEAVY-t','',''),
                 c('$\\lambda$','$a$','$b$','$a$','$b$','$\\eta$',
                   '$a$','$b$','$\\eta$'),
                 pars_lf) 

colnames(df) = c('Model','Parameter','Median','St.Dev.',
                 'Acc.prob.')

print(xtable(df,
             caption = 'Parameter estimation results for the low-frequency models:
             RiskMetrics estimated (RMe), Dynamic Conditional Correlation
             with Gaussian and $t$ copulas (DCC and DCC-t) and DCC-HEAVY with $t$ copula.',
             label = 'table:all_pars_lf_EX', digits = 4),
      file='tables_and_figures/all_pars_lf_EX.tex',
      include.rownames = FALSE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,size="\\footnotesize",
      sanitize.text.function=function(x){x})
df

##------
## xm1
##------

load('temp/results_xm1_EX.Rdata')
M   = dim(res$resc)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth
dm  = dim(stand)[2]

lag = res$resc[ind,1]
nu  = res$resc[ind,2]
b1  = res$resc[ind,3:(dm+2)]
b2  = res$resc[ind,(dm+3):(2*dm+2)]

pars_hf = NULL
pars_hf = rbind(pars_hf,
                c(median(lag),sqrt(var(lag)),mean(res$accl)),
                c(median(nu),sqrt(var(nu)),mean(res$accnu)))

allbs  = matrix(0,ncol=3,nrow=dm*2)
for(i in 1:dm){
  allbs[i,1] = median(b1[,i])
  allbs[i,2] = sqrt(var(b1[,i]))
  allbs[i,3] = mean(res$accB)
  allbs[i+dm,1] = median(b2[,i])
  allbs[i+dm,2] = sqrt(var(b2[,i]))
  allbs[i+dm,3] = mean(res$accB)
}

pars_hf = rbind(pars_hf,allbs)


pdf(file='tables_and_figures/xm_1_EX.pdf',height = 5,width = 10)
par(mfrow=c(2,3))
plot(lag,type='l',main='AIW: l2',ylab='',xlab='',ylim=c(min(lag)-1,max(lag)+1))
acf(lag,lwd=2,main='ACF',xlab='')
pacf(lag,lwd=2,main='PACF',xlab='')
plot(nu,type='l',main='AIW: nu',ylab='',xlab='')
acf(nu,lwd=2,main='ACF',xlab='')
pacf(nu,lwd=2,main='PACF',xlab='')
dev.off()

pdf(file='tables_and_figures/xm_2_EX.pdf',height = 5,width = 12)
par(mfrow=c(2,ceiling(dm/2)))
for (i in 1:dm) plot(b1[,i],type='l',main=paste('b1_',i,sep=''),
                     ylab='',xlab='')
dev.off()

pdf(file='tables_and_figures/xm_3_EX.pdf',height = 5,width = 12)
par(mfrow=c(2,ceiling(dm/2)))
for (i in 1:dm) plot(b2[,i],type='l',main=paste('b2_',i,sep=''),
                     ylab='',xlab='')
dev.off()

bnames = rep(NA,2*dm)

for(i in 1:dm){
  bnames[i]=paste("$b_{1,",i,"}$",sep="")
  bnames[i+5]=paste("$b_{2,",i,"}$",sep="")
}


df = data.frame(c('$l_2$','$\\nu$',bnames),
                pars_hf)

colnames(df) = c('Parameter','Median','St.Dev.',
                 'Acc.prob.')

print(xtable(df,
             caption = 'Parameter estimation results for the high-frequency AIW model.',
             label = 'table:all_pars_hf_EX', digits = 4),
      file='tables_and_figures/all_pars_hf_EX.tex',
      include.rownames = FALSE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,size="\\footnotesize",
      sanitize.text.function=function(x){x})

df


