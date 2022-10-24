rm(list=ls(all=TRUE))
load('data/3data_EX.Rdata')
library(xtable)

post.sample = 1000

pars_lf = NULL

##------
## Rme
##------

load('temp/results_RMe_EX.Rdata')
M   = length(res$resRMe)
ind = round(seq(1,M,length=post.sample)) #thin every xth
lam = res$resRMe[ind]

pars_lf = rbind(pars_lf,
                c(median(lam),sqrt(var(lam)),mean(res$accRMe)))

#pdf(file='tables_and_figures/RMe_EX.pdf',width=12,height=3)
par(mfrow=c(1,3))
plot(lam,type='l',main='RMe: lambda',ylab='',xlab='')
acf(lam,lwd=2,main='ACF',xlab='')
pacf(lam,lwd=2,main='PACF',xlab='')
#dev.off()

##------
## scalar dcc
##------

load('temp/results_scalar_dcc_EX.Rdata')
M   = dim(res$resdcc)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

a      <- res$resdcc[ind,1]
b      <- res$resdcc[ind,2]

pars_lf = rbind(pars_lf,
                c(median(a),sqrt(var(a)),mean(res$accdcc)),
                c(median(b),sqrt(var(b)),mean(res$accdcc)))

#pdf(file='tables_and_figures/dcc_EX.pdf',width=10,height=5)
par(mfrow=c(2,3))
plot(a,type='l',main='DCC: a',ylab='',xlab='')
acf(a,lwd=2,main='ACF',xlab='')
pacf(a,lwd=2,main='PACF',xlab='')
plot(b,type='l',main='DCC: b',ylab='',xlab='')
acf(b,lwd=2,main='ACF',xlab='')
pacf(b,lwd=2,main='PACF',xlab='')
#dev.off()

##------
## scalar dcc t-copula
##------

load('temp/results_scalar_dcc_t_EX.Rdata')
M   = dim(res$restdcc)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

a      <- res$restdcc[ind,1]
b      <- res$restdcc[ind,2]
nu     <- res$restdcc[ind,3]

#pdf(file='tables_and_figures/dcct_EX.pdf',width=12,height=8)
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
#dev.off()

pars_lf = rbind(pars_lf,
                c(median(a),sqrt(var(a)),mean(res$acctdcc)),
                c(median(b),sqrt(var(b)),mean(res$acctdcc)),
                c(median(nu),sqrt(var(nu)),mean(res$acctdcc)))

df <- data.frame(c('RMe','DCC','','DCC-t','',''),
                 c('$\\lambda$','$a$','$b$','$a$','$b$','$\\eta$'),
                 pars_lf) 

colnames(df) = c('Model','Parameter','Median','St.Dev.',
                      'Acc.prob.')

print(xtable(df,
             caption = 'Parameter estimation results for the low-frequency models:
             RiskMetrics estimated (RMe) and Dynamic Conditional Correlation
             with Gaussian and $t$ copulas (DCC and DCC-t).',
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


#pdf(file='tables_and_figures/xm_1_EX.pdf',height = 5,width = 10)
par(mfrow=c(2,3))
plot(lag,type='l',main='AIW: l2',ylab='',xlab='',ylim=c(min(lag)-1,max(lag)+1))
acf(lag,lwd=2,main='ACF',xlab='')
pacf(lag,lwd=2,main='PACF',xlab='')
plot(nu,type='l',main='AIW: nu',ylab='',xlab='')
acf(nu,lwd=2,main='ACF',xlab='')
pacf(nu,lwd=2,main='PACF',xlab='')
#dev.off()

#pdf(file='tables_and_figures/xm_2_EX.pdf',height = 5,width = 12)
par(mfrow=c(1,3))
for (i in 1:dm) plot(b1[,i],type='l',main=paste('b1_',i,sep=''),
                     ylab='',xlab='')
#dev.off()

#pdf(file='tables_and_figures/xm_3_EX.pdf',height = 5,width = 12)
par(mfrow=c(1,3))
for (i in 1:dm) plot(b2[,i],type='l',main=paste('b2_',i,sep=''),
                     ylab='',xlab='')
#dev.off()


df = data.frame(c('$l_2$','$\\nu$',
                  '$b_{1,1}$','$b_{1,2}$','$b_{1,3}$',
                  '$b_{2,1}$','$b_{2,2}$','$b_{2,3}$'),
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


