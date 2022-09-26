rm(list=ls(all=TRUE))
load('temp/portfolio5.Rdata')

esfun=function(x,p){
  es=mean(x[which(x<quantile(x,p))])
  return(es)
}

resm =matrix(NA,ncol=5,nrow=12)

for(i in 1:5){
  resm[1,i] = median(apply(gvm_ret[i,,],1,quantile,0.05))
  resm[2,i] = median(apply(gvm_ret[i,,],1,quantile,0.1))
  resm[3,i] = median(apply(gvm_ret[i,,],1,esfun,0.05))
  resm[4,i] = median(apply(gvm_ret[i,,],1,esfun,0.1))
  resm[5,i] = median(apply(gmv_to[i,,-1],1,mean))
  resm[6,i] = median(apply(gmv_co[i,,],1,mean))
  resm[7,i] = median(apply(gmv_sp[i,,],1,mean))
  
  resm[8,i] = median(apply(gvm_ret[i,,],1,sum))#annualized
  resm[9,i] = median(apply(gvm_ret[i,,],1,sd)*sqrt(252))
  
  resm[10,i] = median(apply(gvm_ret[i,,-1]-0.01/252*gmv_to[i,,-1],1,sum))
  
  resm[11,i] = median(apply(gvm_ret[i,,-1]-0.01/252*gmv_to[i,,-1],1,sum)/
                       (apply(gvm_ret[i,,-1],1,sd)*sqrt(252)))

  resm[12,i] = median(100*sqrt(252)*(apply(gvm_ret[4,,],1,sd)-apply(gvm_ret[i,,],1,sd))/
                        apply(gvm_ret[i,,],1,sd))
}

resm = resm[,c(2,1,3,4,5)]
colnames(resm) = c( 'Geweke','Jore1','Equal','AIW','DCC-t')
rownames(resm) = c('VaR5%','VaR10%','ES5%','ES10%','TO','CO','SP',
                   'return','sd','adj.return(1%)',
                   'adj.Sharpe','G/L')

round(resm,3)


print(xtable(resm,
             caption = "GMV portfolio results based on 1-step-ahead predicitons 
             for  2009/01/02-2009/12/31 out-of-sample period ($K=252$ observations).
             The table reports the posterior medians of various Global Minimum
             Variance portfolio metrics for the pooled models: 
             Geweke's, Jore's and equally weighted, 
             as well as two best individual models, Additive Inverse Wishart (AIW) and 
             Dynamic Conditional Correlation with $t$ copula (DCC-t).",
             label = 'table:gmvfull5', digits = 3),
      file='tables_and_figures/gmvfull5.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,
      hline.after = getOption("xtable.hline.after", c(-1,0,7,nrow(resm))))
