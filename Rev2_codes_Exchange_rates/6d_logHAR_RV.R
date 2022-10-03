


##-------------
## log(RV)
##-------------

RVsforc1_m3 = matrix(NA,ncol=dm,nrow=oos)

for(i in 1:dm){
  rvt   = log(RVs[,i])
  rvt1  = c(NA,rvt[-length(rvt)])
  rvt5  = rollmean(rvt1,5,fill = c(NA,NA,NA),align = 'right')
  rvt22 = rollmean(rvt1,22,fill = c(NA,NA,NA),align = 'right')
  dtf   = data.frame(rvt,rvt1,rvt5,rvt22)
  lrvm4 = lm(rvt~rvt1+rvt5+rvt22,data = dtf[1:(nn-oos),])
  RVsforc1_m3[,i]  = exp(predict(lrvm4, newdata = dtf[(nn-oos+1):nn,]))
}

##