
N     = 10000
mylik = function(A,B,lam) res = lam*A+(1-lam)*B
lam   = runif(N)
A     = dnorm(seq(0,3,length=6))
B     = A[1]
trans = list()

pdf(file='tables_and_figures/DNsim.pdf',width=10,height=6)
par(mfrow=c(2,3))
for(i in 1:6){
  ws    = mylik(A[i],B,lam)
  ind   = sample(N,N,replace = TRUE,prob = ws/sum(ws))
  hist(lam[ind],main=paste('A=',round(A[i],3),
                           ', B=',round(B,3),
                           ', Me(DN)=',round(median(lam[ind]),2),
                           ', J1=',round(A[i]/(A[i]+B),2),sep=''),
       xlab='',ylab='')
  abline(v=A[i]/(A[i]+B),lwd=5)
  abline(v=mean(lam[ind]),lwd=4,col=2)
  trans = c(trans,list(qnorm(lam[ind])))
}
dev.off()

