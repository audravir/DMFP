rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')

nn       = length(date)
end.date = which(zoo::as.yearmon(date)=="jan 2021")[1]-1
#end.date = which(zoo::as.yearmon(date)=="ene 2021")[1]-1
date[end.date]

T0  = end.date  #
T0 # for estimation
K = nn-end.date
K # for oos evaluation

T0+K
nn
# the same v

data = stand[1:T0,]

# the length of the estimation window, 1:T0, which is 

# the second number is MCMC size. The codes have been run for 50k and the results are in
# Rdata files in folder temp. The files are too large for GitHub
# they are stored locally 
# 25k for burn-in+25k 
# thinned every 25th, the resulting posterior samples of 1000k

scalardcc(data,200)


scalartdcc(data,1500)


rmetrics(data,1500)
xm1(Sigma,1500)

