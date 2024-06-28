library(MCS)
rm(list=ls(all=TRUE))

load("empirical/temp/data_mctest_3.RData")
load("empirical/temp/oraclesd.Rdata")

#all.model.sd: predicted portfolio stdevs for all models
#all.portsd: "true" portfolio stdevs
#all.prets: ex-post realized portfolio log returns


########################
### Portfolio  GMV #####
########################

sd_true <- t(oracle.portsd[[1]])
sd_est  <- t(all.portsd$gmvr)

out <- dim(sd_est)[1] # out-of-sample
nm <- dim(sd_est)[2] #number of models
nombre <- c("Jores model",
            "Gewekes model",
            "Equally weighted",
            "CAW",
            "DCC-t",
            "DCC-HEAVY-t")

### Testing Stadnard deviation

# Step 1. Loss function: LossVol(realized, evaluated, which = 'SE1','SE2','QLIKE','R2LOG','AE1','AE2')
# the option 'which' selects the type of loss function; 
# SE2 is squared difference of variances


### just a small test to check if the loss function is working
test1 <- LossVol(realized = sd_true[,1], evaluated = sd_est[,1], which = "SE1")
plot(test1, type = "l")
test2 <- LossVol(realized = sd_true[,1], evaluated = sd_est[,2], which = "SE1")
plot(test2, type = "l")
lines(test1, col=2)
###

loss.sd <- matrix(NA, out, nm)
for (i in 1:6){
  loss.sd[,i] <- LossVol(realized = sd_true[,1], evaluated = sd_est[,i], which = "SE1")
}
colnames(loss.sd) <- nombre

# Step 2. MCS procedure at alpha = 5% conf. level

SSM.sd <- MCSprocedure(Loss = loss.sd, alpha = 0.05, B = 10000, statistic = "TR")

########################
### Portfolio  VaR05 ###
########################
sd_true <- t(oracle.portsd[[2]])
sd_est <- t(all.portsd$CVAR05r)


# 1. 
loss.var5 <- matrix(NA, out, nm)
for (i in 1:6){
  loss.var5[,i] <- LossVol(realized = sd_true[,1], evaluated = sd_est[,i], which = "SE2")
}
colnames(loss.var5) <- nombre

# 2.
SSM.var5 <- MCSprocedure(Loss = loss.var5, alpha = 0.05, B = 10000, statistic = "Tmax")

########################
### Portfolio  VaR10 ###
########################
sd_true <- t(oracle.portsd[[3]])
sd_est <- t(all.portsd$CVAR10r)

# 1. 
loss.var10 <- matrix(NA, out, nm)
for (i in 1:6){
  loss.var10[,i] <- LossVol(realized = sd_true[,1], evaluated = sd_est[,i], which = "SE1")
}
colnames(loss.var10) <- nombre

# 2.
SSM.var10 <- MCSprocedure(Loss = loss.var10, alpha = 0.05, B = 10000, statistic = "Tmax")

apply(sd_est,2,mean)

