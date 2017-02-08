
# Define working directory as path to local project 
wd <- "E:/GitRepos/HerringMSE"

## Load MSE Object ##
Name <- "BaseCase"

MSE <- readRDS(file=paste0(wd, "/SaveObjs/", Name, "_MSE.rds"))

Converge(MSE) # check MSE for convergence

Pplot(MSE)
# Info stored in MSE 
slotNames(MSE)

MSE@MPs 
MSE@nMPs 
MSE@nsim 
dim(MSE@C) # Catch - nsim, nMPs, proyears 
dim(MSE@SSB) # SSB nsim, nMPs, proyears 
dim(MSE@B_BMSY) # B/BMSY nsim, nMPs, proyears 

MSE@OM$RefY # reference yield for each simulation


## Unfished ## 
ind <- which(MSE@MPs == "NFref")
SSB0 <- MSE@OM$SSB0 # 
hist(SSB0, col="gray", breaks=40)
matplot(t(MSE@SSB[,ind,]), type="l") # Unfished SSB projections 


ylim <- c(0, max(MSE@C))
op <- par(mfrow=c(4,4))
for (mm in 1:MSE@nMPs) matplot(t(MSE@C[, mm, ]), type="l", main=MSE@MPs[mm], bty="l", ylim=ylim)
par(op)

ylim <- c(0, max(MSE@C))
op <- par(mfrow=c(4,4))
for (mm in 1:MSE@nMPs) matplot(t(MSE@SSB[, mm, ]), t(MSE@C[, mm, ]), type="l", main=MSE@MPs[mm], bty="l", ylim=ylim)
par(op)

