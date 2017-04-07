
setwd(paste0(wd, "/SaveObjs/"))
MSE <- readRDS("BaseCase_MSE.rds")

Converge(MSE) # check convergence 

# Load libraries 
library(dplyr)

## Performance Metrics and Plots ##

# Limit Reference Point
fstyr <- 1  # first reference year 
lstyr <- 20 # last reference year
Blim <- 0.4 # B/BMSY limit reference
Plim <- 0.1 # maximum risk threshold

lstyr <- min(lstyr, proyears)
nyrs <- length(fstyr:lstyr)
proyears <- MSE@proyears
nMPs <- MSE@nMPs
nsim <- MSE@nsim
rnd <- 2 # round to 2 digits

probs <- round(apply(MSE@B_BMSY[,, fstyr:lstyr] <= Blim, 2, sum)/(nsim*nyrs), rnd)

LimitPM <- data.frame(MP=MSE@MPs, Prob=probs, Yrs=paste0(fstyr, ":", lstyr), Pass=probs <=Plim, Ref=paste0("SSB < ", Blim, "SSB_MSY"))

# MPs that fail to meet the Limit Reference Point 
LimitPM %>% filter(Pass == FALSE) %>% select(MP, Prob)

# MPs that pass the Limit Reference Point 
LimitPM %>% filter(Pass == TRUE) %>% select(MP, Prob)
passMPs <- as.character((LimitPM %>% filter(Pass == TRUE) %>% select(MP))$MP)


# Subset MSE object for only MPs that pass the Limit Reference Point # 
SubMSE <- Sub(MSE, MPs=passMPs)
nMPs <- SubMSE@nMPs
MPs <- SubMSE@MPs



## Biological Performance Metrics ##

# Biomass above 40% SSBMSY over 50 years 
B40_50 <- round(apply(SubMSE@B_BMSY[,, 1:50] >= 0.4, 2, sum)/(nsim*50), rnd)

# Short term target reference point - probability of being at or above 80% SSBMSY within 20 years.
# Biomass above 80% SSBMSY over 20 years 
B80_20 <- round(apply(SubMSE@B_BMSY[,, 1:20] >= 0.8, 2, sum)/(nsim*20), rnd)

# Spawning Stock Biomass relative to 80% SSBMSY 
# IN LAST PROJECTION YEAR? #
# NOT SURE WHAT RELATIVE TO 80% SSBMSY MEANS 
apply(SubMSE@B_BMSY[,,proyears]/0.8, 2, mean) # relative to 80% MSY - but hard to interpret?
meanFinalBio <- apply(SubMSE@B_BMSY[,,proyears], 2, mean) # relative to MSY 


## Age Metrics - age structure currently not returned in MSE object 
## Have added this to Herring development version of DLMtool 
# Mean age of the herring population - protecting the age structure of herring came up as one of the management objectives for this fishery
# MEAN AGE IN LAST PROJECTION YEAR? 

ages <- 1:dim(SubMSE@PopAS)[3]
nfref <- which(SubMSE@MPs == "NFref")
if (length(nfref) > 0) {
  SubMSE@CatchAS[, nfref, ,] <- 0  # make zero catch for no fishing reference
  SubMSE@C[, nfref,] <- 0  # needed because no-fishing method actually sets catch very low
} 

lstAS <- SubMSE@PopAS[,,, proyears]
ind <- as.matrix(expand.grid(1:nsim, 1:nMPs, ages))
temp <- array(NA, dim=dim(lstAS))
temp[ind] <- (lstAS[ind] * ind[,3])
meanAge <- apply(temp, c(1,2), sum)/apply(lstAS, c(1,2), sum) # mean by simulation and MP 
meanAgePop <- apply(meanAge, 2, mean) # average across simulations 


# Mean age of the herring catch
# MEAN AGE IN LAST PROJECTION YEAR? 
lstAS2 <- SubMSE@CatchAS[,,, proyears]
ind2 <- as.matrix(expand.grid(1:nsim, 1:nMPs, ages))
temp2 <- array(NA, dim=dim(lstAS2))
temp2[ind2] <- (lstAS2[ind2] * ind2[,3])
meanAgeC <- apply(temp2, c(1,2), sum)/apply(lstAS2, c(1,2), sum) # mean by simulation and MP 
meanAgeCatch <- apply(meanAgeC, 2, mean, na.rm=TRUE) # average across simulations 
# Average age of catch IF there is catch (i.e TAC != 0)

# Percentage of age 5+ herring in the catch
# IN LAST PROJECTION YEAR? 
lstAS <- SubMSE@CatchAS[,,, proyears]
P5up <- function(x) {
  probs <- x/sum(x)
  sum(probs[5:length(probs)])
}  

temp <- apply(lstAS, c(1,2), P5up) # proportion of catch 5+ by simulation and MP 
mean5up <- apply(temp, 2, mean, na.rm=TRUE) # average across simulations 


# Biological Performance Metrics Table 
BPM <- data.frame(MP=MPs, B40_50=B40_50, B80_20=B80_20, meanB_BMSY=meanFinalBio, meanAgePop=meanAgePop,
  meanAgeCatch=meanAgeCatch, mean5up=mean5up)
BPM 


## Forage Fish Reference Points ## 
# - will add - just variants of above 


## Yield Performance Metrics ##
# WE USUALLY PRESENT RELATIVE TO REFERENCE YIELD BUT THESE ARE ABSOLUTE 
# NEED TO CONFIRM THE UNITS OF CATCH, TAC, WEIGHT-PARAMETERS, ETC 
# ALSO ASSUMING TAC IS CAUGHT - NO-UNDER CATCH 

# Short-term (20 yrs) cumulative catch 
sty <- apply(SubMSE@C[,,1:20], c(1,2), sum)
# Long-term (50 yrs) cumulative catch 
lty <- apply(SubMSE@C, c(1,2), sum)
# Average annual catch (last 10 years)
aay <- apply(SubMSE@C[,,41:50], 2, mean)


## Catch Variability ## 
# average CV of catch in all years
meanCV <- apply(apply(SubMSE@C, c(1,2), sd)/apply(SubMSE@C, c(1,2), mean), 2, mean)

# Proportion of years with non-zero quota
Pnonzero <- apply(SubMSE@C > 0, 2, sum) / (nsim * proyears) 

# Average number of consecutive years of non - zero quota.
consecNonZero <- function(x) { # calculates mean consecutive years with non-zero catch 
  x[x != 0] <- 9999
  temp <- rle(x)
  if (all(temp$values == 0)) return(0)
  mean(temp$lengths[temp$values == 9999])
}

AvgConsecC <- apply(apply(SubMSE@C, c(1,2), consecNonZero), 2, mean) # mean consecutive years with non-zero catch 


# Percentage of years at the maximum exploitation rate (10%).
# NOTE: HCRs are based on SSB but exploitation rate here is calculated from vulnerable biomass
# which is about half SSB. So exploitation rate is HIGHER than that set in the HCR  
exrate <- SubMSE@C/SubMSE@VB # change to SubMSE@SSB to calculate exploitation rate given SSB 
PyrEx <- apply(exrate > 0.1, 2, sum)/(proyears * nsim) # proportion of years exploitation rate > 10%

# Average number of consecutive years of maximum exploitation rate (10%)
consecMaxEx <- function(x, max=0.1) { # calculates mean consecutive years with exploitation rate above maximum ) 
  x[x > max] <- 9999
  x[x <= max] <- 0
  temp <- rle(x)
  if (all(temp$values == 0)) return(0)
  mean(temp$lengths[temp$values == 9999])  
}

AvgConseEx <- apply(apply(exrate, c(1,2), consecMaxEx), 2, mean) # mean consecutive years with non-zero catch 

#################
##### PLOTS #####
#################

# Barplot - Biological Performance Metrics  
BPM

refMP <- which(MPs == "fixed0.05")

par(mfcol=c(2,6), oma=c(10,2,1,1), mar=c(2,4,2,2))
cols <- c("gray", "black")
# Probablity biomass is above 40% BMSY over 50 years 
ylab <- "P B > 0.4BMSY (50 years)"
barplot(BPM$B40_50, ylim=c(0,1), las=2, col=cols, ylab=ylab)
barplot(BPM$B40_50 - BPM$B40_50[refMP], names=MPs, las=2, col=cols, ylab=ylab)

# Probablity biomass is above 80% BMSY over 20 years 
ylab <- "P B > 0.8BMSY (20 years)"
barplot(BPM$B80_20, ylim=c(0,1), las=2, col=cols, ylab=ylab)
barplot(BPM$B80_20 - BPM$B80_20[refMP], names=MPs, las=2, col=cols, ylab=ylab)

# Mean B/BMSY in last year - MAY NEED TO CHANGE THIS TO 0.8BMSY
ylab <- "average B/BMSY in last year"
barplot(BPM$meanB_BMSY, ylim=c(0,2), las=2, col=cols, ylab=ylab)
barplot(BPM$meanB_BMSY - BPM$meanB_BMSY[refMP], names=MPs, las=2, col=cols, ylab=ylab)

# Mean age pop 
ylab <- "average age in population in last year"
barplot(BPM$meanAgePop, ylim=c(0,6), las=2, col=cols, ylab=ylab)
barplot(BPM$meanAgePop - BPM$meanAgePop[refMP], names=MPs, las=2, col=cols, ylab=ylab)

# Mean age catch  
ylab <- "average age in catch in last year"
barplot(BPM$meanAgeCatch, ylim=c(0,6), las=2, col=cols, ylab=ylab)
barplot(BPM$meanAgeCatch - BPM$meanAgeCatch[refMP], names=MPs, las=2, col=cols, ylab=ylab)

# Proportion catch 5+ 
ylab <- "proportion 5+ in catch in last year"
barplot(BPM$mean5up, ylim=c(0,1), las=2, col=cols, ylab=ylab)
barplot(BPM$mean5up - BPM$mean5up[refMP], names=MPs, las=2, col=cols, ylab=ylab)


# Boxplots - can customize 
boxplot(MSE)


# Kite plots 






