

# Define working directory as path to local project 
wd <- "E:/GitRepos/HerringMSE"

##############################
## Load Libraries and Setup ##
##############################

## Uncomment the line below to install the CA herring version of DLMtool from GitHub
## This will be replace with CRAN version once project is complete
## May require installation of rtools or other software to compile CPP code 
# devtools::install_github("DLMtool/DLMtool", ref="CAHerringMSE")

library(DLMtool)
if(!"mvtnorm" %in% rownames(installed.packages())) install.packages("mvtnorm")
source(paste0(wd, "/Code/makeMPs.r")) # source script to create MPs 

## MSE Setup ##

nsim <- 60 # number of simulations - increase later and re-run script 
proyears <- 50 # number of years to project into future 
set.seed(101) # set seed so random parameter draws are consistent

StockName <- "Stock_Base" # Name of Stock CSV 
FleetName <- "Fleet_Base" # Name of Fleet CSV 
ObsName <- "Obs_Base" # Name of Observation CSV 
Name <- "BaseCase" # Name of object to be saved 
Overwrite <- FALSE # should the saved object be written over if it already exists? 

### Load CSV stock object - 99 or 999 code for parameters to be sampled externally 
Stock <- new("Stock", paste0(wd, "/CSVs/", StockName, ".csv"))

### Load CSV fleet object 
Fleet <- new("Fleet", paste0(wd, "/CSVs/", FleetName, ".csv"))

Obs <- Perfect_Info # use this for now - will source ObsName CSV here later 

# SSB index from stock assessment data 
SSBindex <- c(8169, 21389, 15481, 50482, 29361, 3526, 10550, 9236, 11331, 11682,
              10996, 32845, 58789, 144309, 10601, 10435, 4322, 38409, 55356, 
			  59353, 77002, 57428, 16628, 14405)

qssb <- 0.341 			  
SSBabsolute <- SSBindex/qssb # SSB data - for comparison with simulated data 

####################
# Generate Samples #
####################
source(paste0(wd, "/Code/generatePars.r"))

# script sourced above generates:
#	- nsim samples of correlated von Bertalanffy parameters (from vbPars.csv)
# 	- nsim samples of selectivity-at-age for each year (no variability) 
#	- nsim samples of maturity-at-age (assumed fixed for all simulations and years)
# 	- nsim samples of R0 (assuming CV = 0.22)
# 	- fixed maximum age (currently assumed at 6)
#	- lower and upper bounds of historical relative fishing effort for each historical year (from histEffort.csv)
# 	Above all stored in list 'custompars'

############################################### 
# Create the Operating Model Parameter Object #
###############################################

OM <- new("OM", Stock, Fleet, Obs) # 

##################################### 
# Inspect the OM Parameters 		#
# Comment out to skip this section	#
#####################################

### check that SSB0 is close to K (at least right order of magnitude)
MaxAge <- Stock@maxage
R0 <- Stock@R0 
M <- mean(Stock@M)
ages <- 0:(MaxAge-1)
Linf <- mean(vBpars$Linf)
K <- mean(vBpars$K)
t0 <- mean(vBpars$t0)
Ns <- R0 * exp(-M * ages)
Lens <- Linf * (1-exp(-K * ((ages+1) - t0)))
Wgt <- Stock@a * Lens^Stock@b
matAge <- c(0, 0.36, 0.94, 1, 1, 1) # Table A1.2
sum(Ns * Wgt * matAge) # cf. 84799 estimate of Ksp in Table 9 

### Check Simulated Historical Trajectories ###
setup() # set-up parallel processing
runHist <- runMSE(OM, MPs=NA, nsim=nsim, proyears=proyears, Hist=TRUE, 
    useTestCode=TRUE, custompars=custompars)

plot(runHist$SampPars$dep, 	runHist$SampPars$Depletion, 
  xlab="Sampled Depletion", ylab="Simulated Depletion")

# myOM <- OM
# myOM@AC <- c(0, 0)  
# myOM@Perr <- c(0.8, 1.1)
# runHist <- runMSE(myOM, MPs=NA, nsim=nsim, proyears=proyears, Hist=TRUE, 
    # useTestCode=TRUE, custompars=custompars)

# plot(runHist2$SampPars$dep, runHist2$SampPars$Depletion, 
  # xlab="Sampled Depletion", ylab="Simulated Depletion")
  
## Plot Simulated Trajectories 
matplot(EffYears, t(runHist$SampPars$Find), type="l") # Historical Relative Fishing Effort
lines(EffYears, t(runHist$SampPars$Find)[,sample(1:nsim,1)], lwd=3)  # Example single trajectory

matplot(EffYears, runHist$TSdata$C, type="l") # Historical Simulated Catches 
lines(EffYears, runHist$TSdata$C[,sample(1:nsim,1)], lwd=3)  # Example single trajectory

ndata <- length(SSBabsolute)
fitdata <- c(rep(NA, OM@nyears-ndata), SSBabsolute)
matplot(EffYears, runHist$TSdata$SSB, type="l", ylim=c(0, max(SSBabsolute))) # Simulated SSB
lines(EffYears, runHist$TSdata$SSB[,sample(1:nsim,1)], lwd=3)  # Example single trajectory
lines(EffYears, fitdata, type="l",  lwd=4, col="blue") # SSB data 

matplot(EffYears, runHist$TSdata$VB, type="l") # Simulated vulnerable biomass
lines(EffYears, runHist$TSdata$VB[,sample(1:nsim,1)], lwd=3)  # Example single trajectory
lines(EffYears, fitdata, type="l",  lwd=4, col="blue") # SSB data 
# vulnerable biomass is about half of spawning biomass 		

# NOTE: high auto-correlation makes it difficult to reach depletion in some cases
#	 	may have to decrease AC and/or increase Perr to simulate large fluctuations in SSB

## FMSY and MSY Harvest Rate ##
harvRate <- 1-exp(-runHist$MSYs$FMSYb)
hist(harvRate)
mean(harvRate)

# NOTE: MSY harvest rate (of spawning biomass) is about 35% which is higher 
#		than that reported in stock assessment or MacCall et al earlier review
#		Suggests that under current conditions, fixed harvest rate of ~35% 
#		will result in the highest long-term yields


# Check with older version (slower)
# runHist2 <- runMSE(OM, MPs=NA, nsim=nsim, proyears=proyears, Hist=TRUE, 
    # useTestCode=FALSE, custompars=custompars)
# harvRate2 <- 1-exp(-runHist2$MSYs$FMSYb)
# hist(harvRate2)
# mean(harvRate2)	

	
#####################################
## Run MSE and save object to disk ## 
#####################################

# Create some MPs 
slope10_60 <- makeMinMaxMP(minB=10000, maxB=60000)
slope20_80 <- makeMinMaxMP(minB=20000, maxB=80000)
slope40_100 <- makeMinMaxMP(minB=40000, maxB=100000)
esc15_0.07 <- makeFixEsc(15000, 0.07)
esc20_0.1 <- makeFixEsc(20000, 0.1)
esc5_0.05 <- makeFixEsc(5000, 0.05)
swch0.03_0.07_40_5 <- makeSwitch(0.03, 0.07, 40000, 5000)

fixed0.1 <- makeMinMaxMP(minB=0, maxB=1, PercCap=0.1) # fixed harvest rate at 10% 
fixed0.2 <- makeMinMaxMP(minB=0, maxB=1, PercCap=0.2) # fixed harvest rate at 20% 
fixed0.3 <- makeMinMaxMP(minB=0, maxB=1, PercCap=0.3) # fixed harvest rate at 30% 
fixed0.4 <- makeMinMaxMP(minB=0, maxB=1, PercCap=0.4) # fixed harvest rate at 40% 
fixed0.5 <- makeMinMaxMP(minB=0, maxB=1, PercCap=0.5) # fixed harvest rate at 50% 

# NOTE: 
# Assumption -  HCR acts on absolute spawning stock abundance with TAC applied to 
#				vulnerable biomass (which is about half of SSB)
#	Stock assessment assumes that SSB index is about 1/3 of actual SSB 


Outputs <- avail("DLM_output") # All output controls 
Ind <- NULL 
for (X in seq_along(Outputs)) 
  if (!is.null(attributes(get(Outputs[X]))$case))
    if (attributes(get(Outputs[X]))$case == "CAHerring") Ind <- c(Ind, X)

# Select herring MPs, current effort, FMSY and No Fishing reference methods
MPs <- c(Outputs[Ind], "curE", "FMSYref", "NFref")

MSE <- runMSE(OM, MPs=MPs, nsim=nsim, proyears, interval=1, 
    useTestCode=TRUE, custompars=custompars,
	CheckMPs=FALSE)

	
# Save MSE Object to Disk 	
if (file.exists(paste0(wd, "/SaveObjs/", Name, "_MSE.rds")))  # object already exist
  if(!Overwrite) Name <- paste0(Name, "_", paste0(sample(letters, 5), collapse=''))

saveRDS(MSE, file=paste0(wd, "/SaveObjs/", Name, "_MSE.rds"))


# Initial Explortation of MSE Results 
	
	
NOAA_plot(MSE)	 # under current conditions fixed harvest rate of 30% does best
# MPs with higher harvest rate perform the same because TAC is greater than available catch 

runHist$DLM_data@Abun/runHist$DLM_data@SpAbun # vulnerable biomass about 50% of spawning biomass 

ylim <- c(0, max(MSE@C))
op <- par(mfrow=c(4,4))
for (mm in 1:MSE@nMPs) matplot(t(MSE@C[, mm, ]), type="l", main=MSE@MPs[mm], bty="l", ylim=ylim)
par(op)

ylim <- c(0, max(MSE@C))
op <- par(mfrow=c(4,4))
for (mm in 1:MSE@nMPs) matplot(t(MSE@SSB[, mm, ]), t(MSE@C[, mm, ]), type="l", main=MSE@MPs[mm], bty="l", ylim=ylim)
par(op)





 