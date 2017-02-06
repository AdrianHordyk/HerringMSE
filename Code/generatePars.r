
###############################################################################
##																			  ##
## This code produces samples of parameters for the MSE and estimates fishing ##
## effort rates for historical years to fit the model to historical SSB data  ##
##                           												  ##
## The sampled parameters are saved out to an R object  					  ##
##																			  ##
## This code must be re-run whenever any of the following are changed:		  ##
## 		- number of simulations (nsim)										  ##	
##		- number of projection years (proyears)								  ##
##		- ANY of the parameters in the Stock or Fleet objects				  ##
##		- SSB data is changed (model is optimized to fit this)				  ##
##      - parameters in vbPars.csv 											  ##
##                           												  ##
################################################################################

# A. Hordyk 
# Feb 2017 

################################################################################
## Assumptions and Notes:                                                     ##
##                                                                            ##
## - R0 is normally distributed with CV of 0.22 (Table 9)					  ##
##                                                                            ##
## - SSB index is assumed to be relative index of abundance (model 19)		  ##
##   qssb = 0.341 (Table 9) 												  ##
##                                                                            ##
## - maturity@age assumed to be equivalent to current selectivity-at-age      ##
##   This was neccessary because otherwise SSB was about twice as high as     ##
##   vulnerable biomass in the simulations. Model has been optimized to the   ##
##   SSB data from assessment. The assumption is that this is the amount of   ##
##   biomass that can be exploited by the fishery.  This differs from Stock   ## 
##   Assessment Figure 24. Need to look at this in more detail                ##
##                                                                            ##
## - need to confirm that units of SSB data is correct - think that it is     ##
##                                                                            ##
## Issues:                                                                    ##
## - SSB0 from simple equilibrium model is considerably lower than SSB0       ##
##   estimated in stock assessment. May be caused by:                         ##
##		- assumptions of maturity@age                                         ##
## 		- growth parameters                                                   ##
##		- natural mortality                                                   ##
## 		- maximum age                                                         ##
##      - assumptions of R0 or SSB0                                           ##
##                                                                            ##
################################################################################

##############################
## Load Libraries and Setup ##
##############################

## Uncomment the line below to install the CA herring version of DLMtool from GitHub
## This will be replace with CRAN version once project is complete
## May require installation of rtools or other software to compile CPP code 
# devtools::install_github("DLMtool/DLMtool", ref="CAHerringMSE")

library(DLMtool)
setup() # set-up parallel processing 
if(!"mvtnorm" %in% rownames(installed.packages())) install.packages("mvtnorm")

# Define working directory as path to local project 
wd <- "E:/GitRepos/HerringMSE"

# number of simulations 
nsim <- 36 # increase later and re-run script 
proyears <- 50 # number of years to project into future 
set.seed(101) # set seed so random parameter draws are consistent

Name <- "BaseCase" # Name of object to be saved 
StockName <- "Stock_Base"
FleetName <- "Fleet_Base"

# SSB index from stock assessment data 
SSBindex <- c(8169, 21389, 15481, 50482, 29361, 3526, 10550, 9236, 11331, 11682,
              10996, 32845, 58789, 144309, 10601, 10435, 4322, 38409, 55356, 
			  59353, 77002, 57428, 16628, 14405)

qssb <- 0.341 			  
SSBabsolute <- SSBindex/qssb
		  
#####################################
## Load MSE Objects from CSV Files ##
#####################################

### Load CSV stock object - 99 or 999 code for parameters to be sampled externally 
Stock <- new("Stock", paste0(wd, "/CSVs/", StockName, ".csv"))

### Generate correlated samples of von Bert parameters ###
parsIn <- read.csv(paste0(wd, "/CSVs/vbPars.csv"))
means <- parsIn$Mean[1:3]
names(means) <- c("Linf", "K", "t0")
vcovar <- as.matrix(parsIn[1:3, 4:6])
colnames(vcovar) <- NULL; rownames(vcovar) <- NULL 
vBpars <- as.data.frame(mvtnorm::rmvnorm(nsim, mean=means, sigma=vcovar))

##****************************##

################
# PROBLEM HERE #
################

### check that SSB0 is close to K (at least right order of magnitude)
MaxAge <- Stock@maxage
R0 <- Stock@R0 
M <- mean(Stock@M)
ages <- 0:(MaxAge-1)
Linf <- mean(vBpars$Linf)
K <- mean(vBpars$K)
t0 <- mean(vBpars$t0)
Ns <- R0 * exp(-M * ages)
Lens <- Linf * (1-exp(-K * (ages - t0)))
Wgt <- Stock@a * Lens^Stock@b
matAge <- c(0.00459511, 0.0122648, 0.192453, 0.647015, 1, 1) # assume same as sel2 below
matAge <- c(0, 0.36, 0.94, 1, 1, 1) # Table A1.2
sum(Ns * Wgt * matAge) # cf. 84799 estimate of Ksp in Table 9 
sum(Ns * Wgt * matAge)/84799
# - right order of magnitude but considerably lower 

##****************************##

### Load CSV fleet object 
Fleet <- new("Fleet", paste0(wd, "/CSVs/", FleetName, ".csv"))

### Sample historical trend in historical fishing effort
#### Assuming fishery started in 1972
#### Using trends in estimated F from 1992 - 2015 using stock assessment 
histEffs <- read.csv(paste0(wd, "/CSVs/histEffort.csv"))
EffYears <- histEffs[,1]
Esd <- runif(nsim, Fleet@Fsd[1], Fleet@Fsd[2])

# fishing effort can be sampled anywhere between Lower and Upper bounds  
# mean trend sampled and then Esd added as log-normal error
# this is starting range of fishing effort 
# - effort is optimized together with R0 to fit catch time-series below  
Find <- getEffhist(Esd=Esd, EffYears=EffYears, EffLower=histEffs[,4], 
  EffUpper=histEffs[,3], nyears=length(EffYears))[[1]]
EffLower <- histEffs[,4]  
EffUpper <- histEffs[,3]
matplot(EffYears, range01(histEffs[3:4]), type="l", lwd=4) # fishing effort scaled between 0 and 1 
matplot(EffYears, t(Find), add=TRUE, type="l")

### Sample Selectivity at Age - assuming no variability here 
#### easy to add variability later if required 
V <- array(NA, dim=c(nsim, Stock@maxage, Fleet@nyears+proyears))
## Assuming fixed selectivity at age (no variation between simulations)
##  
## 1972 - 1997 - sel1 
## 1998 - future - sel2 
sel1 <- c(0.000740331, 0.0626284, 0.280703, 0.631475, 1, 1)
sel2 <- c(0.00459511, 0.0122648, 0.192453, 0.647015, 1, 1)

years <- 1972:(2015+proyears)
ind <- as.matrix(expand.grid(1:nsim, 1:Stock@maxage, 1:sum(years<1998)))
V[ind] <- sel1[ind[,2]]
ind2 <- as.matrix(expand.grid(1:nsim, 1:Stock@maxage, (sum(years<1998)+1):(sum(years<1998)+sum(years>=1998))))
V[ind2] <- sel2[ind2[,2]]
matplot(V[1,,], type="l") # plot the selectivity-at-age curves
sum(is.na(V)) # check that there are no NAs left 

## Maturity - make maturity-at-age equal to selectivity-at-age ###
## If age@maturity parameters from SA are used SSB is about
## twice vulnerable biomass
Mat_age <- matrix(sel2, nrow=nsim, ncol=Stock@maxage, byrow=TRUE)

## Create Operating Model Object to use for optimization
optOM <- new("OM", Stock, Fleet, Perfect_Info)

# List of externally sampled parameters to supply to MSE 
custompars <- list(Linf=vBpars$Linf, K=vBpars$K, t0=vBpars$t0, V=V,
 maxage=Stock@maxage, Mat_age=Mat_age, EffLower=EffLower, EffUpper=EffUpper)

##############
# Test model # Unfished Historical Simulations 
##############
# NOTE: This hasn't been optimized to fit SSB index yet 
useOM <- optOM 

useOM@Fsd <- c(0,0)
testPars <- custompars
testPars$EffLower <- rep(0, useOM@nyears)
testPars$EffUpper <- rep(0, useOM@nyears) # make fishing effort zero for most years 
testPars$EffUpper[40:44] <- 1 

runHist <- runMSE(useOM, MPs=NA, nsim=nsim, proyears=proyears, Hist=TRUE, 
    useTestCode=TRUE, ntrials=50, fracD=0.05, custompars=testPars)

hist(runHist$DLM_data@OM$SSB0) # distribution of SSB0 generated by model 

## Unfished except for last 5 years ## 
Years <- 1972:2015
par(mfrow=c(2,2), bty="l")
matplot(Years, runHist$TSdata$Bio, type="l", ylab="Biomass (tonnes)")
matplot(Years, runHist$TSdata$SSB, type="l", ylab="Spawning Biomass (tonnes)")
matplot(Years, runHist$TSdata$VB, type="l", ylab="Vulnerable Biomass (tonnes)")
matplot(Years, runHist$TSdata$C, type="l", ylab="Catch (tonnes)")


	
	
	
	
	
	
	
### Optimize when Stock and Fleet Objects are Finalised ###	
	
###########################################################
## Optimize Fishing Effort to fit to historical SSB data ##
###########################################################

# Fit model to SSBabsolute by estimating Find for each year 
optR0_Eff <- function(Pars, SSBabsolute, OM, nsim, proyears, seed, custompars=NULL) {
  custompars$R0 <- exp(Pars[1]) * 1E9  
  custompars$Find <- matrix(exp(Pars[2:(OM@nyears+1)]), nrow=1, ncol=OM@nyears)
  set.seed(seed)
  runHist <- runMSE(OM, MPs, nsim, proyears, custompars=custompars, CheckMPs, 
    Hist=TRUE, useTestCode=TRUE, ntrials=50, fracD=0.05)
    
  ndata <- length(SSBabsolute)  
  
  datout <- runHist$TSdata$SSB[(OM@nyears-ndata+1):OM@nyears,]
  
  Vadv <- 0.608
  sigmaI <- 0.05 # stock assessment pg 63 
  like1 <- (log(datout) - log(SSBabsolute))/(sqrt(sigmaI^2 + Vadv^2))
  obj <- 0.5 * sum(like1^2 + log(2*pi*sigmaI^2+Vadv^2)) # sum((log(datout) - log(SSBabsolute))^2) 
  
  Effs <- custompars$Find[1,]
  if (max(abs(diff(Effs))) <= 0.5) {
    penalty <- 0
  } else {
	# add penalty if inter-annual change in fishing effort greater than 0.5
    penalty <- max(abs(diff(Effs))) * 100 * obj#
	# aim is to stop enormous spikes in effort in single year 
  }
   
  pen2 <- -log(dnorm(custompars$R0, OM@R0, OM@R0*0.22))
  
  obj <- obj + penalty + pen2 
 
  # ndata <- length(SSBabsolute)
  # fitdata <- c(rep(NA, OM@nyears-ndata), SSBabsolute) 	
  # plot(runHist$TSdata$VB[,1], type="l", ylab="Vulnerable Biomass (tonnes)")
  # lines(fitdat, type="l",  lwd=4) 
  return(obj)
}

OptFun <- function(s, seeds, useOM, SSBabsolute, custompars, Find,  proyears=proyears) {
  library(DLMtool)
  ## Make function to generalize this ##
  custompars1sim <- custompars
  custompars1sim$Linf <- custompars$Linf[s]
  custompars1sim$K <- custompars$K[s]
  custompars1sim$t0 <- custompars$t0[s]
  custompars1sim$V <- custompars$V[s,,,drop=FALSE]
  # custompars1sim$R0 <- custompars$R0[s]
  custompars1sim$Mat_age <- custompars$Mat_age[s,,drop=FALSE]
   
  seed <- seeds[s]
  
  Effs <- Find[s, ]
  starts <- log(c(useOM@R0/1E9, Effs)) 
 
  test <- optim(starts, optR0_Eff, SSBabsolute=SSBabsolute, 
    OM=useOM, nsim=1, proyears=proyears, seed=seed, custompars=custompars1sim, 
	control=list(maxit=2000))

  # estR0 <- exp(test$par[1])*1E9
  # estR0/useOM@R0   
  # plot(exp(test$par[2:45]), type="l") 
  custompars1sim$R0 <- exp(test$par)[1] * 1E9
  custompars1sim$Find <- matrix(exp(test$par)[2:(useOM@nyears+1)], nrow=1)

  set.seed(seed)
  runHist <- runMSE(useOM, MPs=NA, nsim=1, proyears=proyears, Hist=TRUE, 
    useTestCode=TRUE, ntrials=50, fracD=0.05, custompars=custompars1sim)
  
  # plot(runHist$TSdata$VB[,1], type="l", ylab="Vulnerable Biomass (tonnes)")
  # lines(fitdat, type="l",  lwd=4) 
 
 runHist$SampPars
}

joinSFOb <- function(dat) {
  nsim <- dim(dat)[2]
  npars <- dim(dat)[1]
  datout <- dat[,1] 

  for (X in 1:npars) {
    if (class(dat[X,1][[1]]) == "matrix") {
      len <- dim(dat[X,1][[1]])[2]
	  datout[[X]] <- matrix(unlist(dat[X,]), ncol = len, byrow = TRUE)
	} 
    if (class(dat[X,1][[1]]) == "numeric") {
       datout[[X]] <- unlist(lapply(dat[X,], "[[", 1))
	}
    if (class(dat[X,1][[1]]) == "array") {
	   dim1 <- dim(dat[X,1][[1]])
	   ind <- which(dim1 == 1)
	   dim2 <- dim1 
	   dim2[ind] <- nsim
	   out <- array(NA, dim=dim2)
	   for (x in 1:nsim) {
	     out[x, ,] <- dat[X,][[x]]
	   }
       datout[[X]] <- out
	}		
  }
 datout
}


 
# test <- OptFun(s=1, seeds, optOM, custompars, Find,  proyears=proyears)
# plot(test$Find[1,], type="l")

seeds <- 1:nsim
snowfall::sfExport(list = c("seeds", "custompars", "Find", "proyears", "optR0_Eff", "SSBabsolute"))
						
## Optimize Fishing Effort for each simulation ##
st <- Sys.time()
Out <- snowfall::sfSapply(1:nsim, OptFun, seeds, optOM, SSBabsolute, custompars, Find, proyears)
el <- Sys.time() - st
el # 

## takes about 3 mins per simulation on my 4 core machine 

Outpars <- joinSFOb(Out) 

## Save out object 
ParsObj <- list(custompars=Outpars, Stock=Stock, Fleet=Fleet, nsim=nsim, 
               proyears=proyears, SSBabsolute=SSBabsolute)

saveRDS(ParsObj, file=paste0(wd, "/SaveObjs/", Name, "_pars.rds"))























