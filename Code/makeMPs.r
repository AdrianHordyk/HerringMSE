
# Management Procedures for CA Herring 
# A. Hordyk 
# Jan 2017

## Still to do: 
# Add maximum tonnage cap where required 
# come up with some better naming conventions 
# add variations including jump to 200 tons when minB threshold is
#     passed and a cap on the total tonnage if biomass is over 100,000 tons 
# add attainment factor? realized catch can be a fraction of TAC


## Things to check:
## is abundance spawning biomass or vulnerable biomass? HCR assumes that abundance is VBiomass
## 
## check that weight units in MSE must match those in MPs - 
## 

# library(DLMtool)
# testdata <- new("DLM_data") # create an empty data object 
# testdata@SpAbun <- 10000 # need this to initialize functions 

###################################
## First Family of Control Rules ##
###################################

#' MP Factory
#' Fixed harvest rate with a linear control rule below maxB 
#' Test various levels of minB and maxB 
#'
#' Closure to create MP function with minB and maxB 
#' Fixed harvest rate of 10%, with linear decrease if 
#' abundance is below maxB down to 0 if below minB 
makeMinMaxMP <- function(minB, maxB, PercCap=0.1) {
  tempFun <- function(x, DLM_data, reps=100) {
    # estimate of current abundance
    # Abun <-  DLMtool::trlnorm(reps, DLM_data@Abun[x], DLM_data@CV_Abun[x]) 
	Abun <-  DLMtool::trlnorm(reps, DLM_data@SpAbun[x], DLM_data@CV_Abun[x]) # spawning stock biomass
    TAC <- rep(NA, reps)
    #Calc slope of line
    m <- PercCap/(maxB - minB)  
    #Calc intercept of line
    b <- PercCap - (m * maxB)
    
    PercQuota <- rep(PercCap, reps)
    TAC <- PercQuota * Abun
    TAC[Abun <= minB] <- 0 
    TAC[Abun > minB & Abun < maxB] <- Abun[Abun > minB & Abun < maxB] * 
                                      (m * Abun[Abun > minB & Abun < maxB] + b)
    DLMtool::TACfilter(TAC)
  }
  class(tempFun) <- "DLM_output"
  attr(tempFun, "case") <- "CAHerring"
  tempFun 
}

# Example MPs with different minB and maxB
# slope10_60 <- makeMinMaxMP(minB=10000, maxB=60000)

# slope20_80 <- makeMinMaxMP(minB=20000, maxB=80000)

# slope40_100 <- makeMinMaxMP(minB=40000, maxB=100000)


# Very easy to create hundreds of MPs # 

# range of alternative values - example only 
# minBs <- seq(5000, 30000, by=1000)
# maxBs <- seq(50000, 100000, by=1000)
# Grid <- expand.grid(minBs=minBs, maxBs=maxBs)

# Create nrow(Grid) MP with different combinations of 
# minB and maxB  
# for (X in 1:nrow(Grid)) {
  # name <- paste0("slope", Grid$minBs[X]/1000, "_", Grid$maxBs[X]/1000)
  # tempFun <- makeMinMaxMP(Grid$minBs[X], Grid$maxBs[X])
  # assign(name, tempFun) 
  # # Need this line to run function once
  # get(name)(1, testdata, reps=1) 
  # # otherwise weird stuff happens! 
# }

####################################
## Second Family of Control Rules ##
####################################

#' MP Factory
#' Fixed escapement and fixed harvest rate with TAC = 0 below minimum biomass 
#' Test various levels of harvest rate, and minimum biomass
#'
#' Closure to create MP function with minB and PercQuota 
#'
makeFixEsc <- function(minB, PercQuota) {
  tempFun <- function(x, DLM_data, reps=100) {
    # estimate of current abundance
    # Abun <- DLMtool::trlnorm(reps, DLM_data@Abun[x], DLM_data@CV_Abun[x]) 
	Abun <-  DLMtool::trlnorm(reps, DLM_data@SpAbun[x], DLM_data@CV_Abun[x]) # spawning stock biomass	
    TAC <- rep(NA, reps)
    
    TAC[Abun <= minB] <- 0 
    TAC[Abun > minB] <- (Abun[Abun > minB] - minB) * PercQuota
                                   
    DLMtool::TACfilter(TAC)
  }
  class(tempFun) <- "DLM_output"
  attr(tempFun, "case") <- "CAHerring"
  tempFun
}

# saradine 
# # Example MPs 
# esc15_0.07 <- makeFixEsc(15000, 0.07)

# esc20_0.1 <- makeFixEsc(20000, 0.1)

# esc5_0.05 <- makeFixEsc(5000, 0.05)

# # range of alternative values - example only 
# esc <- seq(5000, 50000, by=5000)
# hrate <- seq(0.01, 0.1, 0.01)
# Grid <- expand.grid(esc=esc, hrate=hrate)

# Create nrow(Grid) MP with different combinations of 
# escapement and fixed harvest rate 
# for (X in 1:nrow(Grid)) {
  # name <- paste0("esc", Grid$esc[X]/1000, "_", Grid$hrate[X])
  # tempFun <- makeFixEsc(Grid$esc[X], Grid$hrate[X])
  # assign(name, tempFun)
  # # Need this line to run function once
  # get(name)(1, testdata, reps=1) 
  # # otherwise weird stuff happens!  
# }

####################################
## Third Family of Control Rules ##
####################################

#' MP Factory
#' Switch between two fixed harvest rates 
#' 
#' hr1 fixed harvest rate 1 
#' hr2 fixed harvest rate 2 
#' switchB switch biomass
#' minB minimum biomass below which TAC = 0 
makeSwitch <- function(hr1, hr2, switchB, minB) {
  tempFun <- function(x, DLM_data, reps=100) {
    # estimate of current abundance
    # Abun <- DLMtool::trlnorm(reps, DLM_data@Abun[x], DLM_data@CV_Abun[x]) 
	Abun <-  DLMtool::trlnorm(reps, DLM_data@SpAbun[x], DLM_data@CV_Abun[x]) # spawning stock biomass	
    TAC <- rep(NA, reps)  
    TAC[Abun <= minB] <- 0 
	TAC[Abun > switchB] <- Abun[Abun > switchB] * hr2
	TAC[Abun > minB & Abun <= switchB] <- Abun[Abun > minB & Abun <= switchB] * hr1
     DLMtool::TACfilter(TAC)
  }
  class(tempFun) <- "DLM_output"
  attr(tempFun, "case") <- "CAHerring"
  tempFun 
}


# swch0.03_0.07_40_5 <- makeSwitch(0.03, 0.07, 40000, 5000)

# # Make a couple of examples 
# # range of alternative values - example only 
# hr1s <- seq(0.01, 0.05, by=0.025)
# hr2s <- seq(0.05, 0.1, 0.025)
# switchBs <- seq(50000, 60000, by=10000)
# minBs <- seq(5000, 7000, 2000)

# Grid <- expand.grid(hr1s=hr1s, hr2s=hr2s, switchBs=switchBs, minBs=minBs)

# # Create nrow(Grid) MP with different combinations of 
# # escapement and fixed harvest rate 
# for (X in 1:nrow(Grid)) {
  # name <- paste0("swch", Grid$hr1s[X], "_", Grid$hr2s[X], "_", 
    # Grid$switchBs[X]/1000, "_", Grid$minBs[X]/1000)
  # tempFun <- makeSwitch(Grid$hr1s[X], Grid$hr2s[X], Grid$switchBs[X], Grid$minBs[X])
  # assign(name, tempFun)
  # # Need this line to run function once
  # get(name)(1, testdata, reps=1) 
  # # otherwise weird stuff happens!  
# }





