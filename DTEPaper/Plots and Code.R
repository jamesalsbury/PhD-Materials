#Here is the R script used to produce the figures/results from the DTE paper

#This finds the assurance for the example in Section 2.1

library(truncnorm)
library(survival)


# Figure 1 ----------------------------------------------------------------

#png("DTEKM.png", units="in", width=5, height=5, res=700)

controldata <- read.csv(file = "DTEPaper/data/Brahmer/IPD-control.csv")
treatmentdata <- read.csv(file = "DTEPaper/data/Brahmer/IPD-treatment.csv")

combinedData <- data.frame(time = c(controldata$Survival.time, treatmentdata$Survival.time), 
                           status = c(controldata$Status, treatmentdata$Status), 
                           group = c(rep("Control", nrow(controldata)), rep("Treatment", nrow(treatmentdata))))

kmfit <- survfit(Surv(time, status)~group, data = combinedData)
plot(kmfit, conf.int = F, col=c("blue", "red"), xlab = "Time (months)", ylab = "Progression free survival (% of patients)")
legend("topright", legend = c("Control", "Treatment"), col = c("blue", "red"), lty=1)

#dev.off()

# Robustness of the parameterisation ----------------------------------------------------------------

DTEDataSetsFunc <- function(author){
  
  if (author=="Brahmer"){
    trialLength <- 15
    THat <- 3
    ylabel <- "Progression free survival (% of patients)"
  } else if (author=="Yen"){
    trialLength <- 30
    THat <- 3.5
    ylabel <- "Overall survival (%)"
  } else if (author=="Borghaei"){
    trialLength <- 60
    THat <- 6
    ylabel <- "Overall survival (%)"
  }
  
  controldata <- read.csv(file = paste0("DTEPaper/data/", author, "/IPD-control.csv"))
  treatmentdata <- read.csv(file = paste0("DTEPaper/data/", author, "/IPD-treatment.csv"))
  
  combinedData <- data.frame(time = c(controldata$Survival.time, treatmentdata$Survival.time), 
                             status = c(controldata$Status, treatmentdata$Status), 
                             group = c(rep("Control", nrow(controldata)), rep("Treatment", nrow(treatmentdata))))
  
  kmfit <- survfit(Surv(time, status)~group, data = combinedData)
  plot(kmfit, conf.int = F, col=c("blue", "red"), xlab = "Time (months)", ylab = ylabel)
  
  
  #Finding the Weibull parameters and plotting the line
  weibfit <- survreg(Surv(time, status)~1, data = combinedData[combinedData$group=="Control",], dist = "weibull")
  gammac <- as.numeric(exp(-weibfit$icoef[2]))
  lambdac <- as.numeric(1/(exp(weibfit$icoef[1])))
  controltime <- seq(0, trialLength, by=0.01)
  controlsurv <- exp(-(lambdac*controltime)^gammac)
  lines(controltime, controlsurv, col="blue")
  
  #Now we look at the treatment curve
  #Need to find a least squares estimate for this Weibull parameterisation
  
  kmfit <- survfit(Surv(time, status)~1, data = combinedData[combinedData$group=="Treatment",])
  treatmenttime <- seq(0, THat, by=0.01)
  treatmentsurv <- controlsurv <- exp(-(lambdac*treatmenttime)^gammac)
  lines(treatmenttime, treatmentsurv, col="green")
  
  treatmenttime1 <- seq(THat, trialLength, by=0.01)
  
  
  #Scenario 1
  
  optimfunc1 <- function(par){
    diff <- 0
    gammat <- gammac
    for (i in 1:length(treatmenttime1)){
      y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
      treatmentsurv <- exp(-(lambdac*THat)^gammac - par[1]^gammat*(treatmenttime1[i]^gammat-THat^gammat))
      diff <- diff + (y-treatmentsurv)^2 
    }
    return(diff) 
  }
  
  s1gammat <- gammac
  s1lambdat <- optimize(optimfunc1, c(0, 2))$minimum
  treatmentsurv1 <- exp(-(lambdac*THat)^gammac - s1lambdat^s1gammat*(treatmenttime1^s1gammat-THat^s1gammat))
  lines(treatmenttime1, treatmentsurv1, col="red", lty=1)
  
  
  #Scenario 2
  
  
  optimfunc2 <- function(par){
    diff <- 0
    for (i in 1:length(treatmenttime1)){
      y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
      treatmentsurv <- exp(-(lambdac*THat)^gammac - par[1]^par[2]*(treatmenttime1[i]^par[2]-THat^par[2]))
      diff <- diff + (y-treatmentsurv)^2 
    }
    return(diff) 
  }
  
  optimoutput <-  optim(par = c(0.2, 1), fn = optimfunc2)
  s2lambdat <- optimoutput$par[1]
  s2gammat <- optimoutput$par[2]
  treatmentsurv1 <- exp(-(lambdac*THat)^gammac - s2lambdat^s2gammat*(treatmenttime1^s2gammat-THat^s2gammat))
  lines(treatmenttime1, treatmentsurv1, col="red", lty=2)
  
  legend("topright", legend = c("Delay", "Control", expression(paste("Treatment (estimate ", hat(lambda)[t], ")")), expression(paste("Treatment (estimate ", tilde(lambda)[t], ",", tilde(gamma)[t], ")"))), col=c("green", "blue", "red", "red"), lty=c(1, 1, 1,2))
  
  #Now look at calculating power under the two different scenarios
  #Scenario 1
  
  powerfunc <- function(n, lambdat, gammat){
    powervec <- rep(NA, 1000)
    for (i in 1:length(powervec)){
      #Simulating control times
      u <- runif(n)
      controltimes <- (1/lambdac)*(-log(u))^(1/gammac)

      #Simulating treatment times
      u <- runif(n)
      CP <- exp(-(lambdac*THat)^gammac)
      treatmenttimes <- ifelse(u>=CP, (1/lambdac)*(-log(u))^(1/gammac), (1/(lambdat^gammat)*((lambdat*THat)^gammat-log(u)-(lambdac*THat)^gammac))^(1/gammat))

      #combining control and treatment
      combinedDataPower <- data.frame(time = c(controltimes, treatmenttimes),
                                      status = rep(1, 2*n), group = c(rep("Control", n), rep("Treatment", n)))

      test <- survdiff(Surv(time, status)~group, data = combinedDataPower)
      #If the p-value of the test is less than 0.05 then assvec = 1, 0 otherwise
      powervec[i] <- test$chisq > qchisq(0.95, 1)
    }

    return(mean(powervec))
  }

  #png(paste0("DTE", author, "power.png"), units="in", width=5, height=5, res=700)
  nvec <- round(seq(10, 500, length = 30))
  powervecs1 <- mapply(powerfunc, nvec, s1lambdat, s1gammat)
  powervecs2 <- mapply(powerfunc, nvec, s2lambdat, s2gammat)

  powers1smooth <- loess(powervecs1~nvec)
  powers2smooth <- loess(powervecs2~nvec)

  plot(nvec*2, predict(powers1smooth), type="l", col="blue", ylim=c(0, 1), xlab = "Total sample size", ylab = "Power")
  lines(nvec*2, predict(powers2smooth), col="red", lty=2)

  legend("bottomright", legend = c("Scenario 1", "Scenario 2"), col = c("blue", "red"), lty=1:2)
  #dev.off()
}

DTEDataSetsFunc("Brahmer")
DTEDataSetsFunc("Yen")
DTEDataSetsFunc("Borghaei")

# Looking at the hazard function ----------------------------------------------------------------

gammac <- 0.8
lambdac <- 0.08
t <- seq(0, 100, by=0.01)
HFC <- gammac*lambdac^gammac*t^(gammac-1)
plot(t, HFC, type="l", ylim=c(0, 0.5), col="blue", ylab = "Hazard function", xlab = "Time")

gammat <-0.8
lambdat <- 0.04
HFT <- gammat*lambdat^gammat*t^(gammat-1)
lines(t, HFT, col="red")

# HR <- HFT/HFC
# lines(t, HR)

legend("topright", legend = c("Control", "Treatment"), col=c("blue", "red"), lty=1)


# Power/assurance figure for the example ----------------------------------------------------------------

#png("PowerAss.png", units="in", width=8, height=5, res=700)
assFunc <- function(nc, nt, N, recTime, Lmax, lambdac, gammac, massT0, massHR1){
  #Making the simplification
  gammat <- gammac
  #Setting up the assurance vector
  assvec <- rep(NA, N)
  for (i in 1:N){
    #sampling the delay time
    u <- runif(1)
    if (u < massT0){
      bigT <- 0
    } else {
      bigT <- rgamma(1, 5.76, 0.899)
    }
    #sampling the post-delay HR
    u <- runif(1)
    if (u < massHR1){
      HR <- 1
    } else {
      HR <- rbeta(1, 10.8, 6.87)
    }
    lambdat <- exp((log(HR)/gammac)+log(lambdac))
     
    #Simulating the control and treatment data
    combinedData <- simulateDTEWeibullData(nc, nt, gammat, gammac, lambdat, lambdac, bigT, recTime, Lmax)
   
    #Performing a log-rank test on the combined data set
    test <- survdiff(Surv(time, event)~group, data = combinedData)
    assvec[i] <- test$chisq > qchisq(0.95, 1)
  }
  return(mean(assvec))
}

ncvec <- seq(20, 500, by=10)
ntvec <- seq(20, 500, by=10)
#Setting up the assurance vector
finalassvec<- rep(NA, length(ncvec))
#Finding the assurance for varying sample sizes
for (j in 1:length(finalassvec)){
  finalassvec[j] <- assFunc(ncvec[j], ntvec[j], 500, 6, 36, 0.08, 0.8, 0.05, 0.1)
}

samplesizevec <- ncvec+ntvec

#Smoothing the output
asssmooth <- loess(finalassvec~samplesizevec)

#Plotting the output
plot(samplesizevec, predict(asssmooth), ylim=c(0,1), 
     type = "l", xlab="Total sample size", ylab="Power/assurance", col="red", lty=2)


#Calculating the power first
powerFunc <- function(nc, nt, N, recTime, Lmax, lambdac, gammac, massT0, massHR1){
  #Making the simplification
  gammat <- gammac
  #Setting up the assurance vector
  powervec <- rep(NA, N)
  for (i in 1:N){
    
    #sampling the delay time
    u <- runif(1)
    if (u < massT0){
      bigT <- 0
    } else {
      bigT <- 6
    }
    #sampling the post-delay HR
    u <- runif(1)
    if (u < massHR1){
      HR <- 1
    } else {
      HR <- 0.6
    }
    lambdat <- exp((log(HR)/gammac)+log(lambdac))
    
    #Simulating the control and treatment data
    combinedData <- simulateDTEWeibullData(nc, nt, gammat, gammac, lambdat, lambdac, bigT, recTime, Lmax)
    
    #Performing a log-rank test on the combined data set
    test <- survdiff(Surv(time, event)~group, data = combinedData)
    powervec[i] <- test$chisq > qchisq(0.95, 1)
  }
  #The mean of the assurance vector is the estimated assurance for the given parameters
  return(mean(powervec))
}

ncvec <- seq(20, 500, by=10)
ntvec <- seq(20, 500, by=10)
#Setting up the assurance vector
finalpowervec <- rep(NA, length(ncvec))
#Finding the assurance for varying sample sizes
for (j in 1:length(finalpowervec)){
  finalpowervec[j] <- powerFunc(ncvec[j], ntvec[j], 500, 6, 36, 0.08, 0.8, 0.05, 0.1)
}

samplesizevec <- ncvec+ntvec

#Smoothing the output
powersmooth <- loess(finalpowervec~samplesizevec)

#Plotting the output
lines(samplesizevec, predict(powersmooth),
     type = "l", col="black")


#Now we look at the power scenario where we do not allow for any delay

#Calculating the power first
powerNoDelayFunc <- function(nc, nt, N, recTime, Lmax, lambdac, gammac, massHR1){
  #Making the simplification
  gammat <- gammac
  #Setting up the assurance vector
  powervec <- rep(NA, N)
  for (i in 1:N){
      
    bigT <- 0
    
    #sampling the post-delay HR
    u <- runif(1)
    if (u < massHR1){
      HR <- 1
    } else {
      HR <- 0.6
    }
    lambdat <- exp((log(HR)/gammac)+log(lambdac))
    
    #Simulating the control and treatment data
    combinedData <- simulateDTEWeibullData(nc, nt, gammat, gammac, lambdat, lambdac, bigT, recTime, Lmax)
    
    #Performing a log-rank test on the combined data set
    test <- survdiff(Surv(time, event)~group, data = combinedData)
    powervec[i] <- test$chisq > qchisq(0.95, 1)
  }
  #The mean of the assurance vector is the estimated assurance for the given parameters
  return(mean(powervec))
}

ncvec <- seq(20, 500, by=10)
ntvec <- seq(20, 500, by=10)
#Setting up the assurance vector
finalpowernodelayvec <- rep(NA, length(ncvec))
#Finding the assurance for varying sample sizes
for (j in 1:length(finalpowernodelayvec)){
  finalpowernodelayvec[j] <- powerNoDelayFunc(ncvec[j], ntvec[j], 500, 6, 36, 0.08, 0.8, 0.1)
}

samplesizevec <- ncvec+ntvec

#Smoothing the output
nodelaypowersmooth <- loess(finalpowernodelayvec~samplesizevec)

lines(samplesizevec, predict(nodelaypowersmooth),
      type = "l", col="blue", lty=3)

legend("bottomright", legend = c("Power", "Assurance", "Power assuming no delay"), col=c("black", "red", "blue"), lty=1:3)

#dev.off()








