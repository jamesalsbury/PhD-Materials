#Here is the R script used to produce the figures/results from the DTE paper

library(truncnorm)
library(survival)

# An example of DTE in a KM plot ----------------------------------------------------------------

png("DTEKM.png", units="in", width=7, height=6, res=700)

controldata <- read.csv(file = "Papers/DTEPaper/data/Brahmer/IPD-control.csv")
treatmentdata <- read.csv(file = "Papers/DTEPaper/data/Brahmer/IPD-treatment.csv")

combinedData <- data.frame(time = c(controldata$Survival.time, treatmentdata$Survival.time), 
                           status = c(controldata$Status, treatmentdata$Status), 
                           group = c(rep("Control", nrow(controldata)), rep("Treatment", nrow(treatmentdata))))

kmfit <- survfit(Surv(time, status)~group, data = combinedData)
plot(kmfit, conf.int = F, col=c("blue", "red"), xlab = "Time (months)", 
     ylab = "Progression free survival (% of patients)", xaxt = "n", yaxt = "n")
axis(1, at=seq(0, 21, by=3), labels=seq(0, 21, by=3))
axis(2, at=seq(0, 1, by=0.1), labels=seq(0, 100, by=10))
legend("topright", legend = c("Docetaxel", "Nivolumab"), col = c("blue", "red"), lty=1)

dev.off()

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
  
  controldata <- read.csv(file = paste0("Papers/DTEPaper/data/", author, "/IPD-control.csv"))
  treatmentdata <- read.csv(file = paste0("Papers/DTEPaper/data/", author, "/IPD-treatment.csv"))
  
  combinedData <- data.frame(time = c(controldata$Survival.time, treatmentdata$Survival.time), 
                             status = c(controldata$Status, treatmentdata$Status), 
                             group = c(rep("Control", nrow(controldata)), rep("Treatment", nrow(treatmentdata))))
  
  kmfit <- survfit(Surv(time, status)~group, data = combinedData)
  #plot(kmfit, conf.int = F, col=c("blue", "red"), xlab = "Time (months)", ylab = ylabel, yaxt = "n")
  #axis(2, at=seq(0, 1, by=0.1), labels=seq(0, 100, by=10))
  
  
  #Finding the Weibull parameters and plotting the line
  weibfit <- survreg(Surv(time, status)~1, data = combinedData[combinedData$group=="Control",], dist = "weibull")
  gammac <- as.numeric(exp(-weibfit$icoef[2]))
  lambdac <- as.numeric(1/(exp(weibfit$icoef[1])))
  controltime <- seq(0, trialLength, by=0.01)
  controlsurv <- exp(-(lambdac*controltime)^gammac)
  #lines(controltime, controlsurv, col="blue")
  
  #Now we look at the treatment curve
  #Need to find a least squares estimate for this Weibull parameterisation
  
  kmfit <- survfit(Surv(time, status)~1, data = combinedData[combinedData$group=="Treatment",])
  treatmenttime <- seq(0, THat, by=0.01)
  treatmentsurv <- controlsurv <- exp(-(lambdac*treatmenttime)^gammac)
  #lines(treatmenttime, treatmentsurv, col="green")
  
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
  #lines(treatmenttime1, treatmentsurv1, col="red", lty=1)
  
  
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
  #lines(treatmenttime1, treatmentsurv1, col="red", lty=2)
  
  #legend("topright", legend = c("Delay", "Control", expression(paste("Treatment (estimate ", hat(lambda)[t], ")")), expression(paste("Treatment (estimate ", tilde(lambda)[t], ",", tilde(gamma)[t], ")"))), col=c("green", "blue", "red", "red"), lty=c(1, 1, 1,2))
  
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


  nvec <- round(seq(10, 500, length = 30))
  powervecs1 <- mapply(powerfunc, nvec, s1lambdat, s1gammat)
  powervecs2 <- mapply(powerfunc, nvec, s2lambdat, s2gammat)

  powers1smooth <- loess(powervecs1~nvec)
  powers2smooth <- loess(powervecs2~nvec)

  plot(nvec*2, predict(powers1smooth), type="l", col="blue", ylim=c(0, 1), xlab = "Total sample size", ylab = "Power")
  lines(nvec*2, predict(powers2smooth), col="red", lty=2)

  legend("bottomright", legend = c("Scenario 1", "Scenario 2"), col = c("blue", "red"), lty=1:2)
}

png("DTEPower.png", units="in", width=15, height=5, res=700)
par(mfrow=c(1,3))
DTEDataSetsFunc("Brahmer")
DTEDataSetsFunc("Yen")
DTEDataSetsFunc("Borghaei")
dev.off()

# Sampling 10 treatment survival curves and CI ----------------------------------------------------------------
png("simSamples.png", units="in", width=14, height=6, res=700)
par(mfrow=c(1,2))
gammac <- 0.8
gammat <- 0.8
lambdac <- 0.08
trialtime <- seq(0, 100, by=0.01)
controlsurv <- exp(-(lambdac*trialtime)^gammac)
plot(trialtime, controlsurv, type="l", col="blue", xlab = "Time (months)", yaxt = "n", ylab = "Overall survival (%)")
axis(2, at=seq(0, 1, by=0.1), labels=seq(0, 100, by=10))
for (j in 1:10){
  bigT <- truncnorm::rtruncnorm(1, mean = 6, sd = 2.97, a = 0)
  HR <- rbeta(1, 10.8, 6.87)
  lambdat <- lambdac*HR^(1/gammac)
  treatmenttime <- ifelse(trialtime<=bigT, exp(-(lambdac*trialtime)^gammac), exp(-(lambdac*bigT)^gammac-lambdat^gammat*(trialtime^gammat-bigT^gammat)))
  lines(trialtime, treatmenttime, col="red")
}

legend("topright", legend = c("Control", "Treatment"), lty = 1, col=c("blue", "red"))

gammac <- 0.8
gammat <- 0.8
lambdac <- 0.08
trialtime <- seq(0, 100, by=0.01)
controlsurv <- exp(-(lambdac*trialtime)^gammac)
plot(trialtime, controlsurv, type="l", col="blue", xlab = "Time (months)", yaxt = "n", ylab = "Overall survival (%)")
axis(2, at=seq(0, 1, by=0.1), labels=seq(0, 100, by=10))
M <- 500
SimMatrix <- matrix(NA, nrow = M, ncol = length(trialtime))
for (i in 1:M){
  bigT <- truncnorm::rtruncnorm(1, mean = 6, sd = 2.97, a = 0)
  HR <- rbeta(1, 10.8, 6.87)
  lambdat <- lambdac*HR^(1/gammac)
  SimMatrix[i,] <- ifelse(trialtime<=bigT, exp(-(lambdac*trialtime)^gammac), exp(-(lambdac*bigT)^gammac-lambdat^gammat*(trialtime^gammat-bigT^gammat)))
}
lowerbound <- rep(NA, length(trialtime))
medianbound <- rep(NA, length(trialtime))
upperbound <- rep(NA, length(trialtime))
for (j in 1:length(trialtime)){
  lowerbound[j] <- quantile(SimMatrix[,j], 0.1)
  medianbound[j] <- quantile(SimMatrix[,j], 0.5)
  upperbound[j] <- quantile(SimMatrix[,j], 0.9)
}
lines(trialtime, medianbound, col="red")
lines(trialtime, lowerbound, lty=2)
lines(trialtime, upperbound, lty=2)
legend("topright", legend = c("Control", "Treatment", "Treatment CI"), lty = c(1, 1, 2), col=c("blue", "red", "black"))
dev.off()

# Showing an example of the more flexible Weibull ----------------------------------------------------------------
gammac <- 0.8
gammat <- 0.8
lambdac <- 0.08
trialtime <- seq(0, 100, by=0.01)
controlsurv <- exp(-(lambdac*trialtime)^gammac)
M <- 500
SimMatrix <- matrix(NA, nrow = M, ncol = length(trialtime))
for (i in 1:M){
  bigT <- truncnorm::rtruncnorm(1, mean = 6, sd = 2.97, a = 0)
  HR <- rbeta(1, 10.8, 6.87)
  lambdat <- lambdac*HR^(1/gammac)
  SimMatrix[i,] <- ifelse(trialtime<=bigT, exp(-(lambdac*trialtime)^gammac), exp(-(lambdac*bigT)^gammac-lambdat^gammat*(trialtime^gammat-bigT^gammat)))
}
  bigT <- truncnorm::rtruncnorm(1, mean = 6, sd = 2.97, a = 0)
  s1 <- sample(SimMatrix[,round(length(trialtime)*0.25)], 1)
  s2 <- 1
  while (s2>s1){ s2 <- sample(SimMatrix[,round(length(trialtime)*0.6)], 1)}
  estWeib <- function(par){
    t1 <- exp(-(lambdac*bigT)^gammac-par[1]^par[2]*(25^par[2]-bigT^par[2])) - s1
    t2 <- exp(-(lambdac*bigT)^gammac-par[1]^par[2]*(60^par[2]-bigT^par[2])) - s2
    terror <- t1^2+t2^2
    return(terror)
  }
  
  output <- optim(par = c(0, 1), fn = estWeib)
  lambdat <- output$par[1]
  gammat <- output$par[2]
  treatmentsurv <- ifelse(trialtime<=bigT, exp(-(lambdac*trialtime)^gammac), exp(-(lambdac*bigT)^gammac-lambdat^gammat*(trialtime^gammat-bigT^gammat)))

png("FlexSeq.png", units="in", width=14, height=12, res=700)
par(mfrow=c(2,2))
#First plot
plot(trialtime, controlsurv, type="l", col="blue", xlab = "Time (months)", yaxt = "n", ylab = "Overall survival (%)")
axis(2, at=seq(0, 1, by=0.1), labels=seq(0, 100, by=10))
points(bigT, exp(-(lambdac*bigT)^gammac), col = "red", bg = "red", pch = 19)
legend("topright", legend = "Control", lty = 1, col="blue")
#Second plot
plot(trialtime, controlsurv, type="l", col="blue", xlab = "Time (months)", yaxt = "n", ylab = "Overall survival (%)")
axis(2, at=seq(0, 1, by=0.1), labels=seq(0, 100, by=10))
points(bigT, exp(-(lambdac*bigT)^gammac), col = "red", bg = "red", pch = 19)
points(25, s1, col = "red", bg = "red", pch = 19)
legend("topright", legend = "Control", lty = 1, col="blue")
#Third plot
plot(trialtime, controlsurv, type="l", col="blue", xlab = "Time (months)", yaxt = "n", ylab = "Overall survival (%)")
axis(2, at=seq(0, 1, by=0.1), labels=seq(0, 100, by=10))
points(bigT, exp(-(lambdac*bigT)^gammac), col = "red", bg = "red", pch = 19)
points(25, s1, col = "red", bg = "red", pch = 19)
points(60, s2, col = "red", bg = "red", pch = 19)
legend("topright", legend = "Control", lty = 1, col="blue")
#Fourth plot
plot(trialtime, controlsurv, type="l", col="blue", xlab = "Time (months)", yaxt = "n", ylab = "Overall survival (%)")
axis(2, at=seq(0, 1, by=0.1), labels=seq(0, 100, by=10))
points(bigT, exp(-(lambdac*bigT)^gammac), col = "red", bg = "red", pch = 19)
points(25, s1, col = "red", bg = "red", pch = 19)
points(60, s2, col = "red", bg = "red", pch = 19)
lines(trialtime, treatmentsurv, col="red")
legend("topright", legend = c("Control", "Treatment"), lty = 1, col=c("blue", "red"))
dev.off()

# Sampling 10 flexible lines & then the CI ----------------------------------------------------------------
png("simSamplesFlex.png", units="in", width=14, height=6, res=700)
par(mfrow=c(1,2))
gammac <- 0.8
gammat <- 0.8
lambdac <- 0.08
trialtime <- seq(0, 100, by=0.01)
controlsurv <- exp(-(lambdac*trialtime)^gammac)
M <- 500
plot(trialtime, controlsurv, type="l", col="blue", xlab = "Time (months)", yaxt = "n", ylab = "Overall survival (%)")
axis(2, at=seq(0, 1, by=0.1), labels=seq(0, 100, by=10))
SimMatrix <- matrix(NA, nrow = M, ncol = length(trialtime))
for (i in 1:M){
  bigT <- truncnorm::rtruncnorm(1, mean = 6, sd = 2.97, a = 0)
  HR <- rbeta(1, 10.8, 6.87)
  lambdat <- lambdac*HR^(1/gammac)
  SimMatrix[i,] <- ifelse(trialtime<=bigT, exp(-(lambdac*trialtime)^gammac), exp(-(lambdac*bigT)^gammac-lambdat^gammat*(trialtime^gammat-bigT^gammat)))
}
for (j in 1:10){
  bigT <- truncnorm::rtruncnorm(1, mean = 6, sd = 2.97, a = 0)
  s1 <- sample(SimMatrix[,round(length(trialtime)*0.25)], 1)
  s2 <- 1
  while (s2>s1){ s2 <- sample(SimMatrix[,round(length(trialtime)*0.6)], 1)}
  estWeib <- function(par){
    t1 <- exp(-(lambdac*bigT)^gammac-par[1]^par[2]*(25^par[2]-bigT^par[2])) - s1
    t2 <- exp(-(lambdac*bigT)^gammac-par[1]^par[2]*(60^par[2]-bigT^par[2])) - s2
    terror <- t1^2+t2^2
    return(terror)
  }
  
  output <- optim(par = c(0, 1), fn = estWeib)
  lambdat <- output$par[1]
  gammat <- output$par[2]
  treatmentsurv <- ifelse(trialtime<=bigT, exp(-(lambdac*trialtime)^gammac), exp(-(lambdac*bigT)^gammac-lambdat^gammat*(trialtime^gammat-bigT^gammat)))
  lines(trialtime, treatmentsurv, col="red")
}
legend("topright", legend = c("Control", "Treatment"), lty = 1, col=c("blue", "red"))

FlexSimMatrix <- matrix(NA, nrow = M, ncol = length(trialtime))
for (i in 1:M){
  bigT <- truncnorm::rtruncnorm(1, mean = 6, sd = 2.97, a = 0)
  s1 <- sample(SimMatrix[,round(length(trialtime)*0.25)], 1)
  s2 <- 1
  while (s2>s1){ s2 <- sample(SimMatrix[,round(length(trialtime)*0.6)], 1)}
  estWeib <- function(par){
    t1 <- exp(-(lambdac*bigT)^gammac-par[1]^par[2]*(25^par[2]-bigT^par[2])) - s1
    t2 <- exp(-(lambdac*bigT)^gammac-par[1]^par[2]*(60^par[2]-bigT^par[2])) - s2
    terror <- t1^2+t2^2
    return(terror)
  }
  
  output <- optim(par = c(0, 1), fn = estWeib)
  lambdat <- output$par[1]
  gammat <- output$par[2]
  FlexSimMatrix[i,] <- ifelse(trialtime<=bigT, exp(-(lambdac*trialtime)^gammac), exp(-(lambdac*bigT)^gammac-lambdat^gammat*(trialtime^gammat-bigT^gammat)))
}

lowerbound <- rep(NA, length(trialtime))
upperbound <- rep(NA, length(trialtime))
medianbound <- rep(NA, length(trialtime))
for (i in 1:length(trialtime)){
  lowerbound[i] <- quantile(FlexSimMatrix[,i], 0.1, na.rm = T)
  upperbound[i] <- quantile(FlexSimMatrix[,i], 0.9, na.rm = T)
  medianbound[i] <- quantile(FlexSimMatrix[,i], 0.5, na.rm = T)
}

plot(trialtime, controlsurv, type="l", col="blue", xlab = "Time (months)", yaxt = "n", ylab = "Overall survival (%)")
axis(2, at=seq(0, 1, by=0.1), labels=seq(0, 100, by=10))
lines(trialtime, medianbound, col="red")
lines(trialtime, lowerbound, lty=2)
lines(trialtime, upperbound, lty=2)
legend("topright", legend = c("Control", "Treatment", "Treatment CI"), lty = c(1, 1, 2), col = c("blue", "red", "black"))
dev.off()
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

png("PowerAss.png", units="in", width=8, height=5, res=700)


simulateDTEWeibullData <- function(n1, n2, gammat, gammac, lambdat, lambdac, bigT, recTime, censTime){
  #Simulates the treatment data
  CP <- exp(-(lambdac*bigT)^gammac)
  u <- runif(n2)
  
  treatmenttime <- ifelse(u>CP, (1/lambdac)*(-log(u))^(1/gammac), (1/(lambdat^gammat)*((lambdat*bigT)^gammat-log(u)-(lambdac*bigT)^gammac))^(1/gammat))
  
  dataCombined <- data.frame(time = c(rweibull(n1, gammac, 1/lambdac), treatmenttime),
                             group = c(rep("Control", n1), rep("Treatment", n2)))
  
  
  #Adds a random uniformly distributed value, based on the recruitment time
  dataCombined$time <- dataCombined$time + runif(n1+n2, min = 0, max = recTime)
  
  #If the time is less than the total trial length time then the event has happened
  dataCombined$event <- dataCombined$time < censTime
  
  #Making it a binary value (rather than T/F), for ease to read
  dataCombined$event <- dataCombined$event*1
  
  #Need to set the event time to be the censoring time
  if (sum(dataCombined$event)==(n1+n2)){
    
  } else{
    dataCombined[dataCombined$time>censTime,]$time <- censTime
  }
  
  return(dataCombined)
}

powerassFunc <- function(type, n){
  gammac <- 0.8
  lambdac <- 0.08
  recTime <- 6
  Lmax <- 36
  massT0 <- 0.05
  massHR1 <- 0.1
  N <- 500
  #Making the simplification
  gammat <- gammac
  #Setting up the assurance vector
  vec <- rep(NA, N)
  for (i in 1:N){
    #sampling the delay time
    if (type=="assurance"){
      u <- runif(1)
      if (u < massT0){
        bigT <- 0
      } else {
        bigT <- truncnorm::rtruncnorm(1, mean = 6, sd = 2.97, a = 0)
      }
      #sampling the post-delay HR
      u <- runif(1)
      if (u < massHR1){
        HR <- 1
      } else {
        HR <- rbeta(1, 10.8, 6.87)
      }
    } else if (type=="power"){
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
    } else if (type=="powerND"){
      bigT <- 0
      #sampling the post-delay HR
      u <- runif(1)
      if (u < massHR1){
        HR <- 1
      } else {
        HR <- 0.6
      }
    }
    
    lambdat <- exp((log(HR)/gammac)+log(lambdac))
    
    #Simulating the control and treatment data
    combinedData <- simulateDTEWeibullData(n, n, gammat, gammac, lambdat, lambdac, bigT, recTime, Lmax)
    
    #Performing a log-rank test on the combined data set
    test <- survdiff(Surv(time, event)~group, data = combinedData)
    vec[i] <- test$chisq > qchisq(0.95, 1)
  }
  return(mean(vec))
}


calcAssPowerFunc <- function(type){
  nvec <- seq(20, 500, by=10)
  output <- sapply(X = nvec, FUN = powerassFunc, type= type)
  smoothedout <- loess(output~nvec)
  return(list(smoothedout = smoothedout, nvec = nvec))
  
}

ass <- calcAssPowerFunc("assurance")
power <- calcAssPowerFunc("power")
powerND <- calcAssPowerFunc("powerND")

plot(power$nvec*2, predict(power$smoothedout), ylim=c(0,1), type="l", xlab = "Total sample size", ylab = "Power/assurance")
lines(ass$nvec*2, predict(ass$smoothedout), col="red", lty=2)
lines(powerND$nvec*2, predict(powerND$smoothedout), col="blue", lty=3)

legend("bottomright", legend = c("Power", "Assurance", "Power assuming no delay"), col=c("black", "red", "blue"), lty=1:3)

dev.off()


sum(is.na(predict(ass$smoothedout, 1:500)<0.7))+sum(na.omit(predict(ass$smoothedout, 1:500)<0.7))*2
sum(is.na(predict(power$smoothedout, 1:500)<0.7))+sum(na.omit(predict(power$smoothedout, 1:500)<0.7))*2
sum(is.na(predict(powerND$smoothedout, 1:500)<0.7))+sum(na.omit(predict(powerND$smoothedout, 1:500)<0.7))*2


