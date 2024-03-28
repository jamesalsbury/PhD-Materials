#Here is the R script used to produce the figures/results from the DTE paper

library(truncnorm)
library(survival)
library(rjags)
library(nph)

# Code to create Figure 1 ----------------------------------------------------------------

#png("Figure1.png", units="in", width=7, height=6, res=700)

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

#dev.off()

# Code to create Figure 6 ----------------------------------------------------------------

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
 plot(kmfit, conf.int = F, col=c("blue", "red"), xlab = "Time (months)", ylab = ylabel, yaxt = "n")
axis(2, at=seq(0, 1, by=0.2), labels=seq(0, 100, by=20))
  
  
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
 lines(treatmenttime, treatmentsurv, col="red")
  
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
 
 #Plotting horizontal lines to indicate the delay
 y <- seq(0, 1, by=0.01)
 x <- rep(THat, length(y))
 lines(x,y, lty = 2)
 text(x = THat, y = 0.9, labels = paste0("Delay = ", THat, "months"), pos = 4)
  
 legend("topright", legend = c("Control", "Method A", "Method B"), col=c("blue", "red", "red"), lty=c(1, 1,2))
  
}


#png("Figure6.png", units="in", width=15, height=5, res=700)
par(mfrow=c(1,3))
DTEDataSetsFunc("Brahmer")
DTEDataSetsFunc("Yen")
DTEDataSetsFunc("Borghaei")
#dev.off()

# Code to  ----------------------------------------------------------------
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


# Calculating the control distribution (through MCMC) for the example ----------------------------------------------------------------
png("ExampleCombinedControl.png", units="in", width=8, height=5, res=700)

Herbst <- read.csv(file = "Papers/DTEPaper/data/Herbst/Doce2.csv")
kmfit1 <- survfit(Surv(Survival.time, Status)~1, data = Herbst)
plot(kmfit1, col = "blue", conf.int = F, xlab = "Time (months)", ylab = "Overall Survival")

Garon <- read.csv(file = "Papers/DTEPaper/data/Garon/Doce1.csv")
kmfit2 <- survfit(Surv(Survival.time, Status)~1, data = Garon)
lines(kmfit2, col = "red", conf.int = F)

Kim <- read.csv(file = "Papers/DTEPaper/data/Kim/Doce3.csv")
kmfit3 <- survfit(Surv(Survival.time, Status)~1, data = Kim)
lines(kmfit3, col = "yellow", conf.int = F)

legend("topright", legend = c("ZODIAC", "REVEL", "INTEREST"), lty = 1, col = c("blue", "red", "yellow"))
dev.off()
combinedDoce <- rbind(Herbst, Garon, Kim)

kmfit4 <- survfit(Surv(Survival.time, Status)~1, data = combinedDoce)
lines(kmfit4, col = "black", conf.int = F)


#Performing MCMC on this data set

modelstring="

data {
  for (j in 1:n){
    zeros[j] <- 0
  }
}

model {
  C <- 10000
  for (i in 1:n){
    zeros[i] ~ dpois(zeros.mean[i])
    zeros.mean[i] <-  -l[i] + C
    l[i] <- ifelse(datEvent[i]==1, log(gamma2)+gamma2*log(lambda2*datTimes[i])-(lambda2*datTimes[i])^gamma2-log(datTimes[i]), -(lambda2*datTimes[i])^gamma2)
  }
  
    lambda2 ~ dnorm(1,1/10000)T(0,)
    gamma2 ~ dnorm(1,1/10000)T(0,)
    
    }
"

model = jags.model(textConnection(modelstring), data = list(datTimes = combinedDoce$Survival.time, datEvent = combinedDoce$Status, n= nrow(combinedDoce)), quiet = T) 

update(model, n.iter=1000)
output=coda.samples(model=model, variable.names=c("lambda2", "gamma2"), n.iter = 10000)

plot(output)

lambda2sample <- as.numeric(unlist(output[,2]))
gamma2sample <- as.numeric(unlist(output[,1]))

#MCMCSample <- data.frame(shape = gamma2sample, scale = lambda2sample)

#write.csv(MCMCSample, "MCMCSample.csv", row.names=FALSE)



weibfit <- survreg(Surv(Survival.time, Status)~1, data = combinedDoce, dist = "weibull")
fixedgammac <- as.numeric(exp(-weibfit$icoef[2]))
fixedlambdac <- as.numeric(1/(exp(weibfit$icoef[1])))


# Power/assurance figure for the example ----------------------------------------------------------------


simulateDTEWeibullData <- function(n1, n2, gammat, gammac, lambdat, lambdac, bigT, recTime, eventRate){
  #Simulates the treatment data
  CP <- exp(-(lambdac*bigT)^gammac)
  u <- runif(n2)
  
  treatmenttime <- ifelse(u>CP, (1/lambdac)*(-log(u))^(1/gammac), (1/(lambdat^gammat)*((lambdat*bigT)^gammat-log(u)-(lambdac*bigT)^gammac))^(1/gammat))
  
  dataCombined <- data.frame(time = c(rweibull(n1, gammac, 1/lambdac), treatmenttime),
                             group = c(rep("Control", n1), rep("Treatment", n2)))
  
  
  #Samples a recruitment time for each patient, from a Uniform distribution
  dataCombined$recTime <- runif(n1+n2, min = 0, max = recTime)
  
  #Calculates a pseudo event time
  dataCombined$pseudo_time <- dataCombined$time + dataCombined$recTime
  
  #Calculate the number of events required (comes from eventRate)
  numEventsRequired <- floor((n1+n2)*eventRate)
  
  #Work at at what time point this number of events have occurred
  censTime <- sort(dataCombined$pseudo_time)[numEventsRequired]
  
  #Censor the observations which will not have occurred by censTime
  dataCombined$status <- dataCombined$pseudo_time <= censTime
  dataCombined$status <- dataCombined$status*1
  
  #Ensure only the patients enrolled by the censoring time are included in the data set
  dataCombined$enrolled <- dataCombined$recTime <= censTime
  dataCombined <-  dataCombined[dataCombined$enrolled==T,]
  
  #For the patients that are censored, their survival time is the length of time they have been enrolled in the trial for
  dataCombined$survival_time <- ifelse(dataCombined$pseudo_time>censTime, censTime  - dataCombined$recTime, dataCombined$time)
  
  return(list(dataCombined = dataCombined, censTime = censTime))
}

P_E <- 0.9
P_DTE <- 0.7
#Doing the simMatrix (if performing the more flexible assurance)
  LMax <- 100
  trialtime <- seq(0, LMax, by=0.01) 
  SimMatrix <- matrix(NA, ncol = length(trialtime), nrow = 500)
  for (j in 1:500){
    gammac <- sample(gamma2sample, 1)
    lambdac <- sample(lambda2sample, 1)
    
    u <- runif(1)
    if (u>P_E){
      bigT <- 0
      HR <- 1
    } else {
      HR <- rgamma(1, 29.6, 47.8)
      bigT <- ifelse(runif(1)>P_DTE,  0, rgamma(1, 7.29, 1.76))
    }
    
    lambdat <- lambdac*HR^(1/gammac)
    gammat <- gammac
    
    SimMatrix[j,] <- ifelse(trialtime<=bigT, exp(-(lambdac*trialtime)^gammac), exp(-(lambdac*bigT)^gammac-lambdat^gammat*(trialtime^gammat-bigT^gammat)))
  }

powerassFunc <- function(type, n){
  recTime <- 12
  N <- 500
  
  #Setting up the assurance vector
  vec <- rep(NA, N)
  censvec <- rep(NA, N)
  for (i in 1:N){
    if (type=="assurance"){
      #Sampling the control parameters
      gammac <- sample(gamma2sample, 1)
      lambdac <- sample(lambda2sample, 1)
      #Making the simplification
      gammat <- gammac
      #Sampling the elicited parameters (T and HR)
      u <- runif(1)
      if (u>P_E){
        bigT <- 0
        HR <- 1
      } else {
        HR <- rgamma(1, 29.6, 47.8)
        bigT <- ifelse(runif(1)>P_DTE,  0, rgamma(1, 7.29, 1.76))
      }
      lambdat <- lambdac*HR^(1/gammac)
      
    } else if (type=="power"){
      gammac <- fixedgammac
      lambdac <- fixedlambdac
      bigT <- 4
      HR <- 0.6
      #Making the simplification
      gammat <- gammac
      lambdat <- lambdac*HR^(1/gammac)
    } else if (type=="powerND"){
      gammac <- fixedgammac
      lambdac <- fixedlambdac
      bigT <- 0
      HR <- 0.6
      #Making the simplification
      gammat <- gammac
      lambdat <- lambdac*HR^(1/gammac)
    } else if (type=="flexAss"){
      
      gammac <- sample(gamma2sample, 1)
      lambdac <- sample(lambda2sample, 1)
      trialend <- exp((1/gammac)*log(-log(0.01)) -log(lambdac))
      
      time1 <- trialtime[which.min(abs(trialtime-0.35*trialend))]
      time2 <- trialtime[which.min(abs(trialtime-0.55*trialend))]
      
      s1 <- sample(SimMatrix[,which.min(abs(trialtime-0.35*trialend))], 1)
      s2 <- 1
      while (s2>s1){ 
        s2 <- sample(SimMatrix[,which.min(abs(trialtime-0.55*trialend))], 1)
      }
      
      bigT <- ifelse(runif(1)>P_DTE,  0, rgamma(1, 7.29, 1.76))
      
      estWeib <- function(par){
        t1 <- exp(-(lambdac*bigT)^gammac-par[1]^par[2]*(time1^par[2]-bigT^par[2])) - s1
        t2 <- exp(-(lambdac*bigT)^gammac-par[1]^par[2]*(time2^par[2]-bigT^par[2])) - s2
        terror <- t1^2+t2^2
        return(terror)
      }
      
      output <- optim(par = c(0, 1), fn = estWeib)
      lambdat <- output$par[1]
      gammat <- output$par[2]
    
    }
    

    #Simulating the control and treatment data
    
    output <- simulateDTEWeibullData(n, n, gammat, gammac, lambdat, lambdac, bigT, recTime, 0.8)
    
    combinedData <- output$dataCombined
    
    censvec[i] <- output$censTime
    
    #Performing a log-rank test on the combined data set
    test <- logrank.test(combinedData$survival_time, combinedData$status, combinedData$group, rho = 0, gamma = 1)
    
    #Making sure that the HR is less than 1
    coxmodel <- coxph(Surv(survival_time, status)~group, data = combinedData)
    deltad <- as.numeric(exp(coef(coxmodel)))
    
    #Include HR < 1 here
    vec[i] <- (test$test$Chisq > qchisq(0.95, 1) & deltad<1)
  }
  return(list(ass=mean(vec), cens = mean(censvec)))
}


calcAssPowerFunc <- function(type){
  nvec <- ceiling(seq(20, 500, by=10))
  output <- sapply(X = nvec, FUN = powerassFunc, type = type)
  smoothedout <- loess(unlist(output[1,])~nvec)
  return(list(smoothedout = smoothedout, nvec = nvec, cens = output$cens))
  
}

ass <- calcAssPowerFunc("assurance")
power <- calcAssPowerFunc("power")
powerND <- calcAssPowerFunc("powerND")
flexass <- calcAssPowerFunc("flexAss")

#png("PowerAss.png", units="in", width=8, height=5, res=700)

plot(power$nvec*2, predict(power$smoothedout), ylim=c(0,1), type="l", xlab = "Total sample size", ylab = "Power/assurance", lty = 2)
lines(ass$nvec*2, predict(ass$smoothedout), col="red", lty=1)
lines(powerND$nvec*2, predict(powerND$smoothedout), col="blue", lty=3)
lines(flexass$nvec*2, predict(flexass$smoothedout), col = "lightgreen", lty = 4)

legend("bottomright", legend = c("Assurance", "Power", "Power assuming no delay", "Flexible assurance (Section 5)"), col=c("red", "black", "blue", "lightgreen"), lty=1:4)

#dev.off()






