
####
#Looking at stopping for futility in the case we have
####

library(survival)
library(tidyverse)
#Start off with traditional power


#We assume the following parameters for our model
lambda2 <- 0.08
gamma2 <- 0.8
lambda1 <- 0.04
gamma1 <- 0.8

#Draw the control curve
controltime <- seq(0, 100, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)
plot(controltime, controlsurv, type="l", col="blue", ylim=c(0,1))

#We assume that the delay is 6 months
bigT <- 6

#Draw the delay curve
delaytime <- seq(0, bigT, by=0.01)
delaysurv <- exp(-(lambda2*delaytime)^gamma2)
lines(delaytime, delaysurv, type="l", col="green")

#Draw the treatment curve
treatmenttime <- seq(bigT, 100 , by=0.01)
treatmentsurv <- exp(-(lambda2*bigT)^gamma2-lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))
lines(treatmenttime, treatmentsurv, type="l", col="red")

#So we have no uncertainty about our parameters here
#Can we construct a power curve for this setup?
#We only want to run the trial for 60 months

trialLength <- 60

powerfunc <- function(n1, n2){
  powervec <- rep(NA, 250)
  eventvec <- rep(NA, 250)
  for (i in 1:length(powervec)){
    controldata <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))
    
    CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
    u <- runif(n2)
    
    suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
    
    DataCombined <- data.frame(time = c(controldata$time, z),
                               group = c(rep("Control", n1), rep("Treatment", n2)))
    
    DataCombined$cens <- DataCombined$time < trialLength
    
    DataCombined$cens <- DataCombined$cens*1
    
    
    
    if (sum(DataCombined$cens)==(n1+n2)){
      
    } else {
      DataCombined[DataCombined$cens==0,]$time <- trialLength
    }
    
    
    eventvec <- sum(DataCombined$cens==1)
    
    test <- survdiff(Surv(time, cens)~group, data = DataCombined)
    powervec[i] <- test$chisq > qchisq(0.95, 1)
  }
  return(list(powervec = mean(powervec), eventvec=mean(na.omit(eventvec))))
}


n1vec <- seq(10, 500, by=20)
n2vec <- seq(10, 500, by=20)
powervec <- rep(NA, length(n1vec))
eventvec <- rep(NA, length(n1vec))
for (j in 1:length(n1vec)){
  powervec[j] <- powerfunc(n1vec[j], n2vec[j])$powervec
  eventvec[j] <- powerfunc(n1vec[j], n2vec[j])$eventvec
}

nvec <- n1vec+n2vec
powersmooth <- loess(powervec~nvec)
par(mar = c(10, 10, 10, 10))
plot(nvec, predict(powersmooth), ylim=c(0,1), ylab = "Power", xlab="Total sample size", type="l")

eventsmooth <- loess(eventvec~nvec)

par(new = TRUE)		

plot(nvec, predict(eventsmooth), ylab = "", xlab="", type="l", axes=F)
axis(side = 4, at = seq(0, 1000, by=100))
mtext("Events", side = 4, line = 3)


#How would we stop for futility in this scenario?
#We could take the results so far, calculate the hazard ratio post-treatment and simulate data according to this hazard ratio
#We need to have already gone past the delay
#We then calculate conditional power according to these simulated observations
#If CP lower than 0.1, say, we stop for futility
#We base on events rather than sample size

#Lets assume we have an interim look at 10 months

#For 80% power, we need 320 patients
#We expect to see 300 events
#So if we do it event driven

#We need to simulate clinical trials with the following differences: time of delay and post-delay HR
#We look at the data say 10 months in


#We have planned for a HR of 0.574

#We only consider uncertainty in length of delay now

#Length of delay is 8 months, when is the best time to conduct a futility look?
#But why are we looking at futility, when the truth is that it does work?


#Why don't we look at BPP here, as this is better defined than an arbitrary choice for CP

#Before we do this, we need to calculate sample size by using assurance

#We assume the following parameters for our model
lambda2 <- 0.08
gamma2 <- 0.8
gamma1 <- 0.8

#We have elicited the follow distributions
#bigT <- norm(6, 0.5)
#HR <- beta(10, 6)

#Draw the control curve
controltime <- seq(0, 100, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)
plot(controltime, controlsurv, type="l", col="blue", ylim=c(0,1))

#We assume that the delay is 6 months
bigT <- 6

#Draw the delay curve
delaytime <- seq(0, bigT, by=0.01)
delaysurv <- exp(-(lambda2*delaytime)^gamma2)
lines(delaytime, delaysurv, type="l", col="green")

#Draw the treatment curve
#HR mean is 0.625
HR <- 0.625
lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
treatmenttime <- seq(bigT, 100 , by=0.01)
treatmentsurv <- exp(-(lambda2*bigT)^gamma2-lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))
lines(treatmenttime, treatmentsurv, type="l", col="red")

#So we have no uncertainty about our parameters here
#Can we construct a power curve for this setup?
#We only want to run the trial for 60 months

trialLength <- 60

assfunc <- function(n1, n2){
  assvec <- rep(NA, 250)
  eventvec <- rep(NA, 250)
  for (i in 1:length(assvec)){
    bigT <- rnorm(1, 6, 0.5)
    HR <- rbeta(1, 10, 6)
    lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
    
    controldata <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))
    
    CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
    u <- runif(n2)
    
    suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
    
    DataCombined <- data.frame(time = c(controldata$time, z),
                               group = c(rep("Control", n1), rep("Treatment", n2)))
    
    DataCombined$cens <- DataCombined$time < trialLength
    
    DataCombined$cens <- DataCombined$cens*1
    
    
    
    if (sum(DataCombined$cens)==(n1+n2)){
      
    } else {
      DataCombined[DataCombined$cens==0,]$time <- trialLength
    }
    
    
    eventvec <- sum(DataCombined$cens==1)
    
    test <- survdiff(Surv(time, cens)~group, data = DataCombined)
    assvec[i] <- test$chisq > qchisq(0.95, 1)
  }
  return(list(assvec = mean(assvec), eventvec=mean(na.omit(eventvec))))
}


n1vec <- seq(10, 500, by=20)
n2vec <- seq(10, 500, by=20)
assvec <- rep(NA, length(n1vec))
eventvec <- rep(NA, length(n1vec))
for (j in 1:length(n1vec)){
  assvec[j] <- assfunc(n1vec[j], n2vec[j])$assvec
  eventvec[j] <- assfunc(n1vec[j], n2vec[j])$eventvec
}

nvec <- n1vec+n2vec
asssmooth <- loess(assvec~nvec)
par(mar = c(10, 10, 10, 10))
plot(nvec, predict(asssmooth), ylim=c(0,1), ylab = "Assurance", xlab="Total sample size", type="l")

eventsmooth <- loess(eventvec~nvec)

par(new = TRUE)		

plot(nvec, predict(eventsmooth), ylab = "", xlab="", type="l", axes=F)
axis(side = 4, at = seq(0, 1000, by=100))
mtext("Events", side = 4, line = 3)

predict(asssmooth, newdata = 500)

predict(eventsmooth, newdata = 500)

#We require 464 events for 80% power, 500 patients

#Lets simulate some data, we look at 50% through the information fraction, so we look when 232 events have happened
#We simulate the delay to be 8 months, but keep the HR the same
#We can combine the data with the prior to compute BPP

#We simulate 250 patients in each arm

n1 <- 250
n2 <- 250

lambda2 <- 0.08
gamma2 <- 0.8
gamma1 <- 0.8
HR <- 0.625
lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
bigT <- 8

controldata <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))
CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
u <- runif(n2)

suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))

DataCombined <- data.frame(time = c(controldata$time, z),
                           group = c(rep("Control", n1), rep("Treatment", n2)))


DataCombined[order(DataCombined$time),][232,]$time

DataCombined$cens <-  DataCombined$time<DataCombined[order(DataCombined$time),][232,]$time







