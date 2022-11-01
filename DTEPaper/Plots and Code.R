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
  } else if (author=="Yen"){
    trialLength <- 30
    THat <- 3.5
  } else if (author=="Borghaei"){
    trialLength <- 60
    THat <- 6
  }
  
  controldata <- read.csv(file = paste0("DTEPaper/data/", author, "/IPD-control.csv"))
  treatmentdata <- read.csv(file = paste0("DTEPaper/data/", author, "/IPD-treatment.csv"))
  
  
  combinedData <- data.frame(time = c(controldata$Survival.time, treatmentdata$Survival.time), 
                             status = c(controldata$Status, treatmentdata$Status), 
                             group = c(rep("Control", nrow(controldata)), rep("Treatment", nrow(treatmentdata))))
  
  kmfit <- survfit(Surv(time, status)~group, data = combinedData)
  plot(kmfit, conf.int = F, col=c("blue", "red"), xlab = "Time (months)", ylab = "Progression free survival (% of patients)")
  
  
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
    powervec <- rep(NA, 500)
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
  
  nvec <- round(seq(10, 500, length = 20))
  powervecs1 <- mapply(powerfunc, nvec, s1lambdat, s1gammat) 
  powervecs2 <- mapply(powerfunc, nvec, s2lambdat, s2gammat) 
  
  powers1smooth <- loess(powervecs1~nvec)
  powers2smooth <- loess(powervecs2~nvec)
  
  plot(nvec*2, predict(powers1smooth), type="l", col="blue", ylim=c(0, 1), xlab = "Total sample size", ylab = "Power")
  lines(nvec*2, predict(powers2smooth), col="red", lty=2)
  
  legend("bottomright", legend = c("Scenario 1", "Scenario 2"), col = c("blue", "red"), lty=1:2)
}

DTEDataSetsFunc("Brahmer")
DTEDataSetsFunc("Yen")
DTEDataSetsFunc("Borghaei")



n1 <- 1092
n2 <- 1092
Nsim <- 10e4
testvec <- rep(NA, Nsim)
for (i in 1:Nsim){
  theta1 <- rbeta(1, 47, 140)
  rho <- rtruncnorm(1,a = (theta1-1), b = theta1, mean = 0.05, sd = sqrt(0.001))
  theta2 <- theta1 - rho
  control <- rbinom(1,n1, prob=theta1)
  treatment <- rbinom(1, n2, prob=theta2)
  finalData <- data.frame(HF = c(control, treatment),
                          NotHF = c(n1 - control, n2 - treatment),
                          row.names = c("Control", "Treatment"))
  test <- chisq.test(finalData, correct = F)
  testvec[i] <- test$p.value<0.05
}

mean(testvec)

#This produces the power/assurance curve found in Section 2.1

#png("MoxPowerAss.png", units="in", width=5, height=5, res=700)
#Calculates the power
powerFunc <- function(n1, n2){
  powervec <- rep(NA, 200)
  for (i in 1:200){
    theta1 <- 0.25
    theta2 <- 0.2
    control <- rbinom(1,n1, prob=theta1)
    treatment <- rbinom(1, n2, prob=theta2)
    finalData <- data.frame(HF = c(control, treatment),
                            NotHF = c(n1 - control, n2 - treatment),
                            row.names = c("Control", "Treatment"))
    test <- chisq.test(finalData, correct = F)
    powervec[i] <- test$p.value<0.05
  }
  mean(powervec)
}

svec <- seq(30, 2000, by=50)
powervec <- rep(NA, length(svec))
for (j in 1:length(svec)){
  powervec[j] <- powerFunc(svec[j], svec[j])
}

powersmooth <- loess(powervec~svec)

plot(svec*2, predict(powersmooth), ylab="Power/Assurance", xlab="Total sample size", type="l", ylim=c(0,1), col="blue")

#Calculates the assurance
assFunc <- function(n1, n2){
  assvec <- rep(NA, 500)
  for (i in 1:500){
    theta1 <- rbeta(1, 47, 140)
    rho <- rtruncnorm(1,a = (theta1-1), b = theta1, mean = 0.05, sd = sqrt(0.001))
    theta2 <- theta1 - rho
    control <- rbinom(1,n1, prob=theta1)
    treatment <- rbinom(1, n2, prob=theta2)
    finalData <- data.frame(HF = c(control, treatment),
                            NotHF = c(n1 - control, n2 - treatment),
                            row.names = c("Control", "Treatment"))
    test <- chisq.test(finalData, correct = F)
    assvec[i] <- test$p.value<0.05
  }
  mean(assvec)
}

assvec <- rep(NA, length(svec))
for (j in 1:length(svec)){
  assvec[j] <- assFunc(svec[j], svec[j])
}

asssmooth <- loess(assvec~svec)

lines(svec*2, predict(asssmooth), col="red", lty=2)

legend("bottomright", legend = c("Power", "Assurance"), col=c("blue", "red"), lty=1:2)

#This produces the figure found in Section 2.2
lambdac <- 0.08
gammac <- 0.8
HR <- 0.5
lambdat <- exp((log(HR)/gammac)+log(lambdac))
gammat <- gammac
bigT <- 4


controltime <- seq(0, exp((1.527/gammac)-log(lambdac))*1.1, by=0.01)
controlcurve <- exp(-(lambdac*controltime)^gammac)
treatmenttime1 <- seq(0, bigT, by=0.01)
treatmentsurv1 <- exp(-(lambdac*treatmenttime1)^gammac)
treatmenttime2 <- seq(bigT, exp((1.527/gammac)-log(lambdac))*1.1, by=0.01)
treatmentsurv2 <- exp(-(lambdac*bigT)^gammac - lambdat^gammat*(treatmenttime2^gammat-bigT^gammat))



#png("DTE.png", units="in", width=5, height=5, res=700)

plot(controltime, controlcurve, type="l", col="blue", xlab="Time", ylab="Survival")
lines(treatmenttime1, treatmentsurv1, col="green")
lines(treatmenttime2, treatmentsurv2, col="red")
legend("topright", legend = c("Delay", "Control", "Treatment"), col=c("green", "blue", "red"), lty=1)
abline(v = 4, lty=2)
axis(side = 1, at = 4, labels = "T")

#dev.off()

#The following code produces the plots in Section 3.3

#Data sets that Satrajit sent me for DTE
library(survival)
#First data set

#png("DS2.png", units="in", width=5, height=5, res=700)
DTEDataSet1 <- read.csv(file = "DTE/Datasets/Satrajits/CM017PFS.csv")

TreatmentData <- data.frame(time = DTEDataSet1$NIVO_Time[1:135], cens = DTEDataSet1$NIVO_Event_Flag[1:135])

ControlData <- data.frame(time = DTEDataSet1$DOX_Time, cens = DTEDataSet1$DOX_Event_Flag)

TreatmentFit <- survfit(Surv(time, cens)~1, data = TreatmentData)

ControlFit <- survfit(Surv(time, cens)~1, data = ControlData)

plot(TreatmentFit, conf.int = F, col="red", xlim=c(0,18), xlab="Time (months)", ylab="Survival")

lines(ControlFit, conf.int = F, col="blue")

legend("topright", legend = c("Control", "Treatment"), col=c("blue", "red"), lty=1)


##Fitting Weibull model to the data

weibfit <- survreg(Surv(time, cens)~1, data = ControlData, dist = "weibull")

gammac <- as.numeric(exp(-weibfit$icoef[2]))

lambdac <- as.numeric(1/(exp(weibfit$icoef[1])))

controltime <- seq(0, 15, by=0.01)
controlsurv <- exp(-(lambdac*controltime)^gammac)


lines(controltime, controlsurv, col="blue")

#Two scenarios, 1 where we allow gammat to vary, 1 where we do not


#First scenario

#We estimate T to be about 3 here

#Need to find a least squares estimate for this Weibull parameterisation


#gammac = 1.07
#lambdac = 0.21

#Write a function which finds the squared differences

kmfit <- survfit(Surv(time, cens)~1, data = TreatmentData)

treatmenttime <- seq(0, 3, by=0.01)
treatmentsurv <- controlsurv <- exp(-(lambdac*treatmenttime)^gammac)

lines(treatmenttime, treatmentsurv, col="green")


treatmenttime1 <- seq(3, 15, by=0.01)

optimfunc <- function(par){
  diff <- 0
  for (i in 1:length(treatmenttime1)){
    y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
    treatmentsurv <- exp(-(lambdac*3)^gammac - par[1]^par[2]*(treatmenttime1[i]^par[2]-3^par[2]))
    diff <- diff + (y-treatmentsurv)^2 
  }
  return(diff) 
}



optimoutput <-  optim(par = c(0.2, 1), fn = optimfunc)
lambdat <- optimoutput$par[1]
gammat <- optimoutput$par[2]


treatmenttime1 <- seq(3, 15, by=0.01)
treatmentsurv1 <- exp(-(lambdac*3)^gammac - lambdat^gammat*(treatmenttime1^gammat-3^gammat))

lines(treatmenttime1, treatmentsurv1, col="red")

###Now we do it but we do not allow gammat to vary

treatmenttime <- seq(0, 3, by=0.01)
treatmentsurv <- controlsurv <- exp(-(lambdac*treatmenttime)^gammac)

lines(treatmenttime, treatmentsurv, col="green")


treatmenttime1 <- seq(3, 15, by=0.01)

optimfunc <- function(par){
  diff <- 0
  gammat <- gammac
  for (i in 1:length(treatmenttime1)){
    y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
    treatmentsurv <- exp(-(lambdac*3)^gammac - par[1]^gammat*(treatmenttime1[i]^gammat-3^gammat))
    diff <- diff + (y-treatmentsurv)^2 
  }
  return(diff) 
}



lambdatonly <-  optimize(optimfunc, c(0, 2))$minimum


treatmenttime1 <- seq(3, 15, by=0.01)
treatmentsurv1 <- exp(-(lambdac*3)^gammac - lambdatonly^gammac*(treatmenttime1^gammac-3^gammac))

lines(treatmenttime1, treatmentsurv1, col="red", lty=2)

legend("topright", legend = c("Delay", "Control", "Treatment (best fit)", "Treatment (fixed)"), col=c("green", "blue", "red", "red"), lty=c(1, 1, 2))


#This code produces the KM and overlaid lines in Section 3.3
#png("DS2.png", units="in", width=5, height=5, res=700)
DTEDataSet1 <- read.csv(file = "DTE/Datasets/Satrajits/CM017PFS.csv")

TreatmentData <- data.frame(time = DTEDataSet1$NIVO_Time[1:135], cens = DTEDataSet1$NIVO_Event_Flag[1:135])

ControlData <- data.frame(time = DTEDataSet1$DOX_Time, cens = DTEDataSet1$DOX_Event_Flag)

TreatmentFit <- survfit(Surv(time, cens)~1, data = TreatmentData)

ControlFit <- survfit(Surv(time, cens)~1, data = ControlData)

plot(TreatmentFit, conf.int = F, col="red", xlim=c(0,18), xlab="Time (months)", ylab="Survival")

lines(ControlFit, conf.int = F, col="blue")


legend("topright", legend = c("Control", "Treatment"), col=c("blue", "red"), lty=1)




##Fitting Weibull model to the data

weibfit <- survreg(Surv(time, cens)~1, data = ControlData, dist = "weibull")

gammac <- as.numeric(exp(-weibfit$icoef[2]))

lambdac <- as.numeric(1/(exp(weibfit$icoef[1])))

controltime <- seq(0, 15, by=0.01)
controlsurv <- exp(-(lambdac*controltime)^gammac)


lines(controltime, controlsurv, col="blue")

#Two scenarios, 1 where we allow gammat to vary, 1 where we do not


#First scenario

#We estimate T to be about 3 here

#Need to find a least squares estimate for this Weibull parameterisation


#gammac = 1.07
#lambdac = 0.21

#Write a function which finds the squared differences

kmfit <- survfit(Surv(time, cens)~1, data = TreatmentData)

treatmenttime <- seq(0, 3, by=0.01)
treatmentsurv <- controlsurv <- exp(-(lambdac*treatmenttime)^gammac)

lines(treatmenttime, treatmentsurv, col="green")


treatmenttime1 <- seq(3, 15, by=0.01)

optimfunc <- function(par){
  diff <- 0
  for (i in 1:length(treatmenttime1)){
    y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
    treatmentsurv <- exp(-(lambdac*3)^gammac - par[1]^par[2]*(treatmenttime1[i]^par[2]-3^par[2]))
    diff <- diff + (y-treatmentsurv)^2 
  }
  return(diff) 
}



optimoutput <-  optim(par = c(0.2, 1), fn = optimfunc)
lambdat <- optimoutput$par[1]
gammat <- optimoutput$par[2]


treatmenttime1 <- seq(3, 15, by=0.01)
treatmentsurv1 <- exp(-(lambdac*3)^gammac - lambdat^gammat*(treatmenttime1^gammat-3^gammat))

lines(treatmenttime1, treatmentsurv1, col="red")

###Now we do it but we do not allow gammat to vary

treatmenttime <- seq(0, 3, by=0.01)
treatmentsurv <- controlsurv <- exp(-(lambdac*treatmenttime)^gammac)

lines(treatmenttime, treatmentsurv, col="green")


treatmenttime1 <- seq(3, 15, by=0.01)

optimfunc <- function(par){
  diff <- 0
  gammat <- gammac
  for (i in 1:length(treatmenttime1)){
    y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
    treatmentsurv <- exp(-(lambdac*3)^gammac - par[1]^gammat*(treatmenttime1[i]^gammat-3^gammat))
    diff <- diff + (y-treatmentsurv)^2 
  }
  return(diff) 
}



lambdatonly <-  optimize(optimfunc, c(0, 2))$minimum


treatmenttime1 <- seq(3, 15, by=0.01)
treatmentsurv1 <- exp(-(lambdac*3)^gammac - lambdatonly^gammac*(treatmenttime1^gammac-3^gammac))

lines(treatmenttime1, treatmentsurv1, col="red", lty=2)

legend("topright", legend = c("Delay", "Control", "Treatment (best fit)", "Treatment (fixed)"), col=c("green", "blue", "red", "red"), lty=c(1, 1, 2))



#This code produces the power plot in Section 3.3


#Need to calculate power in the scenario when we allow gammat to vary

#So we need to simulate some data according to these parameters


source("DTEPaper/functions.R")

#png("PowerVaryingFixed.png", units="in", width=8, height=5, res=700)
powerFunc <- function(bigT, lambdac, gammac, lambdat, gammat, n1, n2){
  powervec <- rep(NA, 100)
  for (i in 1:100){
    z <- simulateDTEWeibullData(bigT, lambdac, gammac, lambdat, gammat, n1, n2)
    combinedData <- data.frame(time = c(z$controldata, z$treatmentdata), group = c(rep("control", n1), rep("treatment", n2)))
    test <- survdiff(Surv(time)~group, data = combinedData)
    powervec[i] <- test$chisq > qchisq(0.95, 1)
  }
  return(mean(powervec))
}

n1vec <- seq(20, 400, by=20)
n2vec <- seq(20, 400, by=20)
nvec <- n1vec+n2vec
powervary <- rep(NA, length(n1vec))
for (j in 1:length(n1vec)){
  powervary[j] <- powerFunc(3, lambdac, gammac, lambdat, gammat, n1vec[j], n2vec[j])
}

powervarysmooth <- loess(powervary~nvec)

plot(nvec, predict(powervarysmooth), type = "l", ylim=c(0,1), col="red", lty=1, xlab="Total sample size", ylab="Power")

powerfixed <- rep(NA, length(n1vec))
for (j in 1:length(n1vec)){
  powerfixed[j] <- powerFunc(3, lambdac, gammac, lambdatonly, gammac, n1vec[j], n2vec[j])
}

powerfixedsmooth <- loess(powerfixed~nvec)

lines(nvec, predict(powerfixedsmooth), col="blue", lty=2)

legend("bottomright", legend = c(expression(paste("Varying ", gamma[1])), expression(paste("Fixed ", gamma[1]))), lty=1:2, col=c("red", "blue"))

#dev.off()


#This code produces the plot seen in Section 4.1.3
#########################################################
##We have done here, but just with the assurance
##We can look at the power (in both cases) on Monday
########################################################

png("PowerAss.png", units="in", width=8, height=5, res=700)
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

dev.off()








