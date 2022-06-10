#Data sets that Satrajit sent me for DTE
library(survival)
#First data set
DTEDataSet1 <- read.csv(file = "DTE/Datasets/Satrajits/CM017PFS.csv")

TreatmentData <- data.frame(time = DTEDataSet1$NIVO_Time[1:135], cens = DTEDataSet1$NIVO_Event_Flag[1:135])
 
ControlData <- data.frame(time = DTEDataSet1$DOX_Time, cens = DTEDataSet1$DOX_Event_Flag)

TreatmentFit <- survfit(Surv(time, cens)~1, data = TreatmentData)

ControlFit <- survfit(Surv(time, cens)~1, data = ControlData)

plot(TreatmentFit, conf.int = F, col="red")

lines(ControlFit, conf.int = F, col="blue")

legend("topright", legend = c("Control", "Treatment"), col=c("blue", "red"), lty=1)

##Fitting Weibull model to the data

weibfit <- survreg(Surv(time, cens)~1, data = ControlData, dist = "weibull")

gamma2 <- as.numeric(exp(-weibfit$icoef[2]))

lambda2 <- as.numeric(1/(exp(weibfit$icoef[1])))

controltime <- seq(0, 15, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)


lines(controltime, controlsurv, col="blue")

#Two scenarios, 1 where we allow gamma1 to vary, 1 where we do not


#First scenario

#We estimate T to be about 3 here

#Need to find a least squares estimate for this Weibull parameterisation


#gamma2 = 1.07
#lambda2 = 0.21

#Write a function which finds the squared differences

kmfit <- survfit(Surv(time, cens)~1, data = TreatmentData)

treatmenttime <- seq(0, 3, by=0.01)
treatmentsurv <- controlsurv <- exp(-(lambda2*treatmenttime)^gamma2)

lines(treatmenttime, treatmentsurv, col="green")


treatmenttime1 <- seq(3, 15, by=0.01)

optimfunc <- function(par){
  diff <- 0
  for (i in 1:length(treatmenttime1)){
    y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
    treatmentsurv <- exp(-(lambda2*3)^gamma2 - par[1]^par[2]*(treatmenttime1[i]^par[2]-3^par[2]))
    diff <- diff + (y-treatmentsurv)^2 
  }
  return(diff) 
}



optimoutput <-  optim(par = c(0.2, 1), fn = optimfunc)
lambda1 <- optimoutput$par[1]
gamma1 <- optimoutput$par[2]


treatmenttime1 <- seq(3, 15, by=0.01)
treatmentsurv1 <- exp(-(lambda2*3)^gamma2 - lambda1^gamma1*(treatmenttime1^gamma1-3^gamma1))

lines(treatmenttime1, treatmentsurv1, col="red")

###Now we do it but we do not allow gamma1 to vary

treatmenttime <- seq(0, 3, by=0.01)
treatmentsurv <- controlsurv <- exp(-(lambda2*treatmenttime)^gamma2)

lines(treatmenttime, treatmentsurv, col="green")


treatmenttime1 <- seq(3, 15, by=0.01)

optimfunc <- function(par){
  diff <- 0
  gamma1 <- gamma2
  for (i in 1:length(treatmenttime1)){
    y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
    treatmentsurv <- exp(-(lambda2*3)^gamma2 - par[1]^gamma1*(treatmenttime1[i]^gamma1-3^gamma1))
    diff <- diff + (y-treatmentsurv)^2 
  }
  return(diff) 
}



lambda1only <-  optimize(optimfunc, c(0, 2))$minimum


treatmenttime1 <- seq(3, 15, by=0.01)
treatmentsurv1 <- exp(-(lambda2*3)^gamma2 - lambda1only^gamma2*(treatmenttime1^gamma2-3^gamma2))

lines(treatmenttime1, treatmentsurv1, col="red", lty=2)

legend("topright", legend = c("Control", "Treatment (best fit)", "Treatment (fixed)"), col=c("blue", "red", "red"), lty=c(1, 1, 2))


#Second data set
DTEDataSet2 <- read.csv(file = "DTE/Datasets/Satrajits/CM141OS.csv")


TreatmentData <- data.frame(time = DTEDataSet2$Nivolumab_Time, cens = DTEDataSet2$Nivolumab_Event)

ControlData <- data.frame(time = DTEDataSet2$SOC_Time[1:121], cens = DTEDataSet2$SOC_Event[1:121])

TreatmentFit <- survfit(Surv(time, cens)~1, data = TreatmentData)

ControlFit <- survfit(Surv(time, cens)~1, data = ControlData)

plot(TreatmentFit, conf.int = F, col="red")

lines(ControlFit, conf.int = F, col="blue")

##Fitting Weibull model to the data

weibfit <- survreg(Surv(time, cens)~1, data = ControlData, dist = "weibull")

gamma2 <- as.numeric(exp(-weibfit$icoef[2]))

lambda2 <- as.numeric(1/(exp(weibfit$icoef[1])))

controltime <- seq(0, 15, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)


lines(controltime, controlsurv, col="blue")

#We estimate T to be about 3.5 here


#Need to find a least squares estimate for this Weibull parameterisation


#gamma2 = 1.28
#lambda2 = 0.13

#Write a function which finds the squared differences

kmfit <- survfit(Surv(time, cens)~1, data = TreatmentData)

treatmenttime <- seq(0, 3.5, by=0.01)
treatmentsurv <- controlsurv <- exp(-(lambda2*treatmenttime)^gamma2)

lines(treatmenttime, treatmentsurv, col="green")


treatmenttime1 <- seq(3.5, 15, by=0.01)

optimfunc <- function(par){
  diff <- 0
  for (i in 1:length(treatmenttime1)){
    y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
    treatmentsurv <- exp(-(lambda2*3.5)^gamma2 - par[1]^par[2]*(treatmenttime1[i]^par[2]-3.5^par[2]))
    diff <- diff + (y-treatmentsurv)^2 
  }
  return(diff) 
}



optimoutput <-  optim(par = c(0.2, 1), fn = optimfunc)
lambda1 <- optimoutput$par[1]
gamma1 <- optimoutput$par[2]


treatmenttime1 <- seq(3.5, 15, by=0.01)
treatmentsurv1 <- exp(-(lambda2*3.5)^gamma2 - lambda1^gamma1*(treatmenttime1^gamma1-3.5^gamma1))

lines(treatmenttime1, treatmentsurv1, col="red")

###Now we do it but we do not allow gamma1 to vary


treatmenttime1 <- seq(3.5, 15, by=0.01)

optimfunc <- function(par){
  diff <- 0
  gamma1 <- gamma2
  for (i in 1:length(treatmenttime1)){
    y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
    treatmentsurv <- exp(-(lambda2*3.5)^gamma2 - par[1]^gamma1*(treatmenttime1[i]^gamma1-3.5^gamma1))
    diff <- diff + (y-treatmentsurv)^2 
  }
  return(diff) 
}



lambda1only <-  optimize(optimfunc, c(0, 2))$minimum


treatmenttime1 <- seq(3.5, 15, by=0.01)
treatmentsurv1 <- exp(-(lambda2*3.5)^gamma2 - lambda1only^gamma2*(treatmenttime1^gamma2-3.5^gamma2))

lines(treatmenttime1, treatmentsurv1, col="red", lty=2)

legend("topright", legend = c("Control", "Treatment (best fit)", "Treatment (fixed)"), col=c("blue", "red", "red"), lty=c(1, 1, 2))






