
library(survival)
#First data set
DTEDataSetT <- read.csv(file = "DTE/Datasets/IPD_Treatment.csv")

DTEDataSetC <- read.csv(file = "DTE/Datasets/IPD_Control.csv")



TreatmentData <- data.frame(time = DTEDataSetT$Survival.time, cens = DTEDataSetT$Status)

ControlData <- data.frame(time = DTEDataSetC$Survival.time, cens = DTEDataSetC$Status)

TreatmentFit <- survfit(Surv(time, cens)~1, data = TreatmentData)

ControlFit <- survfit(Surv(time, cens)~1, data = ControlData)

plot(TreatmentFit, conf.int = F, col="red")

lines(ControlFit, conf.int = F, col="blue")

##################

weibfit <- survreg(Surv(time, cens)~1, data = ControlData, dist = "weibull")

gamma2 <- as.numeric(exp(-weibfit$icoef[2]))

lambda2 <- as.numeric(1/(exp(weibfit$icoef[1])))

controltime <- seq(0, 36, by=0.01)
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

treatmenttime <- seq(0, 6, by=0.01)
treatmentsurv <- controlsurv <- exp(-(lambda2*treatmenttime)^gamma2)

lines(treatmenttime, treatmentsurv, col="green")


treatmenttime1 <- seq(6, 36, by=0.01)

optimfunc <- function(par){
  diff <- 0
  for (i in 1:length(treatmenttime1)){
    y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
    treatmentsurv <- exp(-(lambda2*6)^gamma2 - par[1]^par[2]*(treatmenttime1[i]^par[2]-6^par[2]))
    diff <- diff + (y-treatmentsurv)^2 
  }
  return(diff) 
}



optimoutput <-  optim(par = c(0.2, 1), fn = optimfunc)
lambda1 <- optimoutput$par[1]
gamma1 <- optimoutput$par[2]


treatmenttime1 <- seq(6, 36, by=0.01)
treatmentsurv1 <- exp(-(lambda2*6)^gamma2 - lambda1^gamma1*(treatmenttime1^gamma1-6^gamma1))

lines(treatmenttime1, treatmentsurv1, col="red")

###Now we do it but we do not allow gamma1 to vary

treatmenttime <- seq(0, 6, by=0.01)
treatmentsurv <- controlsurv <- exp(-(lambda2*treatmenttime)^gamma2)

lines(treatmenttime, treatmentsurv, col="green")


treatmenttime1 <- seq(6, 36, by=0.01)

optimfunc <- function(par){
  diff <- 0
  gamma1 <- gamma2
  for (i in 1:length(treatmenttime1)){
    y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
    treatmentsurv <- exp(-(lambda2*6)^gamma2 - par[1]^gamma1*(treatmenttime1[i]^gamma1-6^gamma1))
    diff <- diff + (y-treatmentsurv)^2 
  }
  return(diff) 
}



lambda1only <-  optimize(optimfunc, c(0, 2))$minimum


treatmenttime1 <- seq(6, 36, by=0.01)
treatmentsurv1 <- exp(-(lambda2*6)^gamma2 - lambda1only^gamma2*(treatmenttime1^gamma2-6^gamma2))

lines(treatmenttime1, treatmentsurv1, col="red", lty=2)

legend("topright", legend = c("Control", "Treatment (best fit)", "Treatment (fixed)"), col=c("blue", "red", "red"), lty=c(1, 1, 2))