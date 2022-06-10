##Plots for my confirmation review report

##A plot to show the power/assurance of the moxonidine trial

#png("AssMoxonidine.png", width = 600, height=400)


N1vec <- floor(seq(50, 500, length=50))
N2vec <- floor(seq(50,500, length = 50))
powervec <- rep(NA, length = length(N1vec))
for (i in 1:length(powervec)){
  power <- power.prop.test(n=N1vec[i],p1=0.45,p2=0.3)
  powervec[i] <- power$power
}


plot(N1vec, powervec, type = "l", ylim = c(0,1), xlab = "Total sample size", xaxt = "n",
     ylab = "Power/Assurance", lty=1)
axis(1, at = seq(0, 500,by = 100), labels = seq(0, 1000, by = 200))

assurancefunc <- function(n1, n2, m, v){
  n <- 40e1
  zvec <- rep(NA, length=n)
  for (i in 1:n){
    theta1 <- rbeta(1, 10.7, 13.1) #Control
    rho <- truncnorm::rtruncnorm(1, mean = m, sd = sqrt(v), a = (theta1-1), b = theta1)
    theta2 <- theta1 - rho #Treatment
    control <- rbinom(1, n1, theta1)
    treatment <- rbinom(1, n2, theta2)
    MoxonidineData <- data.frame(RaisedcTni = c(control, treatment),
                                 NotRaisedcTni = c(n1 - control, n2 - treatment),
                                 row.names = c("Control", "Treatment"))
    test <- chisq.test(MoxonidineData, correct = F)
    zvec[i] <- test$p.value < 0.05
  }
  mean(zvec)
}

#Scenario 1
ass1 <- rep(NA, length=length(N1vec))
for (i in 1:length(ass1)){
  ass1[i] <- assurancefunc(N1vec[i], N2vec[i], m = 0.15, v = 0.0001)
}
lo1 <- loess(ass1~N1vec)
lines(N1vec, predict(lo1), lty=2, col="blue")

#Scenario 2
ass2 <- rep(NA, length=length(N1vec))
for (i in 1:length(ass2)){
  ass2[i] <- assurancefunc(N1vec[i], N2vec[i], m = 0.15, v = 0.01)
}
lo2 <- loess(ass2~N1vec)
lines(N1vec, predict(lo2), lty=3, col="red")


#Scenario 3
ass3 <- rep(NA, length=length(N1vec))
for (i in 1:length(ass3)){
  ass3[i] <- assurancefunc(N1vec[i], N2vec[i], m = 0.1, v  = 0.01)
}

lo3 <- loess(ass3~N1vec)
lines(N1vec, predict(lo3), lty=4, col="green")


legend("bottomright", legend = c("Power", "Scenario 1", "Scenario 2", "Scenario 3"),
       col=c("black", "blue", "red", "green"), lty=1:4)

#dev.off()

#A plot to show DTE

#png("DTE.png", width = 600, height=400)

bigT <- 3

lambda1 <- 0.04
gamma1 <- 1.2

lambda2 <- 0.1
gamma2 <- 0.8


controltime <- seq(0, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)
plot(controltime, controlsurv, type="l", col="blue", ylab="Survival", xlab = "Time")

sharedtime <- seq(0, bigT, by=0.01)
sharedsurv <- exp(-(lambda2*sharedtime)^gamma2)

lines(sharedtime, sharedsurv, col="green")

#Plotting the median treatment lines

#This situation is where gamma1 is allowed to vary (we are using this as a control - for example a clinician might think this is appropriate
#so we are seeing what we are losing by enforcing gamma1 = gamma2)

treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))

lines(treatmenttime, treatmentsurv, col="red")


legend("topright", legend = c("Delay", "Control", "Treatment"), col=c("green", "blue", "red"), lty=1)

axis(1, at = 3, labels = "T")
abline(v = 3, lty=2)


#dev.off()

#The two plots showing DTE from the data sets

#First plot with just the K-M 

library(survival)
#First data set
DTEDataSet1 <- read.csv(file = "DTE/Datasets/Satrajits/CM017PFS.csv")

#png("DTEPlot1.png", width = 600, height=400)

TreatmentData <- data.frame(time = DTEDataSet1$NIVO_Time[1:135], cens = DTEDataSet1$NIVO_Event_Flag[1:135])

ControlData <- data.frame(time = DTEDataSet1$DOX_Time, cens = DTEDataSet1$DOX_Event_Flag)

TreatmentFit <- survfit(Surv(time, cens)~1, data = TreatmentData)

ControlFit <- survfit(Surv(time, cens)~1, data = ControlData)

plot(TreatmentFit, conf.int = F, col="red")

lines(ControlFit, conf.int = F, col="blue")

legend("topright", legend = c("Control", "Treatment"), col=c("blue", "red"), lty=1)

#dev.off()

#Second plot with the K-M and the simplification (or not) overlaid

#png("DTEPlot2.png", width = 600, height=400)

TreatmentData <- data.frame(time = DTEDataSet1$NIVO_Time[1:135], cens = DTEDataSet1$NIVO_Event_Flag[1:135])

ControlData <- data.frame(time = DTEDataSet1$DOX_Time, cens = DTEDataSet1$DOX_Event_Flag)

TreatmentFit <- survfit(Surv(time, cens)~1, data = TreatmentData)

ControlFit <- survfit(Surv(time, cens)~1, data = ControlData)

plot(TreatmentFit, conf.int = F, col="red")

lines(ControlFit, conf.int = F, col="blue")

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

legend("topright", legend = c("Delay", "Control", "Treatment (best fit)", "Treatment (fixed)"), col=c("green", "blue", "red", "red"), lty=c(1, 1, 2))

#dev.off()



