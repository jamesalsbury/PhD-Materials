#Here is the R script used to produce the figures/results from the DTE paper

#This finds the assurance for the example in Section 2.1
library(truncnorm)
library(survival)
n1 <- 1092
n2 <- 1092
Nsim <- 10e4
testvec <- rep(NA, Nsim)
for (i in 1:Nsim){
  theta1 <- rbeta(1, 47, 140)
  rho <- rtruncnorm(1,a = (theta1-1), b = theta1, mean = 0.05, sd = sqrt(0.005))
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

#This produces the figure found in Section 2.2
lambda2 <- 0.08
gamma2 <- 0.8
HR <- 0.5
lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
gamma1 <- gamma2
bigT <- 4


controltime <- seq(0, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
controlcurve <- exp(-(lambda2*controltime)^gamma2)
treatmenttime1 <- seq(0, bigT, by=0.01)
treatmentsurv1 <- exp(-(lambda2*treatmenttime1)^gamma2)
treatmenttime2 <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
treatmentsurv2 <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime2^gamma1-bigT^gamma1))



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

#This code produces the power plot in Section 3.3



#Need to calculate power in the scenario when we allow gamma1 to vary

#So we need to simulate some data according to these parameters

source("functions.R")

png("PowerVaryingFixed.png", units="in", width=8, height=5, res=700)
powerFunc <- function(bigT, lambda2, gamma2, lambda1, gamma1, n1, n2){
  powervec <- rep(NA, 100)
  for (i in 1:100){
    z <- simulateDTEWeibullData(bigT, lambda2, gamma2, x1, y1, n1, n2)
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
  powervary[j] <- powerFunc(3, lambda2, gamma2, x1, y1, n1vec[j], n2vec[j])
}

powervarysmooth <- loess(powervary~nvec)

plot(nvec, predict(powervarysmooth), type = "l", ylim=c(0,1), col="red", lty=1, xlab="Total sample size", ylab="Power")

powerfixed <- rep(NA, length(n1vec))
for (j in 1:length(n1vec)){
  powerfixed[j] <- powerFunc(3, lambda2, gamma2, z1, gamma2, n1vec[j], n2vec[j])
}

powerfixedsmooth <- loess(powerfixed~nvec)

lines(nvec, predict(powerfixedsmooth), col="blue", lty=2)

legend("bottomright", legend = c(expression(paste("Varying ", gamma[1])), expression(paste("Fixed ", gamma[1]))), lty=1:2, col=c("red", "blue"))

dev.off()


#This code produces the plot seen in Section 4.1.3
png("PowerAss.png", units="in", width=8, height=5, res=700)
#Calculating the power first
powerFunc <- function(nc, nt, N, Lmax, lambda2, gamma2){
  #Making the simplification
  gamma1 <- gamma2
  #The parameters are fixed
  bigT <- 6
  HR <- 0.6
  #Can use the HR to calculate lambda1
  lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
  #Setting up the assurance vector
  powervec <- rep(NA, N)
  for (i in 1:N){
    #Generating the control data
    controldata <- data.frame(time = rweibull(nc, gamma2, 1/lambda2))
    
    #Generating the treatment data
    CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
    u <- runif(nt)
    
    suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
    
    #Combining the control and treatment data
    DataCombined <- data.frame(time = c(controldata$time, z),
                               group = c(rep("Control", nc), rep("Treatment", nt)))
    
    #Censoring the observations
    DataCombined$cens <- DataCombined$time < Lmax
    
    DataCombined$cens <- DataCombined$cens*1
    
    
    if (sum(DataCombined$cens)==(nc+nt)){
      
    } else {
      DataCombined[DataCombined$cens==0,]$time <- Lmax
    }
    #Performing a log-rank test on the combined data set
    test <- survdiff(Surv(time, cens)~group, data = DataCombined)
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
  finalpowervec[j] <- powerFunc(ncvec[j], ntvec[j], 500, 40, 0.08, 0.8)
}

samplesizevec <- ncvec+ntvec

#Smoothing the output
powersmooth <- loess(finalpowervec~samplesizevec)

#Plotting the output
plot(samplesizevec, predict(powersmooth), ylim=c(0,1), 
     type = "l", xlab="Total sample size", ylab="Power/assurance", col="blue")



#Now we compute the assurance
assFunc <- function(nc, nt, N, Lmax, lambda2, gamma2){
  #Making the simplification
  gamma1 <- gamma2
  #Setting up the assurance vector
  assvec <- rep(NA, N)
  for (i in 1:N){
    #Generating the control data
    controldata <- data.frame(time = rweibull(nc, gamma2, 1/lambda2))
    #Sampling from the prior distributions for T and HR
    bigT <- rtruncnorm(1, mean = 6, sd = 2.22, a = 0)
    HR <- rbeta(1, 10.8, 6.87)
    #Finding lambda1 from the HR
    lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
    #Generating the treatment data
    CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
    u <- runif(nt)
    
    suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
    
    #Combining the control and treatment data
    DataCombined <- data.frame(time = c(controldata$time, z),
                               group = c(rep("Control", nc), rep("Treatment", nt)))
    
    #Censoring the observations
    DataCombined$cens <- DataCombined$time < Lmax
    
    DataCombined$cens <- DataCombined$cens*1
    
    
    if (sum(DataCombined$cens)==(nc+nt)){
      
    } else {
      DataCombined[DataCombined$cens==0,]$time <- Lmax
    }
    #Performing a log-rank test on the combined data set
    test <- survdiff(Surv(time, cens)~group, data = DataCombined)
    assvec[i] <- test$chisq > qchisq(0.95, 1)
  }
  #The mean of the assurance vector is the estimated assurance for the given parameters
  return(mean(assvec))
}


#Setting up the assurance vector
finalassvec <- rep(NA, length(ncvec))
#Finding the assurance for varying sample sizes
for (j in 1:length(finalassvec)){
  finalassvec[j] <- assFunc(ncvec[j], ntvec[j], 500, 40, 0.08, 0.8)
}

#Smoothing the output
asssmooth <- loess(finalassvec~samplesizevec)

#Plotting the output
lines(samplesizevec, predict(asssmooth), lty=2, col="red")

#We will also calculate power in the situation where we have no delay
powerFuncNoDelay <- function(nc, nt, N, Lmax, lambda2, gamma2){
  #Making the simplification
  gamma1 <- gamma2
  #The parameters are fixed
  bigT <- 0
  HR <- 0.6
  #Can use the HR to calculate lambda1
  lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
  #Setting up the assurance vector
  powernodelayvec <- rep(NA, N)
  for (i in 1:N){
    #Generating the control data
    controldata <- data.frame(time = rweibull(nc, gamma2, 1/lambda2))
    
    #Generating the treatment data
    CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
    u <- runif(nt)
    
    suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
    
    #Combining the control and treatment data
    DataCombined <- data.frame(time = c(controldata$time, z),
                               group = c(rep("Control", nc), rep("Treatment", nt)))
    
    #Censoring the observations
    DataCombined$cens <- DataCombined$time < Lmax
    
    DataCombined$cens <- DataCombined$cens*1
    
    
    if (sum(DataCombined$cens)==(nc+nt)){
      
    } else {
      DataCombined[DataCombined$cens==0,]$time <- Lmax
    }
    #Performing a log-rank test on the combined data set
    test <- survdiff(Surv(time, cens)~group, data = DataCombined)
    powernodelayvec[i] <- test$chisq > qchisq(0.95, 1)
  }
  #The mean of the assurance vector is the estimated assurance for the given parameters
  return(mean(powernodelayvec))
}

#Setting up the assurance vector
finalpowenodelayrvec <- rep(NA, length(ncvec))
#Finding the assurance for varying sample sizes
for (j in 1:length(finalpowenodelayrvec)){
  finalpowenodelayrvec[j] <- powerFuncNoDelay(ncvec[j], ntvec[j], 500, 40, 0.08, 0.8)
}


#Smoothing the output
powernodelaysmooth <- loess(finalpowenodelayrvec~samplesizevec)

lines(samplesizevec, predict(powernodelaysmooth), lty=3, col="green")

legend("bottomright", legend =c("Power", "Assurance", "Power assuming no delay"), col=c("blue", "red", "green"), lty=1:3)

for (j in 40:1000){
  if (predict(powernodelaysmooth, newdata = j)>0.8){
    break
  }
}

predict(powernodelaysmooth, newdata = 174)





