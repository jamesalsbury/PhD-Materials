#Here is the R script used to produce the figures/results from the DTE paper

#This finds the assurance for the example in Section 2.1
library(truncnorm)
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

legend("bottomright", legend =c("Power", "Assurance"), col=c("blue", "red"), lty=1:2)

dev.off()





