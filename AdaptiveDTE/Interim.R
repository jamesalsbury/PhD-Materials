
####
#Looking at stopping for futility in the case we have
####

library(survival)
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

#We assume that the delay is 3 months
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

powerfunc <- function(n1, n2){
  powervec <- rep(NA, 1000)
  for (i in 1:1000){
    controldata <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))
    
    CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
    u <- runif(n2)
    
    suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
    
    DataCombined <- data.frame(time = c(controldata$time, z),
                               group = c(rep("Control", n1), rep("Treatment", n2)))
    
    test <- survdiff(Surv(time)~group, data = DataCombined)
    powervec[i] <- test$chisq > qchisq(0.95, 1)
  }
  return(mean(powervec))
}


n1vec <- seq(10, 500, by=20)
n2vec <- seq(10, 500, by=20)
powervec <- rep(NA, length(n1vec))
for (j in 1:length(n1vec)){
  powervec[j] <- powerfunc(n1vec[j], n2vec[j])
}

nvec <- n1vec+n2vec
powersmooth <- loess(powervec~nvec)

plot(nvec, predict(powersmooth), ylim=c(0,1), ylab = "Power", xlab="Total sample size", type="l")

#How would we stop for futility in this scenario?
#We could take the results so far, calculate the hazard ratio post-treatment and simulate data according to this hazard ratio
#We need to have already gone past the delay
#We then calculate conditional power according to these simulated observations
#If CP lower than 0.1, say, we stop for futility
#We base on events rather than sample size

#Lets assume we have an interim look at 10 











