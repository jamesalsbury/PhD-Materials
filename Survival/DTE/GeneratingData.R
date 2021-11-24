library(survival)
#################################################################################
#Looking at splitting survival curve up
################################################################################# 
#Simulating data
#We have lambda2, gamma2 from historical data

lambda2 <- 0.06
gamma2 <- 0.8
  
#We have elicited bigT, lambda1, gamma1 from questioning our experts
  
bigT <- 10
lambda1 <- 0.015
gamma1 <- 0.5

#Our sample sizes for both groups; 1=control, 2=treatment
n1 <- 500
n2 <- 500

#Simulate data for the control group
simcontrol <- data.frame(time = rweibull(100, gamma2, 1/lambda2), cens = rep(1, 100))
#Plot this on a KM curve
fitcontrolKM <- survfit(Surv(time, cens)~1, data = simcontrol)
plot(fitcontrolKM, conf.int = F, xlim=c(0,max(simcontrol$time)*1.1), ylab="Survival", xlab="Time (months)", col="blue",
       main = "Simulated data")
legend("topright", legend = c("KM curve to the control data","KM curve to the treatment data") , col=c("blue", "red"), lty=1)
  
#Fit a Weibull to the control data
fitcontrol <- survreg(Surv(time, cens)~1, dist="weibull", data = simcontrol)

#Plotting the Weibull distribution found for the control data
# lines(x = predict(fitcontrol, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,],
#        y = rev(seq(0.01, 0.99, by = 0.01)), type="l", xlab="Time (months)", ylab="Survival",
#        col = "blue", xlim=c(0,100))
#   
  

#Fitting hypothetical lines for the treatment; before changepoint
effectt <- seq(0, bigT, by=0.01)
effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
#lines(effectt, effecty, col="green", lty=3)
  
#Fitting hypothetical lines for the treatment; after changepoint
aftereffectt <- seq(bigT, max(simcontrol$time)*1.1, by=0.01)
aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*bigT)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-bigT^gamma1))
#lines(aftereffectt, aftereffecty, col="red", lty=2)
  
#Combining both before and after the changepoint
treatmentcurvet <- c(effectt, aftereffectt)
treatmentcurvey <- c(effecty, aftereffecty)
  

#Simulating data for the treatment group
SimTreatmentFunc <- function(n2){
  
  split <- seq(1, 0, by=-(1/(n2/10)))
  
  events <- rep(NA, n2)
  for (i in 1:(length(split)-1)){
    if (i==1){
      events[((i-1)*10+1):(i*10)] <- runif(10, 0, treatmentcurvet[sum(split[i+1]<treatmentcurvey)])
    } 
    else{
      events[((i-1)*10+1):(i*10)] <- runif(10, treatmentcurvet[sum(split[i]<treatmentcurvey)], treatmentcurvet[sum(split[i+1]<treatmentcurvey)])
    }
  }
  
  data <- data.frame(time = events, cens = rep(1, n2))
}

SimTreatment <- SimTreatmentFunc(100)
fittreatmentKM <- survfit(Surv(time, cens)~1, data = SimTreatment)
lines(fittreatmentKM, conf.int = F, col="red")
  



