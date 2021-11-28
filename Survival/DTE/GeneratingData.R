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
lambda1 <- 0.03
gamma1 <- 0.8

#Our sample sizes for both groups; 1=control, 2=treatment
n1 <- 500
n2 <- 500

#Simulate data for the control group
simcontrol <- data.frame(time = rweibull(n1, gamma2, 1/lambda2), cens = rep(1, n1))
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

SimTreatment <- SimTreatmentFunc(n2)
fittreatmentKM <- survfit(Surv(time, cens)~1, data = SimTreatment)
lines(fittreatmentKM, conf.int = F, col="red")



bigT = 5
lambda2 = 0.06
gamma2 = 0.8
lambda1 = 0.03
gamma1 = 0.8
t = 5
CDF = rep(NA, length(t))
for (i in 1:length(t)){
  if (t[i]>bigT){
    CDF[i] = 1-exp(-(lambda2*bigT)^(gamma2)-lambda1^gamma1*(t[i]^gamma1-bigT^gamma1))
  } else{
    CDF[i] = 1-exp(-(lambda2*t[i])^gamma2)
  }
}



u = runif(1000)
samples = rep(NA, 1000)
for (i in 1:length(u)){
  if (u[i]>CDF){
    samples[i] = (1/lambda1)*exp(log(-log(1-u[i])-(lambda2*T)^gamma2+lambda1^gamma1*T^gamma1)/gamma1)
  } else{
    samples[i] = (1/lambda2)*exp(log(-log(1-u[i]))/gamma2)
  }
}
#plot(samples)
  
data <- data.frame(time = samples, cens = rep(1, 1000))
fittreatmentKM <- survfit(Surv(time, cens)~1, data = data)
lines(fittreatmentKM, conf.int = F, col="red")
library(survival)

