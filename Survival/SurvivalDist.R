
library(survival)

#Looking at how well a Weibull distribution can model our problem;
#cp is changepoint (where the treatment starts to take effect),
#Max is the maximum time for a patient to survive on treatment
WeibullCP <- function(cp, max){
  #Generating data
  control <- data.frame(time = c(runif(100, 0, cp), runif(100, cp, 0.75*max)), cens = rep(1, 100))
  treatment <- data.frame(time = c(runif(100, 0, cp),runif(100, cp, max)), cens = rep(1, 100))
  
  #Fitting a weibull model to both control/treatment groups
  fitcontrol <- survreg(Surv(time, cens)~1, dist="weibull", data = control)
  fittreatment <- survreg(Surv(time, cens)~1, dist="weibull", data = treatment)
  
  #Plotting the distribtions found
  plot(x = predict(fitcontrol, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,],
       y = rev(seq(0.01, 0.99, by = 0.01)), type="l", xlab="Time", ylab="Survival",
       col = "blue")
  
  lines(x = predict(fittreatment, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,],
        y = rev(seq(0.01, 0.99, by = 0.01)),
        col = "black")
  
  abline(v = cp, col="green")
  
  legend("topright", legend=c("Control", "Treatment", "Changepoint"), col=c("blue", "black", "green"),
         cex=0.8, lty=1)
}


#We can try a range of values that cp and max might take

WeibullCP(50, 500)
WeibullCP(100, 500)
WeibullCP(200, 500)
WeibullCP(300, 500)

#We see that it does not seem to be flexible enough to deal with the issues we want: for the curves to track each other and then separate at the changepoint
#When the changepoint is low (1/5) of the maximum then it does an okay job
#But as the changepoint gets larger the curves do not capture the changepoint very well at all


#We can see how the Kaplan-Meier estimate compares, function is very similar
KMCP = function(cp, max){
  control = data.frame(time = c(runif(100, 0, cp), runif(100, cp, 0.75*max)), cens = rep(1, 100))
  treatment = data.frame(time = c(runif(100, 0, cp),runif(100, cp, max)), cens = rep(1, 100))

  fitc = survfit(Surv(time, cens)~1, data = control)
  fitt = survfit(Surv(time, cens)~1, data = treatment)

  plot(fitc, conf.int = F, col="blue", xlim=c(0, max))
  lines(fitt, conf.int = F)
  abline(v = cp, col="green")
  legend("topright", legend=c("Control", "Treatment", "Changepoint"), col=c("blue", "black", "green"),
       cex=0.8, lty=1)
}

#Same range of values

KMCP(50, 500)
KMCP(100, 500)
KMCP(200, 500)
KMCP(300, 500)

#Does such a better job, the separation is nearly always on the changepoint, no matter how large we make the changepoint

