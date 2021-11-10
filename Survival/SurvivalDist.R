
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

  plot(fitc, conf.int = F, col="blue")
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



####Looking at piecewise distributions

GenerateData <- function(cp, max){
  control <- data.frame(time = c(runif(100, 0, cp), runif(100, cp, 0.75*max)), cens = rep(1, 100))
  treatment <- data.frame(time = c(runif(100, 0, cp),runif(100, cp, max)), cens = rep(1, 100))
  moxonidineData <- data.frame(time = c(control$time, treatment$time), cens = c(control$cens, treatment$cens), group = c(rep(1, 200), rep(2, 200)))
  return(moxonidineData)
}

moxonidineData = GenerateData(250, 500)

CPFunc <- function(moxonidineData){
  chivec = rep(NA, 491)
  for (i in 20:500){
    subset = moxonidineData[moxonidineData$time<i,]
    test = survdiff(Surv(time, cens) ~ group, subset)
    chivec[i] = test$chisq>qchisq(0.025, 1, lower.tail = F)
  }
  
  sigvec = rep(NA, 486)
  for (j in 21:490){
    if (chivec[j]&chivec[j+1]&chivec[j+2]){
      sigvec[j] = j
    }
  }
  sigvec = na.omit(sigvec)
  return(sigvec[1])
}

max = 500
cp = 150
outputvec = rep(NA, 10)
for (i in 1:10){
  outputvec[i] <- CPFunc(GenerateData(cp, max))
}

plot(outputvec)
outputvec


moxonidineData = GenerateData(150, 500)
fitc = survfit(Surv(time, cens)~group, data = moxonidineData)
plot(fitc, conf.int = F, col=c("blue", "red"))



#############################################################
#Go from here
#############################################################

library(survival)

par(mfrow=c(1,2))

GenerateData <- function(n, lambda){
  control <- data.frame(time = cumsum(rexp(n, 1)), cens = rep(1,n))
  treatment <- data.frame(time = c(cumsum(rexp(n/2,1)), (n/2)+cumsum(rexp(n/2, lambda))), cens=rep(1,n))
  moxonidineData <- data.frame(time = c(control$time, treatment$time), cens = c(control$cens, treatment$cens), group = c(rep(1, n), rep(2, n)))
  return(moxonidineData)
}

  moxonidineData <- GenerateData(200, 0.5)
  
  
  #Fit this data to a Kaplan-Meier estimate
  fitkm <- survfit(Surv(time, cens)~group, data = moxonidineData)
  plot(fitkm, conf.int = F, col=c("blue", "red"))
  abline(v = 100, col="green")
  
  
  #Fit a Weibull dist. to the control group 
  fitcontrol <- survreg(Surv(time, cens)~1, dist="weibull", data = moxonidineData[moxonidineData$group==1,])
  
  
  #Plot this Weibull dist. for the control group
  plot(x = predict(fitcontrol, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,],
       y = rev(seq(0.01, 0.99, by = 0.01)), type="l", xlab="Time", ylab="Survival",
       col = "blue", xlim=c(0,max(moxonidineData)*1.1))
  
  
  
  #Find out where the predictions are less than our changepoint
  PredictControl <- predict(fitcontrol, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,]
  for (i in 1:length(PredictControl)){
    if ((PredictControl[i]<100)==F){
      break
    }
  }
  
  probs <- seq(0.01, 0.01*i, by=.01)
  
  
  #The treatment line is the same as the control line up until the changepoint
  lines(x = predict(fitcontrol, type = "quantile", p =probs)[1,],
        y = rev(seq(1-0.01*length(probs), 0.99, by = 0.01)),
        col = "green")
  
  
  #Need to find the Weibull parameters that best fit the remaining treatment line
  #So need to extract the Kaplan-Meier estimates for the treatment group after changepoint
  
  TreatmentKMCP <- survfit(Surv(time, cens)~group, data = moxonidineData[moxonidineData$group==2,])
  TreatmentKMCPSurv <- summary(TreatmentKMCP)$surv[101:200]
  TreatmentKMCPTime <- summary(TreatmentKMCP)$time[101:200]
  
  fitcontrol <- survreg(Surv(time, cens)~1, dist="weibull", data = moxonidineData[moxonidineData$group==1,])
  
  PredictControl <- predict(fitcontrol, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,]
  for (i in 1:length(PredictControl)){
    if ((PredictControl[i]<100)==F){
      break
    }
  }
  
  probs <- seq(0.01, 0.01*i, by=.01)
  
  lambda = 0.5
  scale = (100+(100/lambda))*1.1
  ##Use optim instead here!
  
  
  LQFunc <- function(params){
    x = 100:floor(max(TreatmentKMCPTime))
    y = dweibull(x, shape=params[1], scale=scale)
    z = cumsum(y)
    
    diffsum = 0
    for (i in 1:length(x)){
      diffsum = diffsum + ((summary(TreatmentKMCP, times=x[i])$surv)-(1-z[i]-max(probs)))^2
    }
    return(diffsum)
  }
  
  optim1 <- optim(par=c(2.8), fn = LQFunc, method = "Brent", lower=1, upper=5)
  
  
  x = 100:floor(max(TreatmentKMCPTime))
  y = dweibull(x, shape = optim1$par[1], scale = optim1$par[2])
  z = cumsum(y)
  lines(x, 1-z-max(probs), col="green")
  
  optim1$par[1]
  optim1$par[2]






MyFunc <- function(lambda){
  moxonidineData <- GenerateData(200, lambda)
  TreatmentKMCP <- survfit(Surv(time, cens)~group, data = moxonidineData[moxonidineData$group==2,])
  TreatmentKMCPSurv <- summary(TreatmentKMCP)$surv[101:200]
  TreatmentKMCPTime <- summary(TreatmentKMCP)$time[101:200]
  LQFunc <- function(params){
    x = 100:floor(max(TreatmentKMCPTime))
    y = dweibull(x, shape=params[1], scale=params[2])
    z = cumsum(y)
    
    diffsum = 0
    for (i in 1:length(x)){
      diffsum = diffsum + ((summary(TreatmentKMCP, times=x[i])$surv)-(1-z[i]-max(probs)))^2
    }
    return(diffsum)
  }
  
  optim1 <- optim(par=c(2.9, 300), fn = LQFunc)
  return(list(shape = optim1$par[1], scale = optim1$par[2]))
}

shapevec = rep(NA, 10)
scalevec = rep(NA, 10)
for (i in 1:10){
  x = MyFunc()
  shapevec[i] = x$shape
  scalevec[i] = x$scale
}

plot(shapevec)
plot(scalevec)
plot(scalevec/shapevec)

mean(shapevec)
mean(scalevec)
mean(scalevec/shapevec)


lambda = seq(0.25, 0.75, by=0.125)
shape = c(2.11, 2.58, 2.76, 2.77, 3.46)
scale = c(608, 450, 343, 302, 258)
sos = c(294, 184, 128, 112, 77)

plot(lambda, sos)









