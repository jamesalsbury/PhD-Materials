library(survival)

GenerateData <- function(n, lambda){
  control <- data.frame(time = cumsum(rexp(n, 1)), cens = rep(1,n))
  treatment <- data.frame(time = c(cumsum(rexp(n/2,1)), (n/2)+cumsum(rexp(n/2, lambda))), cens=rep(1,n))
  moxonidineData <- data.frame(time = c(control$time, treatment$time), cens = c(control$cens, treatment$cens), group = c(rep(1, n), rep(2, n)))
  return(moxonidineData)
}

MyFunc <- function(lambda){
moxonidineData <- GenerateData(200, lambda)


fitcontrol <- survreg(Surv(time, cens)~1, dist="weibull", data = moxonidineData[moxonidineData$group==1,])

#Find out where the predictions are less than our changepoint
PredictControl <- predict(fitcontrol, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,]
for (i in 1:length(PredictControl)){
  if ((PredictControl[i]<100)==F){
    break
  }
}

probs <- seq(0.01, 0.01*i, by=.01)



TreatmentKMCP <- survfit(Surv(time, cens)~group, data = moxonidineData[moxonidineData$group==2,])
TreatmentKMCPSurv <- summary(TreatmentKMCP)$surv[101:200]
TreatmentKMCPTime <- summary(TreatmentKMCP)$time[101:200]

scale = summary(TreatmentKMCP)$time[190]

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

optim1 <- optimise(f = LQFunc, lower=1, upper=15)
return(optim1$minimum)
}

CalcFunc <- function(lambda){
  outvec <- rep(NA, 50)
  for (i in 1:50){
    outvec[i] = MyFunc(lambda)
  }
  return(mean(outvec))
}

lambdavec = seq(0.10, 0.7, by=0.05)
assvec = rep(NA, 13)
for (j in 1:13){
  assvec[j] = CalcFunc(lambdavec[j])
}


plot(lambdavec[1:9], assvec[1:9], ylim=c(3,7))

lo <- loess(assvec~lambdavec)
plot(predict(lo))

lm(assvec[1:9]~lambdavec[1:9])

#Plots!

par(mfrow=c(1,2))
GenerateData <- function(cp, max){
  control <- data.frame(time = c(runif(100, 0, cp), runif(100, cp, 0.75*max)), cens = rep(1, 100))
  treatment <- data.frame(time = c(runif(100, 0, cp),runif(100, cp, max)), cens = rep(1, 100))
  moxonidineData <- data.frame(time = c(control$time, treatment$time), cens = c(control$cens, treatment$cens), group = c(rep(1, 200), rep(2, 200)))
  return(moxonidineData)
}

cp=200
moxonidineData <- GenerateData(cp, 500)

fitkm <- survfit(Surv(time, cens)~group, data = moxonidineData)
plot(fitkm, conf.int = F, col=c("blue", "red"))
abline(v = cp, col="green")

fitcontrol <- survreg(Surv(time, cens)~1, dist="weibull", data = moxonidineData[moxonidineData$group==1,])


#Plot this Weibull dist. for the control group
plot(x = predict(fitcontrol, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,],
     y = rev(seq(0.01, 0.99, by = 0.01)), type="l", xlab="Time", ylab="Survival",
     col = "blue", xlim=c(0,max(moxonidineData)*1.1))

PredictControl <- predict(fitcontrol, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,]
for (i in 1:length(PredictControl)){
  if ((PredictControl[i]<cp)==F){
    break
  }
}

probs <- seq(0.01, 0.01*i, by=.01)


#The treatment line is the same as the control line up until the changepoint
lines(x = predict(fitcontrol, type = "quantile", p =probs)[1,],
      y = rev(seq(1-0.01*length(probs), 0.99, by = 0.01)),
      col = "green")

TreatmentKMCP <- survfit(Surv(time, cens)~group, data = moxonidineData[moxonidineData$group==2,])
TreatmentKMCPSurv <- summary(TreatmentKMCP)$surv[101:200]
TreatmentKMCPTime <- summary(TreatmentKMCP)$time[101:200]

scale = summary(TreatmentKMCP)$time[200]


LQFunc <- function(params){
  x = cp:floor(max(TreatmentKMCPTime))
  y = dweibull(x, shape=params[1], scale=scale)
  z = cumsum(y)
  
  diffsum = 0
  for (i in 1:length(x)){
    diffsum = diffsum + ((summary(TreatmentKMCP, times=x[i])$surv)-(1-z[i]-max(probs)))^2
  }
  return(diffsum)
}
optim1 <- optimise(f = LQFunc, lower=1, upper=15)

x = cp:floor(max(TreatmentKMCPTime))
y = dweibull(x, shape=optim1$minimum, scale=scale)
z = cumsum(y)
lines(x, 1-z-max(probs))








