GenerateData <- function(n, lambda){
  control <- data.frame(time = cumsum(rexp(n, 1)), cens = rep(1,n))
  treatment <- data.frame(time = c(cumsum(rexp(n/2,1)), (n/2)+cumsum(rexp(n/2, lambda))), cens=rep(1,n))
  moxonidineData <- data.frame(time = c(control$time, treatment$time), cens = c(control$cens, treatment$cens), group = c(rep(1, n), rep(2, n)))
  return(moxonidineData)
}

moxonidineData <- GenerateData(200, 0.5)

fitkm <- survfit(Surv(time, cens)~group, data = moxonidineData)
plot(fitkm, conf.int = F, col=c("blue", "red"))
abline(v = 100, col="green")

fitcontrol <- survreg(Surv(time, cens)~1, dist="weibull", data = moxonidineData[moxonidineData$group==1,])


#Find out where the predictions are less than our changepoint
PredictControl <- predict(fitcontrol, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,]
for (i in 1:length(PredictControl)){
  if ((PredictControl[i]<100)==F){
    break
  }
}

probs <- seq(0.01, 0.01*i, by=.01)



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
  scale = (100+(100/lambda))*1.1
  
  TreatmentKMCP <- survfit(Surv(time, cens)~group, data = moxonidineData[moxonidineData$group==2,])
  TreatmentKMCPSurv <- summary(TreatmentKMCP)$surv[101:200]
  TreatmentKMCPTime <- summary(TreatmentKMCP)$time[101:200]
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
  return(list(shape = optim1$par[1]))
}
MyFunc(lambda=0.5)


optim