library(survival)
lambdac <- 0.08
lambdat <- 0.04
bigT <- 5
time <- seq(0, 60, by=0.001)
controlData <- exp(-lambdac*time)
plot(time, controlData, type="l", col="blue")
treatmentData <- ifelse(time<bigT, exp(-lambdac*time), exp(-(lambdac*bigT)-lambdat*(time-bigT)) )
lines(time, treatmentData, col="red")

controldensity <- lambdac*exp(-lambdac*time)
plot(time, controldensity)

#Can use inversion sampling

u <- runif(10000, 0, 0.08)

correstime <- -(1/lambdac)*log(u/lambdac)

controlData <- data.frame(time= correstime)

lambda <- seq(0, 0.5, by=0.001)
likelihood <- 100*log(lambda)-lambda*sum(correstime)
plot(lambda, likelihood)

10000/sum(correstime)


treatmentDensity <- ifelse(time<=bigT, lambdac*exp(-lambdac*time), lambdat*exp(-lambdat*(time-bigT)-bigT*lambdac))
plot(time,treatmentDensity)

controlfit <- survfit(Surv(time)~1, data = controlData)
lines(controlfit, conf.int = F)

CP <- lambdat*exp(-lambdac*bigT)

u <- runif(10000, 0, 0.08)

correstime <- ifelse(u<CP, -(1/lambdac)*log(u/lambdac), (1/lambdat)*(lambdat*bigT-log(u/lambdat)-bigT*lambdac))


treatmentData <- data.frame(time= correstime)


treatmentfit <-  survfit(Surv(time)~1, data = treatmentData)
lines(treatmentfit, conf.int = F)

u <- 0.04
(1/lambdat)*(lambdat*bigT-log(u/lambdat)-bigT*lambdac)
