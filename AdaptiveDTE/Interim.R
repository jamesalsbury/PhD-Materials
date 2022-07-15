#Libraries
library(plyr)
library(survival)
library(tidyverse)


#We change the underlying true values of T and HR and see how that affects our judgements

#We will keep HR fixed for now, at 0.625

#We will fix T at 4.04

lambda2 <- 0.08
gamma2 <- 0.8
bigT <- 4.04
gamma1 <- 0.8
HR <- 0.625
lambda1 <- exp((log(HR)/gamma2)+log(lambda2))

n1 <- 490/2
n2 <- 490/2

CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
u <- runif(n2)

suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))

DataCombined <- data.frame(time = c(rweibull(n1, gamma2, 1/lambda2),z),
                           group = c(rep("Control", n1), rep("Treatment", n2)))



DataCombined <- DataCombined[order(DataCombined$time),]

IATime <- DataCombined[344,]$time

DataCombined$event <- DataCombined$time<IATime

DataCombined[DataCombined$event==0,]$time <- IATime

kmfit <- survfit(Surv(time, event)~group, data = DataCombined)
plot(kmfit, xlim=c(0,60), col=c("blue", "red"))


#How to estimate the HR

#We need to use the control data to estimate lambda2 and gamma2
controlSample <- DataCombined %>%
  filter(group=="Control")
weibfit <- survreg(Surv(time, event)~1, data = controlSample, dist = "weibull")
lambda2est <- as.numeric(1/(exp(weibfit$icoef[1])))
gamma2est <- as.numeric(exp(-weibfit$icoef[2]))

#We can use these estimated parameters to simulate the remaining control data

#There are 245 patients in the control group
#How many have had an event?

sum(controlSample$time<IATime)

#How many are yet to have an event?
controleventsremaining <- n1 - sum(controlSample$time<IATime)

u <- runif(controleventsremaining, 0, controleventsremaining/n1)


gamma1est <- gamma2est

#We now need to estimate lambda1, which in turn estimates the HR

controltimeest <- seq(0, 60, by=0.01)
controlsurvest <- exp(-(lambda2est*controltimeest)^gamma2est)
lines(controltimeest, controlsurvest)














