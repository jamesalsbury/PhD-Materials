#Here is the R script used to produce the figures/results from the DTE interim paper

library(truncnorm)
library(survival)
library(rjags)

# Posterior updating of the hazard ratio ----------------------------------------------------------------

#Let us show a real example of a clinical trial which exhibits DTE
# controldata <- read.csv(file = paste0("Papers/DTEPaper/data/Brahmer/IPD-control.csv"))
# treatmentdata <- read.csv(file = paste0("Papers/DTEPaper/data/Brahmer/IPD-treatment.csv"))
# 
# combinedData <- data.frame(time = c(controldata$Survival.time, treatmentdata$Survival.time), 
#                            status = c(controldata$Status, treatmentdata$Status), 
#                            group = c(rep("Control", nrow(controldata)), rep("Treatment", nrow(treatmentdata))))
# IATime <- 10
# combinedData$status <- combinedData$time<IATime
# 
# combinedData$time[combinedData$time>IATime] <- IATime
# 
# 
# kmfit <- survfit(Surv(time, status)~group, data = combinedData)
# plot(kmfit, conf.int = F, col=c("blue", "red"), xlab = "Time (months)", ylab = "Progression free survival (% of patients)", yaxt = "n")
# axis(2, at=seq(0, 1, by=0.1), labels=seq(0, 100, by=10))
# 
# n <- nrow(controldata)
# m <- nrow(combinedData)
# 
# #JAGS code which calculates posterior distributions
# 
# modelstring="
# 
# data {
#   for (j in 1:m){
#     zeros[j] <- 0
#   }
# }
# 
# model {
#   C <- 10000
#   for (i in 1:n){
#     zeros[i] ~ dpois(zeros.mean[i])
#     zeros.mean[i] <-  -l[i] + C
#     l[i] <- ifelse(datEvent[i]==1, log(gammac)+gammac*log(lambdac*datTimes[i])-(lambdac*datTimes[i])^gammac-log(datTimes[i]), -(lambdac*datTimes[i])^gammac)
#   }
#   for (i in (n+1):m){                                                                                                             
#     zeros[i] ~ dpois(zeros.mean[i])
#     zeros.mean[i] <-  -l[i] + C
#     l[i] <- ifelse(datEvent[i]==1, ifelse(datTimes[i]<bigT, log(gammac)+gammac*log(lambdac*datTimes[i])-(lambdac*datTimes[i])^gammac-log(datTimes[i]), log(gammac)+gammac*log(lambdat)+(gammac-1)*log(datTimes[i])-lambdat^gammac*(datTimes[i]^gammac-bigT^gammac)-(bigT*lambdac)^gammac), 
#       ifelse(datTimes[i]<bigT, -(lambdac*datTimes[i])^gammac, -(lambdac*bigT)^gammac-lambdat^gammac*(datTimes[i]^gammac-bigT^gammac)))
#   }
#   
#     lambdac ~ dbeta(1,1)T(0,)
#     gammac ~ dbeta(1,1)T(0,)
#     HR ~ dnorm(0.7, 10)T(0,)
#     bigT ~ dnorm(6, 1)T(0,)
#     lambdat <- lambdac*pow(HR, 1/gammac)
#     
#     }
# "
# 
# model = jags.model(textConnection(modelstring), data = list(datTimes = combinedData$time, datEvent = combinedData$status, n= n, m=m), quiet = T) 
# 
# 
# update(model, n.iter=500)
# output=coda.samples(model=model, variable.names=c("bigT", "HR"), n.iter = 5000)
# 
# plot(output)
# 
# output <- as.numeric(output[[1]])
# 
# hist(output[1:5000], xlab = "post-delay HR", freq = F, main = "Histogram of post-delay HR posterior", breaks = 30, xlim = c(0, 1))







#Defining the parameters used to simulate the data
HR <- rbeta(1, 10.8, 6.87)
bigT <- rnorm(1, mean = 6, sd = 1.5)
gammat <- gammac <- 0.8
lambdac <- 0.08
lambdat  <- lambdac*HR^(1/gammac)
n1 <- 500
n2 <- 500
IATime <- 30

#Simulating the control data
controldata <- rweibull(n1, gammac, 1/lambdac)
#Treatment
CP <- exp(-(lambdac*bigT)^gammac)[[1]]
u <- runif(n2)
suppressWarnings(treatmentdata <- ifelse(u>CP, (1/lambdac)*exp(1/gammac*log(-log(u))), ((1/(lambdat^gammac))*(-log(u)-(lambdac*bigT)^gammac)+bigT^gammac)^(1/gammac)))

combinedData <- data.frame(time = c(controldata, treatmentdata), status = rep(0, n1+n2), group = c(rep("Control", n1), rep("Treatment", n2)))

combinedData$status <- combinedData$time<IATime

combinedData$time[combinedData$time>IATime] <- IATime

kmfit <- survfit(Surv(time, status)~group, data = combinedData)

plot(kmfit, col = c("blue", "red"), xlab = "Time", ylab = "Survival")
legend("topright", legend = c("Control", "Treatment"), col = c("blue", "red"), lty=1)

n <- n1
m <- n1+n2

#JAGS code which calculates posterior distributions

modelstring="

data {
  for (j in 1:m){
    zeros[j] <- 0
  }
}

model {
  C <- 10000
  for (i in 1:n){
    zeros[i] ~ dpois(zeros.mean[i])
    zeros.mean[i] <-  -l[i] + C
    l[i] <- ifelse(datEvent[i]==1, log(gammac)+gammac*log(lambdac*datTimes[i])-(lambdac*datTimes[i])^gammac-log(datTimes[i]), -(lambdac*datTimes[i])^gammac)
  }
  for (i in (n+1):m){                                                                                                             
    zeros[i] ~ dpois(zeros.mean[i])
    zeros.mean[i] <-  -l[i] + C
    l[i] <- ifelse(datEvent[i]==1, ifelse(datTimes[i]<bigT, log(gammac)+gammac*log(lambdac*datTimes[i])-(lambdac*datTimes[i])^gammac-log(datTimes[i]), log(gammac)+gammac*log(lambdat)+(gammac-1)*log(datTimes[i])-lambdat^gammac*(datTimes[i]^gammac-bigT^gammac)-(bigT*lambdac)^gammac), 
      ifelse(datTimes[i]<bigT, -(lambdac*datTimes[i])^gammac, -(lambdac*bigT)^gammac-lambdat^gammac*(datTimes[i]^gammac-bigT^gammac)))
  }
  
    lambdac ~ dbeta(1,1)T(0,)
    gammac ~ dbeta(1,1)T(0,)
    HR ~ dnorm(0.7, 10)T(0,)
    bigT ~ dnorm(6, 1)T(0,)
    lambdat <- lambdac*pow(HR, 1/gammac)
    
    }
"

model = jags.model(textConnection(modelstring), data = list(datTimes = combinedData$time, datEvent = combinedData$status, n= n, m=m), quiet = T) 


update(model, n.iter=500)
output=coda.samples(model=model, variable.names=c("bigT", "HR"), n.iter = 5000)

output <- as.numeric(output[[1]])

hist(output[1:5000], xlab = "post-delay HR", freq = F, main = "Histogram of post-delay HR posterior", breaks = 30, xlim = c(0.45, 1))















