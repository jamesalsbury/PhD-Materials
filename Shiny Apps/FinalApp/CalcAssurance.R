
library(nleqslv)
library(survival)
#Need to look at what we actually need to compute assurance??
#We are sampling from the control and treatment curves
#But maybe we only need number of events that occur (or do not occur) by the time-point that the assurance is going to take place>

#We could compute assurance in the way that we have done previously
#Then find a 'quicker' way by just estimating the number of events


#We assume the control is lambda2 = 0.06, gamma2 = 0.8

lambda2 <- 0.06
gamma2 <- 0.8

controltime <- seq(0, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)
plot(controltime, controlsurv, type="l", col="blue", ylab="Survival", xlab = "Time")

#We assume T comes from a Normal(3, 0.01) dist

x <- seq(0, 6, by=0.01)
y <- dnorm(x, mean = 3, sd = sqrt(0.01))
plot(x, y, type="l")

#We assume HR comes from a Beta(6, 4) dist

x <- seq(0, 1, by=0.01)
y <- dbeta(x, 6, 4)
plot(x, y, type="l")

#Plot a sampled treatment curve

bigT <- rnorm(1, 3, sd = sqrt(0.01))
HR <- rbeta(1, 6, 4)
lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
gamma1 <- gamma2

treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))

lines(treatmenttime, treatmentsurv, col="red")

#Calculate assurance the way we have been doing
#Firstly we assume n1 = 100, n2 = 100
#We assume the trial will be stopped after 40 months
#0 = alive
#1 = dead

n1 <- 100
n2 <- 100
assvec <- rep(NA, 10000)
lambda2 <- 0.06
gamma2 <- 0.8
gamma1 <- gamma2
triallength <-15

for (i in 1:10000){
  bigT <- rnorm(1, 3, sd = sqrt(0.1))
  HR <- rbeta(1, 6, 4)
  lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
  CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
  u <- runif(n2)
  
  suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
  
  dataCombined <- data.frame(Time = c(rweibull(n1, gamma2, 1/lambda2), z), Group = c(rep("Control", n1), rep("Treatment", n2)), 
                             Cens = rep(1, n1+n2))
  
  dataCombined$Cens <-  dataCombined$Time < triallength
  dataCombined$Cens <- dataCombined$Cens*1
  test <- survdiff(Surv(Time, Cens)~Group, data = dataCombined)
  pvalue <- 1-pchisq(test$chisq, 1)
  assvec[i] <-pvalue< 0.05
}

mean(assvec)

# 'True' assurance is 0.591

#Can we replicate this by just looking at the survival probabilities?

#We can calculate how many patients we expect to be alive after 40 months in the control group fairly easily





#Have the same process for elicitation but see whether we can make assurance more flexible by allowing gamma1 to vary in the confidence intervals
bigT <- rnorm(1, 3, sd = sqrt(0.1))
HR <- rbeta(1, 6, 4)
lambda1 <- exp((log(HR)/gamma2)+log(lambda2))

lambda2 <- 0.06
gamma2 <- 0.8

bigTMedian <- qnorm(0.5, 3, sqrt(0.1))
HRMedian <- qbeta(0.5, 6, 4)
lambda1Median <- exp((log(HRMedian)/gamma2)+log(lambda2))

controltime <- seq(0, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)
plot(controltime, controlsurv, type="l", col="blue", ylab="Survival", xlab = "Time")

gamma1 <- gamma2

treatmenttime <- seq(bigTMedian, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
treatmentsurv <- exp(-(lambda2*bigTMedian)^gamma2 - lambda1Median^gamma1*(treatmenttime^gamma1-bigTMedian^gamma1))

lines(treatmenttime, treatmentsurv, col="red")


##Drawing sim curves
gamma1 <- gamma2

bigTVec <- rnorm(500, 3, sd = sqrt(0.1))
HRVec <- rbeta(500, 6, 4)

time <- seq(0, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
SimMatrix <- matrix(NA, nrow = 500, ncol=length(time))

for (i in 1:500){
  bigT <- bigTVec[i]
  HR <- HRVec[i]
  lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
  
  controltime <- seq(0, bigT, by=0.01)
  controlsurv <- exp(-(lambda2*controltime)^gamma2)
  
  treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
  treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))
  
  timecombined <- c(controltime, treatmenttime)[1:length(time)]
  survcombined <- c(controlsurv, treatmentsurv)[1:length(time)]
  
  
  SimMatrix[i,] <- survcombined
  
}

lowerbound <- rep(NA, length(time))
upperbound <- rep(NA, length(time))
for (j in 1:length(time)){
  lowerbound[j] <- quantile(SimMatrix[,j], 0.1)
  upperbound[j] <- quantile(SimMatrix[,j], 0.9)
}

lines(time, lowerbound)
lines(time, upperbound)

#Now sample a random y axis from time = 40? For example

timechosen <- 40

#hist(SimMatrix[,which(time==timechosen)], breaks = 50, xlim=c(0,0.8))

sampledpoint <- sample(SimMatrix[,which(time==timechosen)], 1)

bigT <- rnorm(1, 3, sd = sqrt(0.1))

sampledlambda1 <- exp((1/gamma1)*(log((-log(sampledpoint) - (lambda2*bigT)^gamma2)/(timechosen^gamma1-bigT^gamma1))))

points(timechosen, sampledpoint)

treatmenttime <- seq(3, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
treatmentsurv <- exp(-(lambda2*3)^gamma2 - sampledlambda1^gamma1*(treatmenttime^gamma1-3^gamma1))
lines(treatmenttime, treatmentsurv)

#Now we sample two time points, at time = 30 and time  = 60

timechosen1 <- 30
timechosen2 <- 60

sampledpoint1 <- sample(SimMatrix[,which(time==timechosen1)], 1)
sampledpoint2 <- sample(SimMatrix[,which(time==timechosen2)], 1)

points(timechosen1, sampledpoint1)
points(timechosen2, sampledpoint2)

bigT <- rnorm(1, 3, sd = sqrt(0.1))

dslnex <- function(x) {
  y <- numeric(2)
  y[1] <- exp(-(lambda2*bigT)^gamma2-x[1]^x[2]*(timechosen1^x[2]-bigT^x[2])) - sampledpoint1
  y[2] <- exp(-(lambda2*bigT)^gamma2-x[1]^x[2]*(timechosen2^x[2]-bigT^x[2])) - sampledpoint2
  y
}


xstart <- c(0.05,1)

output <- nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta="newton"))

treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - output$x[1]^output$x[2]*(treatmenttime^output$x[2]-bigT^output$x[2]))

lines(treatmenttime, treatmentsurv, col="red")


#Now time to actually calculate assurance

#We have 5 simulated data sets, calculate assurance in each case

#Firstly we will calculate assurance with not letting gamma1 vary in both elicitation and assurance

#This will change depending on what data set we are working with
inputtedData <- simData1

n1 <- 100
n2 <- 100
assvec <- rep(NA, 100)
controldata <- inputtedData[inputtedData$group=="control",]
weibcontrolfit <- survreg(Surv(time, cens)~1, data = controldata, dist = "weibull")
gamma2 <- as.numeric(exp(-weibcontrolfit$icoef[2]))
lambda2 <- as.numeric(1/(exp(weibcontrolfit$icoef[1])))
gamma1 <- gamma2
triallength <-40

for (i in 1:100){
  bigT <- rnorm(1, 3, sd = sqrt(0.1))
  HR <- rbeta(1, 6, 4)
  lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
  CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
  u <- runif(n2)
  
  suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
  
  dataCombined <- data.frame(Time = c(rweibull(n1, gamma2, 1/lambda2), z), Group = c(rep("Control", n1), rep("Treatment", n2)), 
                             Cens = rep(1, n1+n2))
  
  dataCombined$Cens <-  dataCombined$Time < triallength
  dataCombined$Cens <- dataCombined$Cens*1
  test <- survdiff(Surv(Time, Cens)~Group, data = dataCombined)
  pvalue <- 1-pchisq(test$chisq, 1)
  assvec[i] <-pvalue< 0.05
}

mean(assvec)


#Or maybe we could put uncertainty around gamma1 and lambda1 and draw the same simulation curves to see how that would look? 
#Will probably look similar to the stuff I have been looking at 





bigT <- 3
lambda2 <- 0.06
gamma2 <- 0.8

lambda1 <- 0.02
gamma1 <- 1.2

chosenLength <- 40

controltime <- seq(0, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)
plot(controltime, controlsurv, type="l", col="blue", ylab="Survival", xlab = "Time")

treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))

lines(treatmenttime, treatmentsurv, col="red")

#Drawing sim curves
lambda2 <- 0.06
gamma2 <- 0.8

bigTVec <- rnorm(500, 3, sd = sqrt(0.1))
lambda1vec <- rnorm(500, 0.02, sd = sqrt(0.00005))
gamma1vec <- rnorm(500, 1.2, sd = sqrt(0.01))

time <- seq(0, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
SimMatrix <- matrix(NA, nrow = 500, ncol=length(time))

for (i in 1:500){
  bigT <- bigTVec[i]
  lambda1 <- lambda1vec[i]
  gamma1 <- gamma1vec[i]
  
  controltime <- seq(0, bigT, by=0.01)
  controlsurv <- exp(-(lambda2*controltime)^gamma2)
  
  treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
  treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))
  
  timecombined <- c(controltime, treatmenttime)[1:length(time)]
  survcombined <- c(controlsurv, treatmentsurv)[1:length(time)]
  
  
  SimMatrix[i,] <- survcombined
  
}

lowerbound <- rep(NA, length(time))
upperbound <- rep(NA, length(time))
for (j in 1:length(time)){
  lowerbound[j] <- quantile(SimMatrix[,j], 0.1, na.rm = T)
  upperbound[j] <- quantile(SimMatrix[,j], 0.9, na.rm = T)
}

lines(time, lowerbound)
lines(time, upperbound)


#How can we best approximate this with setting gamma1 = gamma2?
#We can use least squares to approximate this?

gamma1 <- 1.2
lambda1 <- 0.02 


treatmenttime <- seq(3, chosenLength)
treatmentsurv <- exp(-(lambda2*3)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-3^gamma1))

optimfunc1 <- function(par){
  diff <- 0
  for (i in 1:length(treatmenttime)){
    y <- exp(-(lambda2*3)^gamma2 - par^gamma2*(treatmenttime[i]^gamma2-3^gamma2))
    diff <- diff + (y - treatmentsurv[i])^2
  }
  return(diff)
}


optimoutput1 <-  optimize(f = optimfunc1, interval = c(0,1))$minimum

treatmenttime1 <- seq(3, 120, by=0.01)
treatmentsurv1 <- exp(-(lambda2*3)^gamma2 - 0.02^gamma2*(treatmenttime1^gamma2-3^gamma2))

#plot(controltime, controlsurv, type="l", col="blue", ylab="Survival", xlab = "Time")

lines(treatmenttime1, treatmentsurv1, col="green")

#Draw sim lines
lambda2 <- 0.06
gamma2 <- 0.8

bigTVec <- rnorm(500, 3, sd = sqrt(0.1))
lambda1vec <- rnorm(500, optimoutput1, sd = sqrt(0.0001))

gamma1 <- gamma2

time <- seq(0, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
SimMatrix <- matrix(NA, nrow = 500, ncol=length(time))

for (i in 1:500){
  bigT <- bigTVec[i]
  lambda1 <- lambda1vec[i]
  
  controltime <- seq(0, bigT, by=0.01)
  controlsurv <- exp(-(lambda2*controltime)^gamma2)
  
  treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
  treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))
  
  timecombined <- c(controltime, treatmenttime)[1:length(time)]
  survcombined <- c(controlsurv, treatmentsurv)[1:length(time)]
  
  
  SimMatrix[i,] <- survcombined
  
}

lowerbound <- rep(NA, length(time))
upperbound <- rep(NA, length(time))
for (j in 1:length(time)){
  lowerbound[j] <- quantile(SimMatrix[,j], 0.1, na.rm = T)
  upperbound[j] <- quantile(SimMatrix[,j], 0.9, na.rm = T)
}

lines(time, lowerbound, col= "orange")
lines(time, upperbound,col= "orange")


#Now do assurance for both of the situations

#Do assurance for when we allow gamma1 to vary first

lambda2 <- 0.06
gamma2 <- 0.8

chosenLength <- 40

assFunc <- function(n1, n2){
  assvec <- rep(NA, 500)
  for (i in 1:500){
    bigT <- rnorm(1, 3, sd = sqrt(0.1))
    lambda1 <- rnorm(1, 0.02, sd = sqrt(0.00005))
    gamma1 <- rnorm(1, 1.2, sd = sqrt(0.01))
    
    controldata <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))
    
    CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
    u <- runif(n2)
    
    suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
    
    DataCombined <- data.frame(time = c(controldata$time, z),
                               group = c(rep("Control", n1), rep("Treatment", n2)))
    
    
    DataCombined$cens <- DataCombined$time < chosenLength
    
    DataCombined$cens <- DataCombined$cens*1
    
    test <- survdiff(Surv(time, cens)~group, data = DataCombined)
    assvec[i] <- test$chisq > qchisq(0.95, 1)
  }
  mean(assvec)
}


n1vec <- floor(seq(20, 200, length=40))
n2vec <- floor(seq(20, 200, length=40))

assVec <- rep(NA, length(n1vec))
for (i in 1:length(assVec)){
  assVec[i] <- assFunc(n1vec[i], n2vec[i])
}

ssvec <- n1vec + n2vec

asssmooth <- loess(assVec~ssvec)

plot(ssvec, predict(asssmooth), type="l", ylim=c(0,1), xlab="Sample size", ylab="Assurance")

#Now do the situation where we do not allow gamma1 to vary
gamma1 <- gamma2

assFunc <- function(n1, n2){
  assvec <- rep(NA, 500)
  for (i in 1:500){
    bigT <- rnorm(1, 3, sd = sqrt(0.1))
    lambda1 <- rnorm(1, optimoutput1, sd = sqrt(0.00005))
    
    
    controldata <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))
    
    CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
    u <- runif(n2)
    
    suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
    
    DataCombined <- data.frame(time = c(controldata$time, z),
                               group = c(rep("Control", n1), rep("Treatment", n2)))
    
    
    DataCombined$cens <- DataCombined$time < chosenLength
    
    DataCombined$cens <- DataCombined$cens*1
    
    test <- survdiff(Surv(time, cens)~group, data = DataCombined)
    assvec[i] <- test$chisq > qchisq(0.95, 1)
  }
  mean(assvec)
}

n1vec <- floor(seq(20, 200, length=40))
n2vec <- floor(seq(20, 200, length=40))

assVec <- rep(NA, length(n1vec))
for (i in 1:length(assVec)){
  assVec[i] <- assFunc(n1vec[i], n2vec[i])
}

ssvec <- n1vec + n2vec

asssmooth <- loess(assVec~ssvec)

lines(ssvec, predict(asssmooth), col = "red")







