
library(nleqslv)
library(survival)


# Uncertainty -------------------------------------------------------------


#Or maybe we could put uncertainty around gamma1 and lambda1 and draw the same simulation curves to see how that would look? 
#Will probably look similar to the stuff I have been looking at 

bigT <- 3
lambda2 <- 0.06
gamma2 <- 0.8

lambda1 <- 0.02
gamma1 <- 1.2

chosenLength <- 60

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
SimMatrix1 <- matrix(NA, nrow = 500, ncol=length(time))

for (i in 1:500){
  bigT <- bigTVec[i]
  lambda1 <- lambda1vec[i]
  
  controltime <- seq(0, bigT, by=0.01)
  controlsurv <- exp(-(lambda2*controltime)^gamma2)
  
  treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
  treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))
  
  timecombined <- c(controltime, treatmenttime)[1:length(time)]
  survcombined <- c(controlsurv, treatmentsurv)[1:length(time)]
  
  
  SimMatrix1[i,] <- survcombined
  
}

lowerbound <- rep(NA, length(time))
upperbound <- rep(NA, length(time))
for (j in 1:length(time)){
  lowerbound[j] <- quantile(SimMatrix1[,j], 0.1, na.rm = T)
  upperbound[j] <- quantile(SimMatrix1[,j], 0.9, na.rm = T)
}

lines(time, lowerbound, col= "orange")
lines(time, upperbound,col= "orange")


# Assurance ---------------------------------------------------------------


#Now do assurance for both of the situations

#Do assurance for when we allow gamma1 to vary first

lambda2 <- 0.06
gamma2 <- 0.8

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


# Flexible ----------------------------------------------------------------


##Now we can look at assurance when we only elicit lambda1 (through the HR), but we allow gamma1 to vary in the assurance 
#This will hopefully give us a more flexible family of curves


timechosen1 <- 20
timechosen2 <- 50

assFunc <- function(n1, n2){
  assvec <- rep(NA, 100)
  for (i in 1:100){
    sampledpoint1 <- sample(na.omit(SimMatrix[,which(time==timechosen1)]), 1)
    sampledpoint2 <- sample(na.omit(SimMatrix[,which(time==timechosen2)]), 1)
    
    bigT <- rnorm(1, 3, sd = sqrt(0.1))
    
    dslnex <- function(x) {
      y <- numeric(2)
      y[1] <- exp(-(lambda2*bigT)^gamma2-x[1]^x[2]*(timechosen1^x[2]-bigT^x[2])) - sampledpoint1
      y[2] <- exp(-(lambda2*bigT)^gamma2-x[1]^x[2]*(timechosen2^x[2]-bigT^x[2])) - sampledpoint2
      y
    }
    
    
    xstart <- c(0.05,1)
    
    output <- nleqslv(xstart, dslnex)
    
    lambda1 <- output$x[1]
    gamma1 <- output$x[2]
    
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

lines(ssvec, predict(asssmooth), col = "green")







