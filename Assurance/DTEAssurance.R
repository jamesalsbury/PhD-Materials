##We are eliciting as normal
#Then comparing the two scenarios
#Scenario 1: we do not allow gamma1 to vary
#Scenario 2: we are sampling two timepoints and fixing   

library(nleqslv)
library(survival)

# Plotting the control and median treatment curve -------------------------
bigTdist <- qnorm(0.5, mean = 3, sd = 0.741)
lambda2 <- 0.1
gamma2 <- 1.5


HRdist <- qbeta(0.5, 6.63, 4.5)
gamma1 <- 1.5

chosenLength <- 20

controltime <- seq(0, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)
plot(controltime, controlsurv, type="l", col="blue", ylab="Survival", xlab = "Time")

#Plotting the median treatment curve

treatmenttime <- seq(bigTdist, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)

lambda1 <- exp((1/gamma2)*log(HRdist)+log(lambda2))
treatmentsurv <- exp(-(lambda2*bigTdist)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigTdist^gamma1))

lines(treatmenttime, treatmentsurv, col="red")

legend("topright", legend = c("Control", "Treatment"), lty=1, col=c("blue", "red"))

# Showing how assurance works in the first scenario -------------------------
bigT <- rnorm(1, mean = 3, sd = 0.741)
HR <- rbeta(1, 6.63, 4.5)
lambda1 <- exp((1/gamma2)*log(HR)+log(lambda2))

treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))
lines(treatmenttime, treatmentsurv, col="red", lty=2)

# Showing how assurance works in the second scenario -------------------------
bigTVec <- rnorm(500, mean = 3, sd = 0.741)
HRVec <- rbeta(500, 6.63, 4.5)
lambda1vec <- exp((1/gamma2)*log(HRVec)+log(lambda2))

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

lines(time, lowerbound, col="red", lty=3)
lines(time, upperbound, col="red", lty=3)


#SHowing what goes on behind the scenes
bigTdist <- qnorm(0.5, mean = 3, sd = 0.741)
lambda2 <- 0.1
gamma2 <- 1.5

HRdist <- qbeta(0.5, 6.63, 4.5)
gamma1 <- 1.5

chosenLength <- 20

timechosen1 <- floor(0.4*chosenLength)
timechosen2 <- floor(0.75*chosenLength)

sampledpoint1 <- sample(na.omit(SimMatrix[,which(time==timechosen1)]), 1)
sampledpoint2 <- sample(na.omit(SimMatrix[,which(time==timechosen2)]), 1)

if (sampledpoint1>sampledpoint2){
  points(timechosen1, sampledpoint1, pch=19)
  points(timechosen2, sampledpoint2, pch=19)
  
  dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- exp(-(lambda2*bigTdist)^gamma2-x[1]^x[2]*(timechosen1^x[2]-bigTdist^x[2])) - sampledpoint1
    y[2] <- exp(-(lambda2*bigTdist)^gamma2-x[1]^x[2]*(timechosen2^x[2]-bigTdist^x[2])) - sampledpoint2
    y
  }
  
  
  xstart <- c(0.05,1)
  
  output <- nleqslv(xstart, dslnex)
  
  lambda1 <- output$x[1]
  gamma1 <- output$x[2]
  
  
  treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
  treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))
  
  lines(treatmenttime, treatmentsurv, col="yellow")
}



# Performing the assurance calculation for the first scenario -------------------------

lambda2 <- 0.1
gamma2 <- 1.5


gamma1 <- 1.5

chosenLength <- 20

assFunc <- function(n1, n2){
  assvec <- rep(NA, 100)
  for (i in 1:100){
    bigT <- rnorm(1, mean = 3, sd = 0.741)
    HR <- rbeta(1, 6.63, 4.5)
    lambda1 <- exp((1/gamma2)*log(HR)+log(lambda2))
    
    
    controldata <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))
    
    CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
    u <- runif(n2)
    
    suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
    
    DataCombined <- data.frame(time = c(controldata$time, z),
                               group = c(rep("Control", n1), rep("Treatment", n2)))
    
    
    DataCombined$cens <- DataCombined$time < chosenLength
    
    DataCombined$cens <- DataCombined$cens*1
    
    if (sum(DataCombined$cens)==(n1+n2)){
      
    } else {
      DataCombined[DataCombined$cens==0,]$time <- chosenLength
    }
    
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

plot(ssvec, predict(asssmooth), type="l",col = "green", ylim=c(0,1), xlab = "Total sample size", ylab="Assurance")


# Performing the assurance calculation for the second scenario -------------------------

timechosen1 <- floor(0.4*chosenLength)
timechosen2 <- floor(0.75*chosenLength)

bigTdist <- qnorm(0.5, mean = 3, sd = 0.741)
lambda2 <- 0.1
gamma2 <- 1.5

HRdist <- qbeta(0.5, 6.63, 4.5)
gamma1 <- 1.5

chosenLength <- 20


assFunc <- function(n1, n2){
  assvec <- rep(NA, 100)
  for (i in 1:100){
    sampledpoint1 <- sample(na.omit(SimMatrix[,which(time==timechosen1)]), 1)
    sampledpoint2 <- sample(na.omit(SimMatrix[,which(time==timechosen2)]), 1)
    
    if (sampledpoint1>sampledpoint2){
      bigT <- rnorm(1, mean = 3, sd = 0.741)
      
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
      
      if (lambda1<0){
        lambda1 <- 0.000001
      }
      
      controldata <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))
      
      CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
      u <- runif(n2)
      
      suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
      
      DataCombined <- data.frame(time = c(controldata$time, z),
                                 group = c(rep("Control", n1), rep("Treatment", n2)))
      
      DataCombined$cens <- DataCombined$time < chosenLength
      
      DataCombined$cens <- DataCombined$cens*1
      
      if (sum(DataCombined$cens)==(n1+n2)){
        
      } else {
        DataCombined[DataCombined$cens==0,]$time <- chosenLength
      }
      
      test <- survdiff(Surv(time, cens)~group, data = DataCombined)
      assvec[i] <- test$chisq > qchisq(0.95, 1)
    } 
    
    
  }
  mean(na.omit(assvec))
}

n1vec <- floor(seq(20, 200, length=40))
n2vec <- floor(seq(20, 200, length=40))

assVec <- rep(NA, length(n1vec))
for (i in 1:length(assVec)){
  assVec[i] <- assFunc(n1vec[i], n2vec[i])
}

ssvec <- n1vec + n2vec

asssmooth <- loess(assVec~ssvec)

lines(ssvec, predict(asssmooth), col = "orange")



