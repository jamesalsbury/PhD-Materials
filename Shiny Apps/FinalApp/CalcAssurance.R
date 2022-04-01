
library(nleqslv)
library(survival)


# Uncertainty -------------------------------------------------------------

#Plotting the control 
bigT <- 3
lambda2 <- 0.06
gamma2 <- 0.8

lambda1 <- 0.02
gamma1 <- 1.5

chosenLength <- 60

controltime <- seq(0, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)
plot(controltime, controlsurv, type="l", col="blue", ylab="Survival", xlab = "Time")

#Plotting the median treatment lines

#This situation is where gamma1 is allowed to vary (we are using this as a control - for example a clinician might think this is appropriate 
#so we are seeing what we are losing by enforcing gamma1 = gamma2)

treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))

lines(treatmenttime, treatmentsurv, col="red")

#Calculating the confidence intervals in this situation
lambda2 <- 0.06
gamma2 <- 0.8

bigTVec <- rnorm(500, 3, sd = sqrt(0.1))
lambda1vec <- rnorm(500, 0.02, sd = sqrt(0.00005))
gamma1vec <- rnorm(500, gamma1, sd = sqrt(0.01))

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

lines(time, lowerbound, col="red", lty=3)
lines(time, upperbound, col="red", lty=3)


#Resetting gamma1 and lambda1
#This is now the situation where we are optimising the treatment curve - but we are not allowing gamma1 to vary anymore

gamma1 <- 1.5
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

#We are finding the lambda1 which corresponds to the closest median treatment curve as above
optimoutput1 <-  optimize(f = optimfunc1, interval = c(0,1))$minimum

treatmenttime1 <- seq(3, 120, by=0.01)
treatmentsurv1 <- exp(-(lambda2*3)^gamma2 - 0.02^gamma2*(treatmenttime1^gamma2-3^gamma2))

lines(treatmenttime1, treatmentsurv1, col="green")

#Calculating the confidence intervals for this situation where gamma1 = gamma2
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

lines(time, lowerbound, col= "green", lty=3)
lines(time, upperbound,col= "green", lty=3)

legend("topright", legend = c("Control", "Gamma1 not varying", "Gamma1 varying"), col=c("blue", "green", "red"), lty=1)

# Assurance ---------------------------------------------------------------

#The only one we will have access to is the green simulation lines
#SimMatrix1

#We calculate assurance for the situation where gamma1 is allowed to vary
#We would not normally know this - but we are calculating as a benchmark

lambda2 <- 0.06
gamma2 <- 0.8

assFunc <- function(n1, n2){
  assvec <- rep(NA, 100)
  for (i in 1:100){
    bigT <- rnorm(1, 3, sd = sqrt(0.1))
    lambda1 <- truncnorm::rtruncnorm(1, a = 0, mean = 0.02, sd = sqrt(0.00005))
    gamma1 <- rnorm(1, 1.5, sd = sqrt(0.01))
    
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

plot(ssvec, predict(asssmooth), type="l", ylim=c(0,1), xlab="Sample size", ylab="Assurance", col="red")

#Now do the situation where we do not allow gamma1 to vary in the elicitation
#This is the calculation we have been doing thus far

gamma1 <- gamma2

assFunc <- function(n1, n2){
  assvec <- rep(NA, 100)
  for (i in 1:100){
    bigT <- rnorm(1, 3, sd = sqrt(0.1))
    lambda1 <- truncnorm::rtruncnorm(1, a = 0, mean = optimoutput1, sd = sqrt(0.00005))

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

lines(ssvec, predict(asssmooth), col = "green")


# Flexible ----------------------------------------------------------------

#We can look at more flexible models to calculate assurance
#In this situation we have enforced gamma1 = gamma2 in the elicitation
#But we can pick two time-points and sample the treatment curves from these
#Then fit a Weibull parameterisation where gamma1 does not equal gamma2

#This section shows what goes on behind the scenes

# timechosen1 <- 20
# timechosen2 <- 50
# 
# sampledpoint1 <- sample(na.omit(SimMatrix1[,which(time==timechosen1)]), 1)
# sampledpoint2 <- sample(na.omit(SimMatrix1[,which(time==timechosen2)]), 1)
# 
# points(timechosen1, sampledpoint1, pch=19)
# points(timechosen2, sampledpoint2, pch=19)
# 
# dslnex <- function(x) {
#   y <- numeric(2)
#   y[1] <- exp(-(lambda2*bigT)^gamma2-x[1]^x[2]*(timechosen1^x[2]-bigT^x[2])) - sampledpoint1
#   y[2] <- exp(-(lambda2*bigT)^gamma2-x[1]^x[2]*(timechosen2^x[2]-bigT^x[2])) - sampledpoint2
#   y
# }
# 
# 
# xstart <- c(0.05,1)
# 
# output <- nleqslv(xstart, dslnex)
# 
# lambda1 <- output$x[1]
# gamma1 <- output$x[2]
# 
# 
# treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
# treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))
# 
# lines(treatmenttime, treatmentsurv, col="yellow")

#This section does the same as above but calculates assurance in this case
#Picks two random y points at t=20 and t=50 and then draws a survival curve according to this
#Performs a log-rank test with these data


timechosen1 <- floor(0.4*chosenLength)
timechosen2 <- floor(0.75*chosenLength)


assFunc <- function(n1, n2){
  assvec <- rep(NA, 100)
  for (i in 1:100){
    sampledpoint1 <- sample(na.omit(SimMatrix1[,which(time==timechosen1)]), 1)
    sampledpoint2 <- sample(na.omit(SimMatrix1[,which(time==timechosen2)]), 1)
    
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

lines(ssvec, predict(asssmooth), col = "orange")

legend("bottomright", legend = c("Benchmark", "Gamma1 not varying", "Two timepoints"), col=c("red", "green", "orange"), lty=1)


##Now we need to do the non parametric case
#We split the x-axis up into 20 slices, from 20 to 60?
#Uniformly sample from the slice between the 0.1 and 0.9 confidence intervals
#Make it a monotonically decreasing function by bounding it above by the previous y value sampled


#Show it on a graph to ensure we are doing the right thing throughout

xvec <- seq(20, 100, by=10)
yvec <- rep(NA, length(xvec))

yvec[1] <- runif(1, min = quantile(na.omit(SimMatrix1[,which(time==xvec[1])]), 0.1), max = quantile(na.omit(SimMatrix1[,which(time==xvec[1])]), 0.9))

for (i in 2:length(xvec)){
  yvec[i] <- runif(1, min = quantile(na.omit(SimMatrix1[,which(time==xvec[i])]), 0.1), max =  yvec[i-1])
}


lines(xvec, yvec)

#Looks like we are doing the right thing, but we are much more likely to hit the bottom of the confidence interval
#It only takes 1 random y-axis draw to be low and then we are stuck being low

#Maybe we can sample the y-axis in the same way as when we allow gamma1 to vary in the assurance, but plot a non-parametric line instead of a Weibull line

timechosen1 <- floor(0.4*chosenLength)
timechosen2 <- chosenLength

sampledpoint1 <- sample(na.omit(SimMatrix1[,which(time==timechosen1)]), 1)
sampledpoint2 <- sample(na.omit(SimMatrix1[,which(time==timechosen2)]), 1)

if (sampledpoint1>sampledpoint2){
  points(timechosen1, sampledpoint1)
  points(timechosen2, sampledpoint2)
  xvec <- seq(timechosen1, timechosen2, by=5)
  yvec <- rep(NA, length(xvec))
  yvec[1] <- sampledpoint1
  
  for (i in 2:(length(xvec)-1)){
    yvec[i] <- runif(1, min = sampledpoint2, max =  yvec[i-1])
  }
  yvec[length(xvec)] <- sampledpoint2
  lines(xvec, yvec)
}



