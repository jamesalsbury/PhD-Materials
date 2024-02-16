library(survival)


P_S <- 0.7
P_DTE <- 0.6

# T ~ Gamma(2.84, 0.519)
# HR* ~ Gamma(25.6, 33.8)

n <- 340
num_Events <- 512
lambdac <- log(2)/12
recTime <- 34


#Now to simulate some treatment curves

NRep <- 5e4

# Sample probabilities
u1 <- runif(NRep)
u2 <- runif(NRep)

# Sample HRStar
HRStar <- ifelse(u1 > P_S, 1, rgamma(NRep, 25.6, 33.8))

# Sample bigT based on HRStar values
bigT <- ifelse(u1 <= P_S & u2 <= P_DTE & HRStar != 1, rgamma(NRep, 2.84, 0.519), 0)
treatmentSamplesDF <- data.frame(bigT, HRStar)

#Calculate assurance

NRep <- 5000
assvec <- rep(NA, NRep)

for (i in 1:NRep){
  # Sample control times
  u <- runif(n)
  controlTime <- -log(u)/lambdac
  
  #Compute treatment times
  HRStar <- sample(treatmentSamplesDF[,2], 1)
  bigT <- sample(treatmentSamplesDF[,1], 1)
  
  CP <- exp(-lambdac*bigT)
  u <- runif(n)
  treatmentTime <- ifelse(u>CP, -log(u)/lambdac, (1/(HRStar*lambdac))*(HRStar*lambdac*bigT-log(u)-lambdac*bigT))
  
  combinedData <- data.frame(time = c(controlTime, treatmentTime), 
                             group = c(rep("Control", n), rep("Treatment", n)))
  
  combinedData$recTime <- runif(n*2, 0, recTime)
  
  combinedData$pseudoTime <- combinedData$time + combinedData$recTime
  
  
  
  # Censor this data set at 512 events
  
  # Assuming combinedData is already defined
  combinedData <- combinedData[order(combinedData$pseudoTime), ]
  
  censTime <- combinedData$pseudoTime[num_Events]
  
  combinedData$status <- combinedData$pseudoTime <= censTime
  combinedData$status <- combinedData$status * 1
  combinedData$enrolled <- combinedData$recTime < censTime
  combinedData <- combinedData[combinedData$enrolled, ]
  combinedData$survival_time <- ifelse(combinedData$pseudoTime > censTime,
                                       censTime - combinedData$recTime,
                                       combinedData$time)
  
  #Do a log-rank test on these data
  
 test <- survdiff(Surv(survival_time, status)~group, data = combinedData)
 
 coxmodel <- coxph(Surv(survival_time, status)~group, data = combinedData)
 deltad <- as.numeric(exp(coef(coxmodel)))
 
 assvec[i] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
 
}

mean(assvec)



# 
# x <- seq(0, 15, by=0.01)
# y <- dgamma(x, 2.84, 0.519)
# 
# plot(x,y,type = "l")













