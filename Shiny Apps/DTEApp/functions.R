SimDTEDataSet <- function(n1, n2, gamma1, gamma2, lambda1, lambda2, bigT, recTime, censTime){
  #Simulates the treatment data
  CP <- exp(-(lambda2*bigT)^gamma2)
  u <- runif(n2)
  
  treatmenttime <- ifelse(u>CP, (1/lambda2)*(-log(u))^(1/gamma2), (1/(lambda1^gamma1)*((lambda1*bigT)^gamma1-log(u)-(lambda2*bigT)^gamma2))^(1/gamma1))
  
  dataCombined <- data.frame(time = c(rweibull(n1, gamma2, 1/lambda2), treatmenttime),
                             group = c(rep("Control", n1), rep("Treatment", n2)))
  
  
  #Adds a random uniformly distributed value, based on the recruitment time
  dataCombined$time <- dataCombined$time + runif(n1+n2, min = 0, max = recTime)
  
  #If the time is less than the total trial length time then the event has happened
  dataCombined$event <- dataCombined$time < censTime
  
  #Making it a binary value (rather than T/F), for ease to read
  dataCombined$event <- dataCombined$event*1
  
  #Need to set the event time to be the censoring time
  if (sum(dataCombined$event)==(n1+n2)){
    
  } else{
    dataCombined[dataCombined$time>censTime,]$time <- censTime
  }
  
  return(dataCombined)
}

# n1 <- 10000
# n2 <- 10000
# gamma1 <- 0.8
# gamma2 <- 0.8
# lambda1 <- 0.04
# lambda2 <- 0.08
# bigT <- 6
# recTime <- 6
# censTime <- 60
# 
# y <- SimDTEDataSet(n1, n2, gamma1, gamma2, lambda1, lambda2, bigT, recTime, censTime)
# 
# coxmodel <- coxph(Surv(time, event)~group, data = y)
# coxmodel
# 
# exp(coef(coxmodel)) 
# 
# CI <- exp(confint(coxmodel))
# 
# CI[2]
# 
# HR <- (lambda1/lambda2)^gamma2
# 
# (sum(y$time<6)+(n1+n2-sum(y$time<6))*HR)/(n1+n2)

