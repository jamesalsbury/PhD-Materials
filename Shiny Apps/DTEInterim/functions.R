SimDTEDataSet <- function(n, lambdac, bigT, HRStar, recTime) {
  
  #Simulate the control times
  u <- runif(n)
  controlTime <- -log(u)/lambdac
  
  #Simulate the treatment times
  CP <- exp(-lambdac*bigT)
  u <- runif(n)
  treatmentTime <- ifelse(u>CP, -log(u)/lambdac, (1/(HRStar*lambdac))*(HRStar*lambdac*bigT-log(u)-lambdac*bigT))

  #Combine the two groups
  dataCombined <- data.frame(time = c(controlTime, treatmentTime),
                             group = c(rep("Control", n), rep("Treatment", n)))
  
  #Add on a random recruitment time
  dataCombined$recTime <- runif(n*2, 0, recTime)
  
  #Calculate the pseudo time of the event
  dataCombined$pseudoTime <- dataCombined$time + dataCombined$recTime
  
  return(dataCombined)
}

CensFunc <- function(dataCombined, numEvents) {
  
  dataCombined <- dataCombined[order(dataCombined$pseudoTime), ]
  
  censTime <- dataCombined$pseudoTime[numEvents]
  
  dataCombined$status <- dataCombined$pseudoTime <= censTime
  dataCombined$status <- dataCombined$status * 1
  dataCombined$enrolled <- dataCombined$recTime < censTime
  dataCombined <- dataCombined[dataCombined$enrolled, ]
  dataCombined$survival_time <- ifelse(dataCombined$pseudoTime > censTime,
                                       censTime - dataCombined$recTime,
                                       dataCombined$time)
  
  return(list(dataCombined = dataCombined, censTime = censTime, SS = nrow(dataCombined)))
}

interimLookFunc <- function(dataCombined, observedHR){
  
  coxmodel <- coxph(Surv(survival_time, status)~group, data = dataCombined)
  deltad <- as.numeric(exp(coef(coxmodel)))
  Outcome <- "Continue"
  if (deltad>observedHR){ Outcome <- "Stop"}
  return(Outcome)
}
