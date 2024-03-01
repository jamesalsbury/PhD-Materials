
lambdac <- log(2)/12
HR1 <- 1
T1 <- 3
HR2 <- 0.75
numPatients <- 340
recTime <- 34


generateData <- function(lambdac, HR1, T1, HR2, numPatients, recTime) {
  CP <- exp(-lambdac*HR1*T1)
  u <- runif(numPatients)
  controlData <- -log(u)/lambdac
  u <- runif(numPatients)
  treatmentData <- ifelse(u>CP, -log(u)/(lambdac*HR1), (1/(lambdac*HR2))*(lambdac*HR2*T1-log(u)-lambdac*HR1*T1))
  dataCombined <- data.frame(time = c(controlData, treatmentData), group = c(rep("Control", numPatients), rep("Treatment", numPatients)))
  dataCombined$recTime <- runif(numPatients*2, min = 0, max = recTime)
  dataCombined$pseudo_time <- dataCombined$time + dataCombined$recTime
  return(dataCombined)
}

# Create a function for the inner loop to avoid unnecessary copying
censFuncInner <- function(data, threshold) {
  newDataCombined <- censFunc(data, threshold)$dataCombined
  newDataCombined <- newDataCombined[newDataCombined$status == 1, ]
  mean(newDataCombined$survival_time > 3)
}


# Function to perform censoring and analysis
censFunc <- function(dataset, censTime) {
  
  # Censor the observations
  dataset$status <- as.integer(dataset$pseudo_time <= censTime)
  
  # Only include patients enrolled by the censoring time
  dataset <- dataset[dataset$recTime <= censTime, ]
  
  # Calculate survival time
  dataset$survival_time <- ifelse(dataset$status == 1, dataset$time, censTime - dataset$recTime)
  
  sampleSize <- nrow(dataset)
  
  return(list(dataCombined = dataset, censTime = censTime, sampleSize = sampleSize))
}



propEventFunc <- function(lambdac, HR1, T1, HR2, numPatients, recTime){
  # Set parameters
  numSimulations <- 1
  eventVec <- 1:(numPatients*2)
  
  # Initialize an empty matrix for storing results
  eventMatrix <- matrix(NA, nrow = numSimulations, ncol = length(eventVec))
  
  # Simulate data and fill in the matrix
  for (k in 1:numSimulations) {
    dataCombined <- generateData(lambdac, HR1, T1, HR2, numPatients, recTime)
  }
  
  dataCombined <- dataCombined[order(dataCombined$pseudo_time), ]
  
  calTime <- dataCombined$pseudo_time
  
  
  propMatrix <- matrix(NA, nrow = numSimulations, ncol = length(calTime))
  
  for (k in 1:numSimulations){
    
    # Use sapply for vectorized calculations
    propMatrix[k,] <- sapply(calTime, censFuncInner, data = dataCombined)
    
  }
  
  # Calculate mean for each column
  eventProp <- colMeans(propMatrix)
  
  return(list(eventProp = eventProp, eventVec = eventVec, calTime = calTime))
  
}



x1 <- propEventFunc(lambdac, HR1, T1, HR2, numPatients, recTime)
plotTime <- seq(1, max(x1$calTime), by = 6)
plotEvents <- rep(NA, length(plotTime))
for (i in 1:length(plotTime)){
  plotEvents[i] <-    sum(x1$calTime<plotTime[i])
}
plot(x1$eventVec, x1$eventProp, type = "l", ylim = c(0, 1), xlab = "Number of events", 
     ylab = "Proportion of events > 3 months", col = "red", main = paste0("Recruitment: ", recTime, " months"))
points(x1$eventVec[plotEvents], x1$eventProp[plotEvents])
abline(h = 2/3, lty = 2)
 abline(v = 512*0.5, lty = 3)
 abline(v = 512*0.75, lty = 3)
 print(plotTime)
 
 
 
# abline(v = 512*0.75, lty = 3)
# abline(v = 512, lty = 3)
# 
# 
# 
# plot(x1$calTime, x1$eventVec, type = "l")
# abline(v = 14)
# abline(h = 512)
# 
# x1 <- propEventFunc(lambdac, HR1, T1, HR2, numPatients, 0)
# x2 <- propEventFunc(lambdac, HR1, T1, HR2, numPatients, 12)
# x3 <- propEventFunc(lambdac, HR1, T1, HR2, numPatients, 34)
# x4 <- propEventFunc(lambdac, HR1, T1, HR2, numPatients, 60)
# 
# plot(x1$calTime, x1$eventVec, type = "l", xlab = "Calendar time", ylab = "Number of events")
# lines(x2$calTime, x2$eventVec, lty = 2, col = "blue")
# lines(x3$calTime, x3$eventVec, lty = 3, col = "red")
# lines(x4$calTime, x4$eventVec, lty = 4, col = "green")
# legend("bottomright", legend = c("R = 0", "R = 12", "R = 34", "R = 60"), col = c("black", "blue", "red", "green"), lty = 1:4)
# 
# plot(x1$calTime, x1$eventProp, type = "l", xlab = "Calendar time", ylab = "Proportion of events > 3 months")
# lines(x2$calTime, x2$eventProp, lty = 2, col = "blue")
# lines(x3$calTime, x3$eventProp, lty = 3, col = "red")
# lines(x4$calTime, x4$eventProp, lty = 4, col = "green")
# legend("bottomright", legend = c("R = 0", "R = 12", "R = 34", "R = 60"), col = c("black", "blue", "red", "green"), lty = 1:4)
# abline(h = 2/3, lty = 2)
# 
# 
# plot(x1$eventVec, x1$eventProp, type = "l", ylim = c(0, 1), xlab = "Number of events", ylab = "Proportion of events > 3 months")
# lines(x2$eventVec, x2$eventProp, col = "blue", lty = 2)
# lines(x3$eventVec, x3$eventProp, col = "red", lty = 3)
# lines(x4$eventVec, x4$eventProp, col = "green", lty = 4)
# legend("bottomright", legend = c("R = 0", "R = 12", "R = 34", "R = 60"), col = c("black", "blue", "red", "green"), lty = 1:4)
# abline(h = 2/3, lty = 2)
# abline(v = 512*0.5, lty = 3)
# abline(v = 512*0.75, lty = 3)
# abline(v = 512, lty = 3)













