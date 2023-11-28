
lambdac <- log(2)/12
HR1 <- 0.75
T1 <- 0
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


# Set parameters
numSimulations <- 50
calTime <- seq(0, 100, by = 0.1)

# Initialize an empty matrix for storing results
eventMatrix <- matrix(NA, nrow = numSimulations, ncol = length(calTime))

# Simulate data and fill in the matrix
for (k in 1:numSimulations) {
  dataCombined <- generateData(lambdac, HR1, T1, HR2, numPatients, recTime)
  eventMatrix[k,] <- sapply(calTime, function(t) sum(dataCombined$pseudo_time < t))
}

# Calculate mean for each column
eventVec <- colMeans(eventMatrix)


numSimulations <- 50
calTime <- seq(0, 100, by = 0.1)
propMatrix <- matrix(NA, nrow = numSimulations, ncol = length(calTime))

for (k in 1:numSimulations){
  
  # Assuming dataCombined$pseudo_time is a vector of pseudo times
  dataCombined <- generateData(lambdac, HR1, T1, HR2, numPatients, recTime)
  
  # Create a function for the inner loop to avoid unnecessary copying
  censFuncInner <- function(data, threshold) {
    newDataCombined <- censFunc(data, threshold)$dataCombined
    newDataCombined <- newDataCombined[newDataCombined$status == 1, ]
    mean(newDataCombined$survival_time > 3)
  }
  
  # Use sapply for vectorized calculations
  propMatrix[k,] <- sapply(calTime, censFuncInner, data = dataCombined)
  
}


# Calculate mean for each column
eventProp <- colMeans(propMatrix)







plot(eventVec, eventProp, type = "l", xlab = "Number of events", ylab = "Proportion of events > 3 months", ylim = c(0,1))


 abline(h = 2/3, lty = 2)
# abline(v = 512*0.5)
# abline(v = 512*0.75)
# abline(v = 512)


#plot(calTime, eventProp, type = "l", xlab = "Calendar time (months)", ylab = "Proportion of events > 3 months", ylim = c(0,1))
# abline(h = 2/3)
# abline(v = 10)
# abline(v = 20)
# abline(v = 30)
# abline(v = 40)
# abline(v = 50)






