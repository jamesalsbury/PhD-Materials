library(survival)
library(dplyr)
library(foreach)
library(doParallel)

# Set the number of CPU cores you want to use
num_cores <- 16 # Change this to the number of cores you want to use

# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Parameters
lambdac <- 1/17.5
lambdat <- lambdac * 0.75
numPatients <- 340
numEventsRequired <- 512
NSims <- 1e5
powerVec <- numeric(NSims)
censVec <- numeric(NSims)

# Parallelize the simulation loop
results <- foreach(j = 1:NSims, .packages = c("survival", "dplyr")) %dopar% {
  
  # Generate control and treatment data
  generateData <- function(lambda, group) {
    u <- runif(numPatients)
    data.frame(
      time = -log(u) / lambda,
      group = group
    )
  }
  
  controlData <- generateData(lambdac, "Control")
  treatmentData <- generateData(lambdat, "Treatment")
  
  # Combine control and treatment data
  dataCombined <- bind_rows(controlData, treatmentData)
  
  # Sample recruitment time for each patient from a Uniform distribution
  dataCombined$recTime <- runif(numPatients, min = 0, max = 34)
  
  # Calculate pseudo event time
  dataCombined$pseudo_time <- dataCombined$time + dataCombined$recTime
  
  # Sort pseudo times and determine censoring time
  sortedPseudoTimes <- sort(dataCombined$pseudo_time)
  censTime <- sortedPseudoTimes[numEventsRequired]
  
  # Censor the observations
  dataCombined$status <- as.integer(dataCombined$pseudo_time <= censTime)
  
  # Only include patients enrolled by the censoring time
  dataCombined <- dataCombined[dataCombined$recTime <= censTime, ]
  
  # Calculate survival time
  dataCombined$survival_time <- pmin(censTime - dataCombined$recTime, dataCombined$time)
  
  # Making sure that the HR is less than 1
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = dataCombined)
  deltad <- as.numeric(exp(coef(coxmodel)))
  
  LRT <- survdiff(Surv(time, status) ~ group, data = dataCombined)
  power <- (LRT$chisq > qchisq(0.95, 1) & deltad < 1)
  
  # Return power as part of the result
  list(power = power, censTime = censTime)
}

# Combine results from parallel runs

powerVec <- sapply(results, function(result) result$power)
censTime <- sapply(results, function(result) result$censTime)

# Clean up parallel resources
stopCluster(cl)

# Calculate the mean power and mean censVec
mean_power <- mean(powerVec)
mean_censTime <- mean(censTime)

cat("Mean Power:", mean_power, "\n")
cat("Mean censVec:", mean_censTime, "\n")
