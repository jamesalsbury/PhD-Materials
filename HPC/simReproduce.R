library(survival)
library(dplyr)

#set.seed(123)  # Set a seed for reproducibility

# Parameters
lambdac <- 1/17.5
lambdat <- lambdac * 0.75
numPatients <- 340
numEventsRequired <- 512
NSims <- 1e4
powerVec <- rep(NA, NSims)
censVec <- rep(NA, NSims)

for (j in 1:NSims){
  
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
  
  #Making sure that the HR is less than 1
  coxmodel <- coxph(Surv(survival_time, status)~group, data = dataCombined)
  deltad <- as.numeric(exp(coef(coxmodel)))
  
  LRT <- survdiff(Surv(time, status)~group, data = dataCombined)
  powerVec[j] <- (LRT$chisq > qchisq(0.95, 1) & deltad<1)
  censVec[j] <- censTime
  
}

mean(powerVec)
mean(censVec)

