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
HR1 <- 1
HR1Duration <- 0
HR2 <- 0.75
HR2Duration <- 100
lambdac <- log(0.5)/-12
lambdat <- lambdac*0.75
numPatients <- 340
numEventsRequired <- 512
NSims <- 1e5
recTime <- 12
#recNum <- (numPatients*2)/recTime

# Initialize result vectors
NoIApowerVec <- numeric(NSims)
WpowerVec <- numeric(NSims)
NoIAcensTimeVec <- numeric(NSims)
NoIASSVec <- numeric(NSims)
WCensTimeVec <- numeric(NSims)
WSSVec <- numeric(NSims)

# Generate control and treatment data
generateData <- function(lambda, group) {
  u <- runif(numPatients)
  data.frame(
    time = -log(u) / lambda,
    group = group
  )
}

# Function to perform censoring and analysis
censFunc <- function(dataset, numObs) {
  # Sort pseudo times and determine censoring time
  sortedPseudoTimes <- sort(dataset$pseudo_time)
  censTime <- sortedPseudoTimes[numObs]
  
  # Censor the observations
  dataset$status <- as.integer(dataset$pseudo_time <= censTime)
  
  # Only include patients enrolled by the censoring time
  dataset <- dataset[dataset$recTime <= censTime, ]
  
  # Calculate survival time
  dataset$survival_time <- ifelse(dataset$status == 1, dataset$time, censTime - dataset$recTime)
  
  sampleSize <- nrow(dataset)
  
  return(list(dataCombined = dataset, censTime = censTime, sampleSize = sampleSize))
}

results <- foreach(j = 1:NSims, .packages = c("survival", "dplyr")) %dopar% {
  controlData <- generateData(lambdac, "Control")
  treatmentData <- generateData(lambdat, "Treatment")
  dataCombined <- bind_rows(controlData, treatmentData)
  # result_vector <- c()
  # start_range <- seq(0, recTime-1)
  # end_range <- seq(1, recTime)
  # for (i in 1:length(start_range)) {
  #   random_values <- c(runif(recNum, start_range[i], end_range[i]))
  #   result_vector <- c(result_vector, random_values)
  # }
  #dataCombined$recTime <- sample(result_vector)
  dataCombined$recTime <- runif(numPatients*2, min = 0, max = recTime)
  dataCombined$pseudo_time <- dataCombined$time + dataCombined$recTime
  
  # Do it with no interim analysis first
  NoIAOutcome <- censFunc(dataCombined, numEventsRequired)
  NoIA <- NoIAOutcome$dataCombined
  NoIAcensTime <- NoIAOutcome$censTime
  NoIASS <- NoIAOutcome$sampleSize
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = NoIA)
  deltad <- as.numeric(exp(coef(coxmodel)))
  NoIALRT <- survdiff(Surv(time, status) ~ group, data = NoIA)
  NoIApower <- (NoIALRT$chisq > qchisq(0.95, 1) & deltad < 1)
  
  # Wieand rule
  WieandOutcome <- "Continue"
  W1Outcome <- censFunc(dataCombined, numEventsRequired * 0.5)
  W1 <- W1Outcome$dataCombined
  W1censTime <- W1Outcome$censTime
  W1SS <- W1Outcome$sampleSize
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = W1)
  deltad <- as.numeric(exp(coef(coxmodel)))
  if (deltad > 1) WieandOutcome <- "Stop1"
  
  W2Outcome <- censFunc(dataCombined, numEventsRequired * 0.75)
  W2 <- W2Outcome$dataCombined
  W2censTime <- W2Outcome$censTime
  W2SS <- W2Outcome$sampleSize
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = W2)
  deltad <- as.numeric(exp(coef(coxmodel)))
  if (deltad > 1) WieandOutcome <- "Stop2"
  
  WFinalOutcome <- censFunc(dataCombined, numEventsRequired)
  WFinal <- WFinalOutcome$dataCombined
  WFinalcensTime <- WFinalOutcome$censTime
  WFinalSS <- WFinalOutcome$sampleSize
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = WFinal)
  deltad <- as.numeric(exp(coef(coxmodel)))
  WLRT <- survdiff(Surv(time, status) ~ group, data = WFinal)
  Wpower <- (WLRT$chisq > qchisq(0.95, 1) & deltad < 1 & WieandOutcome == "Continue")
  WCensTime <- ifelse(WieandOutcome == "Stop1", W1censTime, ifelse(WieandOutcome == "Stop2", W2censTime, WFinalcensTime))
  WSS <- ifelse(WieandOutcome == "Stop1", W1SS, ifelse(WieandOutcome == "Stop2", W2SS, WFinalSS))
  
  #OBF rule
  OBFOutcome <- "Continue"
  OBF1Outcome <- censFunc(dataCombined, ceiling(numEventsRequired/3))
  OBF1 <- OBF1Outcome$dataCombined
  OBF1censTime <- OBF1Outcome$censTime
  OBF1SS <- OBF1Outcome$sampleSize
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = OBF1)
  sumcoxmodel <- summary(coxmodel)
  zscore <- -sumcoxmodel$coefficients[4]
  if (zscore < 0.011) OBFOutcome <- "Stop1"
  
  OBF2Outcome <- censFunc(dataCombined, ceiling(2*numEventsRequired/3))
  OBF2 <- OBF2Outcome$dataCombined
  OBF2censTime <- OBF2Outcome$censTime
  OBF2SS <- OBF2Outcome$sampleSize
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = OBF2)
  sumcoxmodel <- summary(coxmodel)
  zscore <- -sumcoxmodel$coefficients[4]
  if (zscore < 0.864) OBFOutcome <- "Stop2"
  
  OBFFinalOutcome <- censFunc(dataCombined, numEventsRequired)
  OBFFinal <- OBFFinalOutcome$dataCombined
  OBFFinalcensTime <- OBFFinalOutcome$censTime
  OBFFinalSS <- OBFFinalOutcome$sampleSize
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = OBFFinal)
  deltad <- as.numeric(exp(coef(coxmodel)))
  OBFLRT <- survdiff(Surv(time, status) ~ group, data = OBFFinal)
  OBFpower <- (OBFLRT$chisq > qchisq(0.95, 1) & deltad < 1 & OBFOutcome == "Continue")
  OBFCensTime <- ifelse(OBFOutcome == "Stop1", OBF1censTime, ifelse(OBFOutcome == "Stop2", OBF2censTime, OBFFinalcensTime))
  OBFSS <- ifelse(OBFOutcome == "Stop1", OBF1SS, ifelse(OBFOutcome == "Stop2", OBF2SS, OBFFinalSS))
  
  
  
  # Return results
  list(NoIApower = NoIApower, NoIAcensTime = NoIAcensTime, NoIASS = NoIASS, 
       Wpower = Wpower, WCensTime = WCensTime, WSS = WSS, OBFpower = OBFpower, OBFCensTime = OBFCensTime, OBFSS = OBFSS)
}

# Extract results
NoIApowerVec <- sapply(results, function(result) result$NoIApower)
NoIAcensTimeVec <- sapply(results, function(result) result$NoIAcensTime)
NoIASSVec <- sapply(results, function(result) result$NoIASS)
WpowerVec <- sapply(results, function(result) result$Wpower)
WCensTimeVec <- sapply(results, function(result) result$WCensTime)
WSSVec <- sapply(results, function(result) result$WSS)
OBFpowerVec <- sapply(results, function(result) result$OBFpower)
OBFCensTimeVec <- sapply(results, function(result) result$OBFCensTime)
OBFSSVec <- sapply(results, function(result) result$OBFSS)

# Clean up parallel resources
stopCluster(cl)

# Calculate the mean power and other statistics
mean_NoIApower <- mean(NoIApowerVec)
mean_NoIAcensTime <- mean(NoIAcensTimeVec)
mean_NoIASS <- mean(NoIASSVec)
mean_Wpower <- mean(WpowerVec)
mean_WCensTime <- mean(WCensTimeVec)
mean_WSS <- mean(WSSVec)
mean_OBFpower <- mean(OBFpowerVec)
mean_OBFCensTime <- mean(OBFCensTimeVec)
mean_OBFSS <- mean(OBFSSVec)

# Print results
cat("No IA power:", mean_NoIApower, "\n")
cat("No IA duration:", mean_NoIAcensTime, "\n")
cat("No IA sample size:", mean_NoIASS, "\n")
cat("Wieand power:", mean_Wpower, "\n")
cat("Wieand duration:", mean_WCensTime, "\n")
cat("Wieand sample size:", mean_WSS, "\n")
cat("OBF power:", mean_OBFpower, "\n")
cat("OBF duration:", mean_OBFCensTime, "\n")
cat("OBF sample size:", mean_OBFSS, "\n")

##############
#Calculating the changing HRs

lambdac <- log(0.5)/-12
HR1 <- 1.3
HR1Duration <- 3
HR2 <- 0.628
HR2Duration <- 100

controlTime <- seq(0, 100, by=0.01)
controlSurv <- exp(-lambdac*controlTime)

plot(controlTime, controlSurv, type = "l", ylim = c(0,1), col = "blue")

CP <- exp(-lambdac*HR1*HR1Duration)

treatmentTime1 <- seq(0, HR1Duration, by=0.01)
treatmentSurv1 <- exp(-lambdac*HR1*treatmentTime1)

treatmentTime2 <- seq(HR1Duration, HR1Duration+HR2Duration, by=0.01)
treatmentSurv2 <- exp(-lambdac*HR1Duration*HR1-lambdac*HR2*(treatmentTime2-HR1Duration))


lines(treatmentTime1, treatmentSurv1)
lines(treatmentTime2, treatmentSurv2)

