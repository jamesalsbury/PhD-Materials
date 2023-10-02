library(survival)
library(dplyr)
library(foreach)
library(doParallel)

KFScenarioList <- list(
  A = list(
    HR1 = 0.75,
    T1 = 1000,
    HR2 = 0.75,
    T2 = 1000
),
  B = list(
    HR1 = 1,
    T1 = 1000,
    HR2 = 1,
    T2 = 1000
),
  C = list(
    HR1 = 1.3,
    T1 = 1000,
    HR2 = 1.3,
    T2 = 1000
),
  D = list(
    HR1 = 1,
    T1 = 3,
    HR2 = 0.693,
    T2 = 1000
),
  E = list(
    HR1 = 1,
    T1 = 6,
    HR2 = 0.62,
    T2 = 1000
),
  F = list(
    HR1 = 1.3,
    T1 = 3,
    HR2 = 0.628,
    T2 = 1000
))

ScenarioList <- list(
  A = list(
    HR1 = 0.55,
    T1 = 1000,
    HR2 = 0.55,
    T2 = 1000
  ),
  B = list(
    HR1 = 0.65,
    T1 = 1000,
    HR2 = 0.65,
    T2 = 1000
  ),
  C = list(
    HR1 = 0.75,
    T1 = 1000,
    HR2 = 0.75,
    T2 = 1000
  ),
  D = list(
    HR1 = 0.85,
    T1 = 1000,
    HR2 = 0.85,
    T2 = 1000
  ),
  E = list(
    HR1 = 0.95,
    T1 = 1000,
    HR2 = 0.95,
    T2 = 1000
  ),
  F = list(
    HR1 = 1,
    T1 = 1000,
    HR2 = 1,
    T2 = 1000
  ),
  G = list(
    HR1 = 1.1,
    T1 = 1000,
    HR2 = 1.1,
    T2 = 1000),
  H = list(
    HR1 = 1.2,
    T1 = 1000,
    HR2 = 1.2,
    T2 = 1000),
  I = list(
    HR1 = 1.3,
    T1 = 1000,
    HR2 = 1.3,
    T2 = 1000),
  J = list(
    HR1 = 1,
    T1 = 3,
    HR2 = 0.75,
    T2 = 1000),
  K = list(
    HR1 = 1,
    T1 = 6,
    HR2 = 0.75,
    T2 = 1000),
  L = list(
    HR1 = 1.3,
    T1 = 3,
    HR2 = 0.628,
    T2 = 1000)
)


paramsList <- list(
  recTime = 12,
  numPatients = 340,
  lambdac = -log(0.5)/12,
  numEventsRequired = 512,
  NSims = 1e3
)

# Generate control and treatment data
generateData <- function(lambdac, HR1, T1, HR2, T2, numPatients, recTime) {
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

for (i in 1:length(ScenarioList)){
  
  # Set the number of CPU cores you want to use
  num_cores <- 8 # Change this to the number of cores you want to use
  
  # Register parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  results <- foreach(j = 1:paramsList$NSims, .packages = c("survival", "dplyr")) %dopar% {
    
    dataCombined <- generateData(paramsList$lambdac, ScenarioList[[i]]$HR1, 
                                 ScenarioList[[i]]$T1, ScenarioList[[i]]$HR2, ScenarioList[[i]]$T2, paramsList$numPatients, paramsList$recTime)
    
    # Do it with no interim analysis first
    NoIAOutcome <- censFunc(dataCombined, paramsList$numEventsRequired)
    NoIA <- NoIAOutcome$dataCombined
    NoIAcensTime <- NoIAOutcome$censTime
    NoIASS <- NoIAOutcome$sampleSize
    coxmodel <- coxph(Surv(survival_time, status) ~ group, data = NoIA)
    deltad <- as.numeric(exp(coef(coxmodel)))
    NoIALRT <- survdiff(Surv(time, status) ~ group, data = NoIA)
    NoIApower <- (NoIALRT$chisq > qchisq(0.95, 1) & deltad < 1)
    
    # Wieand rule
    WieandOutcome <- "Continue"
    W1Outcome <- censFunc(dataCombined, paramsList$numEventsRequired * 0.5)
    W1 <- W1Outcome$dataCombined
    W1censTime <- W1Outcome$censTime
    W1SS <- W1Outcome$sampleSize
    coxmodel <- coxph(Surv(survival_time, status) ~ group, data = W1)
    deltad <- as.numeric(exp(coef(coxmodel)))
    if (deltad > 1) WieandOutcome <- "Stop1"
    
    W2Outcome <- censFunc(dataCombined, paramsList$numEventsRequired * 0.75)
    W2 <- W2Outcome$dataCombined
    W2censTime <- W2Outcome$censTime
    W2SS <- W2Outcome$sampleSize
    coxmodel <- coxph(Surv(survival_time, status) ~ group, data = W2)
    deltad <- as.numeric(exp(coef(coxmodel)))
    if (deltad > 1 & WieandOutcome=="Continue") WieandOutcome <- "Stop2"
    
    WFinalOutcome <- censFunc(dataCombined, paramsList$numEventsRequired)
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
    OBF1Outcome <- censFunc(dataCombined, ceiling(paramsList$numEventsRequired/3))
    OBF1 <- OBF1Outcome$dataCombined
    OBF1censTime <- OBF1Outcome$censTime
    OBF1SS <- OBF1Outcome$sampleSize
    coxmodel <- coxph(Surv(survival_time, status) ~ group, data = OBF1)
    sumcoxmodel <- summary(coxmodel)
    zscore <- -sumcoxmodel$coefficients[4]
    if (zscore < 0.011) OBFOutcome <- "Stop1"
    
    OBF2Outcome <- censFunc(dataCombined, ceiling(2*paramsList$numEventsRequired/3))
    OBF2 <- OBF2Outcome$dataCombined
    OBF2censTime <- OBF2Outcome$censTime
    OBF2SS <- OBF2Outcome$sampleSize
    coxmodel <- coxph(Surv(survival_time, status) ~ group, data = OBF2)
    sumcoxmodel <- summary(coxmodel)
    zscore <- -sumcoxmodel$coefficients[4]
    if (zscore < 0.864 & OBFOutcome=="Continue") OBFOutcome <- "Stop2"
    
    OBFFinalOutcome <- censFunc(dataCombined, paramsList$numEventsRequired)
    OBFFinal <- OBFFinalOutcome$dataCombined
    OBFFinalcensTime <- OBFFinalOutcome$censTime
    OBFFinalSS <- OBFFinalOutcome$sampleSize
    coxmodel <- coxph(Surv(survival_time, status) ~ group, data = OBFFinal)
    deltad <- as.numeric(exp(coef(coxmodel)))
    OBFLRT <- survdiff(Surv(time, status) ~ group, data = OBFFinal)
    OBFpower <- (OBFLRT$chisq > qchisq(0.95, 1) & deltad < 1 & OBFOutcome == "Continue")
    OBFCensTime <- ifelse(OBFOutcome == "Stop1", OBF1censTime, ifelse(OBFOutcome == "Stop2", OBF2censTime, OBFFinalcensTime))
    OBFSS <- ifelse(OBFOutcome == "Stop1", OBF1SS, ifelse(OBFOutcome == "Stop2", OBF2SS, OBFFinalSS))
    
    #Proposed rule
    sortedPseudoTimes <- sort(dataCombined$pseudo_time)
    propVec <- numeric(paramsList$numEventsRequired)
    cutoff_index <- ceiling(paramsList$numEventsRequired/2)
    
    for (k in cutoff_index:paramsList$numEventsRequired){
      censTime <- sortedPseudoTimes[k]
      status <- as.integer(dataCombined$pseudo_time <= censTime)
      filteredData <- dataCombined[dataCombined$recTime <= censTime, ]
      survival_time <- ifelse(status == 1, filteredData$time, censTime - filteredData$recTime)
      propVec[k] <- mean(survival_time[status==1]>3)
    }
    
    propVec[1:(cutoff_index-1)] <- 0
    
    Stop1 <- max(paramsList$numEventsRequired*0.5, which.min(propVec<(2/3)))
    Stop2 <- max(paramsList$numEventsRequired*0.75, which.min(propVec<(2/3)))
    
    PropOutcome <- "Continue"
    Prop1Outcome <- censFunc(dataCombined, Stop1)
    Prop1 <- Prop1Outcome$dataCombined
    Prop1censTime <- Prop1Outcome$censTime
    Prop1SS <- Prop1Outcome$sampleSize
    coxmodel <- coxph(Surv(survival_time, status) ~ group, data = Prop1)
    deltad <- as.numeric(exp(coef(coxmodel)))
    if (deltad > 1) PropOutcome <- "Stop1"
    
    Prop2Outcome <- censFunc(dataCombined, Stop2)
    Prop2 <- Prop2Outcome$dataCombined
    Prop2censTime <- Prop2Outcome$censTime
    Prop2SS <- Prop2Outcome$sampleSize
    coxmodel <- coxph(Surv(survival_time, status) ~ group, data = Prop2)
    deltad <- as.numeric(exp(coef(coxmodel)))
    if (deltad > 1 & PropOutcome=="Continue") PropOutcome <- "Stop2"
    
    PropFinalOutcome <- censFunc(dataCombined, paramsList$numEventsRequired)
    PropFinal <- PropFinalOutcome$dataCombined
    PropFinalcensTime <- PropFinalOutcome$censTime
    PropFinalSS <- PropFinalOutcome$sampleSize
    coxmodel <- coxph(Surv(survival_time, status) ~ group, data = PropFinal)
    deltad <- as.numeric(exp(coef(coxmodel)))
    PropLRT <- survdiff(Surv(time, status) ~ group, data = PropFinal)
    Proppower <- (PropLRT$chisq > qchisq(0.95, 1) & deltad < 1 & PropOutcome == "Continue")
    PropCensTime <- ifelse(PropOutcome == "Stop1", Prop1censTime, ifelse(PropOutcome == "Stop2", Prop2censTime, PropFinalcensTime))
    PropSS <- ifelse(PropOutcome == "Stop1", Prop1SS, ifelse(PropOutcome == "Stop2", Prop2SS, PropFinalSS))
    
    
    # Return results
    list(NoIApower = NoIApower, NoIAcensTime = NoIAcensTime, NoIASS = NoIASS, 
         Wpower = Wpower, WCensTime = WCensTime, WSS = WSS, 
         OBFpower = OBFpower, OBFCensTime = OBFCensTime, OBFSS = OBFSS,
         Proppower = Proppower, PropCensTime = PropCensTime, PropSS = PropSS)
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
  ProppowerVec <- sapply(results, function(result) result$Proppower)
  PropCensTimeVec <- sapply(results, function(result) result$PropCensTime)
  PropSSVec <- sapply(results, function(result) result$PropSS)
  
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
  mean_Proppower <- mean(ProppowerVec)
  mean_PropCensTime <- mean(PropCensTimeVec)
  mean_PropSS <- mean(PropSSVec)
  
  
  if (i==1){
    outcomeDF <- data.frame(NoIAPower = mean_NoIApower, NoIADuration = mean_NoIAcensTime, NoIASS = mean_NoIASS,
                            WieandPower = mean_Wpower, WieandDuration = mean_WCensTime, WieandSS = mean_WSS,
                            OBFPower = mean_OBFpower, OBFDuration = mean_OBFCensTime, OBFSS = mean_OBFSS,
                            PropPower = mean_Proppower, PropDuration = mean_PropCensTime, PropSS = mean_PropSS)
    
    colnames(outcomeDF) <- c("No IA Power", "No IA Duration", "No IA SS", 
                             "Wieand Power", "Wieand Duration", "Wieand SS",
                             "OBF Power", "OBF Duration", "OBF SS",
                             "Prop Power", "Prop Duration", "Prop SS")
  } else {
    
    tempVec <- c(mean_NoIApower, mean_NoIAcensTime, mean_NoIASS,
                 mean_Wpower, mean_WCensTime, mean_WSS,
                 mean_OBFpower, mean_OBFCensTime, mean_OBFSS,
                 mean_Proppower, mean_PropCensTime, mean_PropSS)
    
    outcomeDF <- rbind(outcomeDF, tempVec)
  }
  
}

outcomeDF



par(mfrow = c(3,4))

for (i in 1:length(ScenarioList)){
  plot(outcomeDF$`No IA Power`[i], outcomeDF$`No IA Duration`[i], pch = 19, col = "red", ylim= c(10, 60), xlab = "Power", ylab = "Duration", xlim = c(0, 1))
  points(outcomeDF$`Wieand Power`[i], outcomeDF$`Wieand Duration`[i], pch = 19, col = "blue")
  points(outcomeDF$`OBF Power`[i], outcomeDF$`OBF Duration`[i], pch = 19, col = "yellow")
  points(outcomeDF$`Prop Power`[i], outcomeDF$`Prop Duration`[i], pch = 19, col = "green")
  #legend("bottomright", legend = c("No IA", "Wieand", "OBF", "Proposed"), col = c("red", "blue", "yellow", "green"), pch = 19)
}

par(mfrow = c(3,4))

for (i in 1:length(ScenarioList)){
  plot(outcomeDF$`No IA Power`[i], outcomeDF$`No IA SS`[i], pch = 19, col = "red", ylim= c(300, 700), xlab = "Power", ylab = "Duration", xlim = c(0, 1))
  points(outcomeDF$`OBF Power`[i], outcomeDF$`OBF SS`[i], pch = 19, col = "yellow")
  
  points(outcomeDF$`Prop Power`[i], outcomeDF$`Prop SS`[i], pch = 19, col = "green")
  points(outcomeDF$`Wieand Power`[i], outcomeDF$`Wieand SS`[i], pch = 19, col = "blue")
  
  #legend("bottomright", legend = c("No IA", "Wieand", "OBF", "Proposed"), col = c("red", "blue", "yellow", "green"), pch = 19)
}

















                        