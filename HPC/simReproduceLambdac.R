library(survival)
library(dplyr)
library(foreach)
library(doParallel)
library(rjags)


ScenarioList <- list(
  A = list(
    HR1 = 0.75,
    T1 = 1000,
    HR2 = 0.75,
    T2 = 1000,
    recTime = 60
  ),
  B = list(
    HR1 = 1,
    T1 = 1000,
    HR2 = 1,
    T2 = 1000,
    recTime = 60
  ),
  C = list(
    HR1 = 1.3,
    T1 = 1000,
    HR2 = 1.3,
    T2 = 1000,
    recTime = 60
  ),
  D = list(
    HR1 = 1,
    T1 = 3,
    HR2 = 0.693,
    T2 = 1000,
    recTime = 60
  ),
  E = list(
    HR1 = 1,
    T1 = 6,
    HR2 = 0.62,
    T2 = 1000,
    recTime = 60
  ),
  F = list(
    HR1 = 1.3,
    T1 = 3,
    HR2 = 0.628,
    T2 = 1000,
    recTime = 60
  )
)

# ScenarioList <- list(
#   A = list(
#     HR1 = 0.75,
#     T1 = 1000,
#     HR2 = 0.75,
#     T2 = 1000,
#     recTime = 34
#   ),
#   B = list(
#     HR1 = 1,
#     T1 = 1000,
#     HR2 = 1,
#     T2 = 1000,
#     recTime = 34
#   ),
#   C = list(
#     HR1 = 1.3,
#     T1 = 1000,
#     HR2 = 1.3,
#     T2 = 1000,
#     recTime = 34
#   ),
#   D = list(
#     HR1 = 1,
#     T1 = 3,
#     HR2 = 0.693,
#     T2 = 1000,
#     recTime = 34
#   ),
#   E = list(
#     HR1 = 1,
#     T1 = 6,
#     HR2 = 0.62,
#     T2 = 1000,
#     recTime = 34
#   ),
#   F = list(
#     HR1 = 1.3,
#     T1 = 3,
#     HR2 = 0.628,
#     T2 = 1000,
#     recTime = 34
#   ),
#   G = list(
#     HR1 = 0.75,
#     T1 = 1000,
#     HR2 = 0.75,
#     T2 = 1000,
#     recTime = 12
#   ),
#   H = list(
#     HR1 = 1,
#     T1 = 1000,
#     HR2 = 1,
#     T2 = 1000,
#     recTime = 12
#   ),
#   I = list(
#     HR1 = 1.3,
#     T1 = 1000,
#     HR2 = 1.3,
#     T2 = 1000,
#     recTime = 12
#   ),
#   J = list(
#     HR1 = 1,
#     T1 = 3,
#     HR2 = 0.693,
#     T2 = 1000,
#     recTime = 12
#   ),
#   K = list(
#     HR1 = 1,
#     T1 = 6,
#     HR2 = 0.62,
#     T2 = 1000,
#     recTime = 12
#   ),
#   L = list(
#     HR1 = 1.3,
#     T1 = 3,
#     HR2 = 0.628,
#     T2 = 1000,
#     recTime = 12
#   )
# )

# HR1Vec <- c(0.6, 0.75, 0.9, 1, 1.3)
# T1Vec <- c(0, 3, 6, 9)
# HR2Vec <- c(0.75, 1, 1.3)
# recTimeVec <- seq(0, 40, by=10)

# 
# HR1Vec <- 1.2
# T1Vec <- 0
# HR2Vec <-  1.2
# recTimeVec <- seq(0, 30, by=5)

# count <- 1
# 
# n <- length(HR1Vec)*length(T1Vec)*length(HR2Vec)*length(recTimeVec)  # Number of elements you want in the list
# ScenarioList <- vector("list", length = n)
# 
# 
# for (i in 1:length(HR1Vec)){
#   for (j in 1:length(T1Vec)){
#     for (k in 1:length(HR2Vec)){
#       for (l in 1:length(recTimeVec)){
#         ScenarioList[[count]]$HR1 = HR1Vec[i]
#         ScenarioList[[count]]$T1 = T1Vec[j]
#         ScenarioList[[count]]$HR2 = HR2Vec[k]
#         ScenarioList[[count]]$T2 = 1000
#         ScenarioList[[count]]$recTime = recTimeVec[l]
#         count <- count + 1
#       }
#     }
#   }
# }

# ScenarioList <- list(
#   A = list(
#     HR1 = 1,
#     T1 = 1000,
#     HR2 = 1,
#     T2 = 1000,
#     recTime = 12
#   ),
#   B = list(
#     HR1 = 1.1,
#     T1 = 1000,
#     HR2 = 1.1,
#     T2 = 1000,
#     recTime = 12
#   ),
#   C = list(
#     HR1 = 0.75,
#     T1 = 1000,
#     HR2 = 0.75,
#     T2 = 1000,
#     recTime = 12
#   ),
#   D = list(
#     HR1 = 1,
#     T1 = 3,
#     HR2 = 0.75,
#     T2 = 1000,
#     recTime = 12
#   )
# )


paramsList <- list(
  numPatients = 340,
  lambdac = log(2)/10,
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

NoIAFunc <- function(dataset){
  NoIAOutcome <- censFunc(dataset, paramsList$numEventsRequired)
  NoIA <- NoIAOutcome$dataCombined
  NoIAcensTime <- NoIAOutcome$censTime
  NoIASS <- NoIAOutcome$sampleSize
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = NoIA)
  deltad <- as.numeric(exp(coef(coxmodel)))
  NoIALRT <- survdiff(Surv(time, status) ~ group, data = NoIA)
  NoIApower <- (NoIALRT$chisq > qchisq(0.95, 1) & deltad < 1)
  return(list(NoIApower = NoIApower, NoIAcensTime = NoIAcensTime, NoIASS = NoIASS))
}

WieandFunc <- function(dataset){
  WieandOutcome <- "Continue"
  W1Outcome <- censFunc(dataset, paramsList$numEventsRequired * 0.5)
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
  return(list(Wpower = Wpower, WCensTime = WCensTime, WSS = WSS))
}

OBFFunc <- function(dataset){
  OBFOutcome <- "Continue"
  OBF1Outcome <- censFunc(dataset, ceiling(paramsList$numEventsRequired/3))
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
  return(list(OBFpower = OBFpower, OBFCensTime = OBFCensTime, OBFSS = OBFSS))
}

PropFunc <- function(dataset, monthsDelay, propEvents){
  
  sortedPseudoTimes <- sort(dataset$pseudo_time)
  propVec <- numeric(paramsList$numEventsRequired)
  cutoff_index <- ceiling(paramsList$numEventsRequired/2)
  
  for (k in 1:paramsList$numEventsRequired){
    censTime <- sortedPseudoTimes[k]
    
    newDataset <- dataset
    newDataset$status <- as.integer(newDataset$pseudo_time <= censTime)
    newDataset <- newDataset[newDataset$recTime <= censTime, ]
    newDataset$survival_time <- ifelse(newDataset$status == 1, newDataset$time, censTime - newDataset$recTime)
    propVec[k] <- mean(newDataset[newDataset$status==1, ]$survival_time>monthsDelay)
    
    # status <- as.integer(dataset$pseudo_time <= censTime)
    # filteredData <- dataset[dataset$recTime <= censTime, ]
    # survival_time <- ifelse(status == 1, filteredData$time, censTime - filteredData$recTime)
    #propVec[k] <- mean(survival_time[status==1]>monthsDelay)
  }
  
  #propVec[1:(cutoff_index-1)] <- 0
  
  if (sum(propVec>propEvents)==0){
    PropOutcome <- "Continue"
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
  } else{
    
    Stop1 <- max(paramsList$numEventsRequired*0.5, which(propVec>propEvents)[1])
    Stop2 <- max(paramsList$numEventsRequired*0.75, which(propVec>propEvents)[1])
    
    PropOutcome <- "Continue"
    Prop1Outcome <- censFunc(dataset, Stop1)
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
    
  }
  
  
  return(list(Proppower = Proppower, PropCensTime = PropCensTime, PropSS = PropSS))
  
}

for (i in 1:length(ScenarioList)){
  
  
  # Set the number of CPU cores you want to use
  num_cores <- 32 # Change this to the number of cores you want to use
  
  # Register parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  results <- foreach(j = 1:paramsList$NSims, .packages = c("survival", "dplyr", "rjags")) %dopar% {
    
    dataCombined <- generateData(paramsList$lambdac, ScenarioList[[i]]$HR1, 
                                 ScenarioList[[i]]$T1, ScenarioList[[i]]$HR2, ScenarioList[[i]]$T2, paramsList$numPatients, ScenarioList[[i]]$recTime)
    
    # Do it with no interim analysis first
    NoIAOutcome <- NoIAFunc(dataCombined)
    NoIApower <-  NoIAOutcome$NoIApower
    NoIAcensTime <- NoIAOutcome$NoIAcensTime
    NoIASS <- NoIAOutcome$NoIASS
    
    # Wieand rule
    WieandOutcome <- WieandFunc(dataCombined)
    Wpower <-  WieandOutcome$Wpower
    WCensTime <- WieandOutcome$WCensTime
    WSS <- WieandOutcome$WSS
    
    
    #OBF rule
    OBFOutcome <- OBFFunc(dataCombined)
    OBFpower <-  OBFOutcome$OBFpower
    OBFCensTime <- OBFOutcome$OBFCensTime
    OBFSS <- OBFOutcome$OBFSS
    
    
    PropOutcome <- PropFunc(dataCombined, 3, 2/3)
    Proppower <- PropOutcome$Proppower
    PropCensTime <- PropOutcome$PropCensTime
    PropSS <- PropOutcome$PropSS
    
    # Return results
    list(NoIApower = NoIApower, NoIAcensTime = NoIAcensTime, NoIASS = NoIASS, 
         Wpower = Wpower, WCensTime = WCensTime, WSS = WSS, 
         OBFpower = OBFpower, OBFCensTime = OBFCensTime, OBFSS = OBFSS,
         Proppower = Proppower, PropCensTime = PropCensTime, PropSS = PropSS)
  }
  
  # Extract results
  NoIApowerVec <- mean(sapply(results, function(result) result$NoIApower))
  NoIAcensTimeVec <- mean(sapply(results, function(result) result$NoIAcensTime))
  NoIASSVec <- mean(sapply(results, function(result) result$NoIASS))
  WpowerVec <- mean(sapply(results, function(result) result$Wpower))
  WCensTimeVec <- mean(sapply(results, function(result) result$WCensTime))
  WSSVec <- mean(sapply(results, function(result) result$WSS))
  OBFpowerVec <- mean(sapply(results, function(result) result$OBFpower))
  OBFCensTimeVec <- mean(sapply(results, function(result) result$OBFCensTime))
  OBFSSVec <- mean(sapply(results, function(result) result$OBFSS))
  ProppowerVec <- mean(sapply(results, function(result) result$Proppower))
  PropCensTimeVec <- mean(sapply(results, function(result) result$PropCensTime))
  PropSSVec <- mean(sapply(results, function(result) result$PropSS))
  
  # Clean up parallel resources
  stopCluster(cl)
  
  if (i==1){
    outcomeDF <- data.frame(NoIApowerVec = NoIApowerVec, NoIAcensTimeVec = NoIAcensTimeVec, NoIASSVec = NoIASSVec,
                            WpowerVec = WpowerVec, WCensTimeVec = WCensTimeVec, WSSVec = WSSVec,
                            OBFpowerVec = OBFpowerVec, OBFCensTimeVec = OBFCensTimeVec, OBFSSVec = OBFSSVec,
                            ProppowerVec = ProppowerVec, PropCensTimeVec = PropCensTimeVec, PropSSVec = PropSSVec
    )
    
    colnames(outcomeDF) <- c("No IA Power", "No IA Duration", "No IA SS", 
                             "Wieand Power", "Wieand Duration", "Wieand SS",
                             "OBF Power", "OBF Duration", "OBF SS",
                             "Prop Power", "Prop Duration", "Prop SS"
    )
  } else {
    
    tempVec <- c(NoIApowerVec, NoIAcensTimeVec, NoIASSVec,
                 WpowerVec, WCensTimeVec, WSSVec, 
                 OBFpowerVec, OBFCensTimeVec, OBFSSVec, 
                 ProppowerVec, PropCensTimeVec, PropSSVec)
    
    outcomeDF <- rbind(outcomeDF, tempVec)
  }
  
}


HR1Vec <- rep(NA, length(ScenarioList))
T1Vec <- rep(NA, length(ScenarioList))
HR2Vec <- rep(NA, length(ScenarioList))
recTimeVec <- rep(NA, length(ScenarioList))
# NoIAPowerRank <- rep(NA, length(ScenarioList))
# WieandPowerRank <- rep(NA, length(ScenarioList))
# OBFPowerRank <- rep(NA, length(ScenarioList))
# PropPowerRank <- rep(NA, length(ScenarioList))
# NoIADurationRank <- rep(NA, length(ScenarioList))
# WieandDurationRank <- rep(NA, length(ScenarioList))
# OBFDurationRank <- rep(NA, length(ScenarioList))
# PropDurationRank <- rep(NA, length(ScenarioList))
# NoIASSRank <- rep(NA, length(ScenarioList))
# WieandSSRank <- rep(NA, length(ScenarioList))
# OBFSSRank <- rep(NA, length(ScenarioList))
# PropSSRank <- rep(NA, length(ScenarioList))


for (i in 1:length(ScenarioList)){
  HR1Vec[i] <- ScenarioList[[i]]$HR1
  T1Vec[i] <- ScenarioList[[i]]$T1
  HR2Vec[i] <- ScenarioList[[i]]$HR2
  recTimeVec[i] <- ScenarioList[[i]]$recTime
  
  # PowerRank <- rank(-c(outcomeDF[i,]$`No IA Power`, outcomeDF[i,]$`Wieand Power`, outcomeDF[i,]$`OBF Power`, outcomeDF[i,]$`Prop Power`))
  # NoIAPowerRank[i] <- PowerRank[1]
  # WieandPowerRank[i] <- PowerRank[2]
  # OBFPowerRank[i] <- PowerRank[3]
  # PropPowerRank[i] <- PowerRank[4]
  # 
  # DurationRank <- rank(c(outcomeDF[i,]$`No IA Duration`, outcomeDF[i,]$`Wieand Duration`, outcomeDF[i,]$`OBF Duration`, outcomeDF[i,]$`Prop Duration`))
  # NoIADurationRank[i] <- DurationRank[1]
  # WieandDurationRank[i] <- DurationRank[2]
  # OBFDurationRank[i] <- DurationRank[3]
  # PropDurationRank[i] <- DurationRank[4]
  # 
  # SSRank <- rank(c(outcomeDF[i,]$`No IA SS`, outcomeDF[i,]$`Wieand SS`, outcomeDF[i,]$`OBF SS`, outcomeDF[i,]$`Prop SS`))
  # NoIASSRank[i] <- SSRank[1]
  # WieandSSRank[i] <- SSRank[2]
  # OBFSSRank[i] <- SSRank[3]
  # PropSSRank[i] <- SSRank[4]
  
}

outcomeDF <- cbind(outcomeDF, HR1Vec, T1Vec, HR2Vec, recTimeVec)
# outcomeDF <- cbind(outcomeDF, NoIAPowerRank, WieandPowerRank, OBFPowerRank, PropPowerRank)
# outcomeDF <- cbind(outcomeDF, NoIADurationRank, WieandDurationRank, OBFDurationRank, PropDurationRank)
# outcomeDF <- cbind(outcomeDF, NoIASSRank, WieandSSRank, OBFSSRank, PropSSRank)

outcomeDF




