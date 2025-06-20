library(survival)
library(dplyr)
library(foreach)
library(doParallel)
library(rjags)

# ScenarioList <- list(
#   A = list(
#     HR1 = 0.75,
#     T1 = 1000,
#     HR2 = 0.75,
#     T2 = 1000,
#     recTime = 34
# ),
#   B = list(
#     HR1 = 1,
#     T1 = 1000,
#     HR2 = 1,
#     T2 = 1000,
#     recTime = 34
# ),
#   C = list(
#     HR1 = 1.3,
#     T1 = 1000,
#     HR2 = 1.3,
#     T2 = 1000,
#     recTime = 34
# ),
#   D = list(
#     HR1 = 1,
#     T1 = 3,
#     HR2 = 0.693,
#     T2 = 1000,
#     recTime = 34
# ),
#   E = list(
#     HR1 = 1,
#     T1 = 6,
#     HR2 = 0.62,
#     T2 = 1000,
#     recTime = 34
# ),
#   F = list(
#     HR1 = 1.3,
#     T1 = 3,
#     HR2 = 0.628,
#     T2 = 1000,
#     recTime = 34
# ),
#   G = list(
#     HR1 = 0.75,
#     T1 = 1000,
#     HR2 = 0.75,
#     T2 = 1000,
#     recTime = 12
# ),
#   H = list(
#     HR1 = 1,
#     T1 = 1000,
#     HR2 = 1,
#    T2 = 1000,
#     recTime = 12
# ),
#   I = list(
#    HR1 = 1.3,
#   T1 = 1000,
#    HR2 = 1.3,
#    T2 = 1000,
#    recTime = 12
# ),
#   J = list(
#     HR1 = 1,
#     T1 = 3,
#     HR2 = 0.693,
#     T2 = 1000,
#     recTime = 12
# ),
#   K = list(
#     HR1 = 1,
#     T1 = 6,
#     HR2 = 0.62,
#     T2 = 1000,
#     recTime = 12
# ),
#   L = list(
#     HR1 = 1.3,
#     T1 = 3,
#     HR2 = 0.628,
#     T2 = 1000,
#     recTime = 12
# )
# )


ScenarioList <- list(
  A = list(
    HR1 = 0.75,
    T1 = 1000,
    HR2 = 0.75,
    T2 = 1000,
    recTime = 8
  ),
  B = list(
    HR1 = 1,
    T1 = 1000,
    HR2 = 1,
    T2 = 1000,
    recTime = 8
  ),
  C = list(
    HR1 = 1.3,
    T1 = 1000,
    HR2 = 1.3,
    T2 = 1000,
    recTime = 8
  ),
  D = list(
    HR1 = 1,
    T1 = 3,
    HR2 = 0.693,
    T2 = 1000,
    recTime = 8
  ),
  E = list(
    HR1 = 1,
    T1 = 6,
    HR2 = 0.62,
    T2 = 1000,
    recTime = 8
  ), 
  F = list(
    HR1 = 1.3,
     T1 = 3,
    HR2 = 0.628,
    T2 = 1000,
    recTime = 8
    )
)

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
  lambdac = log(2)/8,
  numEventsRequired = 512,
  NSims = 5e1
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
  
  for (k in cutoff_index:paramsList$numEventsRequired){
    censTime <- sortedPseudoTimes[k]
    status <- as.integer(dataset$pseudo_time <= censTime)
    filteredData <- dataset[dataset$recTime <= censTime, ]
    survival_time <- ifelse(status == 1, filteredData$time, censTime - filteredData$recTime)
    propVec[k] <- mean(survival_time[status==1]>monthsDelay)
  }
  
  propVec[1:(cutoff_index-1)] <- 0
  
  Stop1 <- max(paramsList$numEventsRequired*0.5, which.min(propVec<propEvents))
  Stop2 <- max(paramsList$numEventsRequired*0.75, which.min(propVec<propEvents))
  
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
  
  return(list(Proppower = Proppower, PropCensTime = PropCensTime, PropSS = PropSS))
  
}

BPPFunc <- function(dataset, numEvents, recTime){
  
  BPPOutcome <- censFunc(dataset, numEvents)
  
  #survmodel <- survfit(Surv(survival_time, status)~group, data = NoIAOutcome$dataCombined)
  
  #plot(survmodel, col = c("blue", "red"))
  
  dataCombined <- BPPOutcome$dataCombined
  
  
  #JAGS code which calculates posterior distributions
  
  modelstring="

data {
  for (j in 1:m){
    zeros[j] <- 0
  }
}

model {
  C <- 10000
  for (i in 1:n){
    zeros[i] ~ dpois(zeros.mean[i])
    zeros.mean[i] <-  -l[i] + C
    l[i] <- ifelse(datEvent[i]==1, log(lambda2)-(lambda2*datTimes[i]), -(lambda2*datTimes[i]))
  }
  for (i in (n+1):m){                                                                                                             
    zeros[i] ~ dpois(zeros.mean[i])
    zeros.mean[i] <-  -l[i] + C
    l[i] <- ifelse(datEvent[i]==1, ifelse(datTimes[i]<bigT, log(lambda2)-(lambda2*datTimes[i]), log(lambda1)-lambda1*(datTimes[i]-bigT)-(bigT*lambda2)), 
      ifelse(datTimes[i]<bigT, -(lambda2*datTimes[i]), -(lambda2*bigT)-lambda1*(datTimes[i]-bigT)))
  }
  
    lambda2 ~ dbeta(1, 1)T(0,)
    HR ~ dnorm(0.75, 1)T(0,)
    bigT ~ dnorm(3, 1)T(0,)
    lambda1 <- lambda2*HR
    
}
    "

model = jags.model(textConnection(modelstring), data = list(datTimes = dataCombined$survival_time, 
                                                            datEvent = dataCombined$status, n = sum(dataCombined$group=="Control"), 
                                                            m=nrow(dataCombined)), quiet = T) 

update(model, n.iter=100)
output=coda.samples(model=model, variable.names=c("HR", "bigT", "lambda2"), n.iter = 250)

#plot(output)

#The number of unenrolled patients in each group
cPatientsLeft <- paramsList$numPatients - sum(dataCombined$group=="Control") 
tPatientsLeft <- paramsList$numPatients - sum(dataCombined$group=="Treatment") 

BPPVec <- rep(NA, 200)

for (j in 1:200){
  
  #Sampling the recruitment times for the unenrolled patients
  unenrolledRecTimes <- runif(cPatientsLeft+tPatientsLeft, BPPOutcome$censTime, recTime)
  
  #Extract realisations from the MCMC
  HRoutput <- as.numeric(unlist(output[,1]))
  bigToutput <- as.numeric(unlist(output[,2]))
  lambda2output <- as.numeric(unlist(output[,3]))
  
  #Sample values from the MCMC output
  sampledHR <- sample(HRoutput, 1)
  sampledbigT <- sample(bigToutput, 1)
  sampledlambda2 <- sample(lambda2output, 1)
  sampledlambda1 <- sampledlambda2*sampledHR
  
  
  #For the unenrolled data, we can sample the remaining data according to the updated (sampled) parameters
  CP <- exp(-(sampledlambda2*sampledbigT))
  u <- runif(tPatientsLeft)
  
  unenrolledData <- data.frame(time = c(rexp(cPatientsLeft, rate = sampledlambda2), ifelse(u>CP, (-log(u))/sampledlambda2, (1/sampledlambda1)*(sampledbigT*sampledlambda1-log(u)-sampledbigT*sampledlambda2))), group = c(rep("Control", cPatientsLeft),
                                                                                                                                                                                                                          rep("Treatment", tPatientsLeft)), recTime = unenrolledRecTimes)
  
  unenrolledData$pseudo_time <- unenrolledData$time + unenrolledData$recTime
  
  
  #Extracting the observations that were censored at the IA  
  censoredData <- dataCombined[dataCombined$status==0,]
  
  #Number of censored observations in each group
  cCensored <- sum(censoredData$group=="Control")
  tCensored <- sum(censoredData$group=="Treatment")
  
  #Extracting the censored observations in the control group
  cCensoredData <- censoredData %>%
    filter(group=="Control")
  
  #Adding a exp(sampledlambda2) value to the censored value
  cCensoredData$finalsurvTime <- cCensoredData$survival_time + rexp(cCensored, rate = sampledlambda2)
  
  #Calculating the psuedo time
  cCensoredData$finalPsuedoTime <- cCensoredData$recTime + cCensoredData$finalsurvTime
  
  #Extacting the observations in the treatment group which may still be influenced by the delay (their observation time is smaller than the sampled delay time)
  tBeforeDelay <- censoredData %>%
    filter(group=="Treatment") %>%
    filter(survival_time < sampledbigT)
  
  #Extracting the observations in the treatment group which will not be influenced by the delay (their observation time is bigger than the sampled delay time)
  tAfterDelay <- censoredData %>%
    filter(group=="Treatment") %>%
    filter(survival_time > sampledbigT)
  
  #As these observations are still subject to a delay, we add on a Exp(lambda2) (lambdac) time
  tBeforeDelay$IASurv <- tBeforeDelay$survival_time + rexp(nrow(tBeforeDelay), rate = sampledlambda2)
  
  #Extracting the observations in which the survival time is smaller than the sampled delay time
  tBeforeDelay1 <- tBeforeDelay %>%
    filter(IASurv < sampledbigT)
  
  #Extracting the observations in which the survival time is bigger than the sampled delay time
  tBeforeDelay2 <- tBeforeDelay %>%
    filter(IASurv > sampledbigT)
  
  #For the observations in which the survival time is bigger, we sample a Exp(lambda1) and add it to the sampled delay time
  tBeforeDelay2$IASurv2 <- sampledbigT + rexp(nrow(tBeforeDelay2), rate = sampledlambda1)
  
  #For the observations not influenced by the delay, we sample a Exp(lambda1) time and add it to the current survival time
  tAfterDelay$IASurv <- tAfterDelay$survival_time + rexp(nrow(tAfterDelay), rate = sampledlambda1)
  
  #Calculate the pseudo time for all the data frames
  tBeforeDelay1$IApsuedoTime <- tBeforeDelay1$IASurv + tBeforeDelay1$recTime
  tBeforeDelay2$IApsuedoTime <- tBeforeDelay2$IASurv2 + tBeforeDelay2$recTime
  tAfterDelay$IApsuedoTime <- tAfterDelay$IASurv + tAfterDelay$recTime
  
  #Only keeping the columns of interest
  cCensoredData <- cCensoredData[,c(2:3, 7:8)]
  tBeforeDelay1 <- tBeforeDelay1[,c(2:3, 7:8)]
  tBeforeDelay2 <- tBeforeDelay2[,c(2:3, 8:9)]
  tAfterDelay <- tAfterDelay[,c(2:3, 7:8)]
  
  #Keeping the column names consistent
  colnames(cCensoredData) <- c("group", "recTime", "time", "pseudo_time")
  colnames(tBeforeDelay1) <- c("group", "recTime", "time", "pseudo_time")
  colnames(tBeforeDelay2) <- c("group", "recTime", "time", "pseudo_time")
  colnames(tAfterDelay) <- c("group", "recTime", "time", "pseudo_time")
  
  #Only keeping observations from the censored data set which are dead
  finalDataset <- dataCombined %>%
    filter(status==1)
  
  finalDataset <- finalDataset[,1:4]
  
  #Combining all the above data sets 
  finalDataset <- rbind(finalDataset, tBeforeDelay1)
  finalDataset <- rbind(finalDataset, tBeforeDelay2)
  finalDataset <- rbind(finalDataset, tAfterDelay)
  finalDataset <- rbind(finalDataset, unenrolledData)
  finalDataset <- rbind(finalDataset, cCensoredData)
  
  #Making sure the final data set is correct
  censTime1 <- sort(finalDataset$pseudo_time)[paramsList$numEventsRequired]
  finalDataset$status <- finalDataset$pseudo_time <= censTime1
  finalDataset$status <- finalDataset$status*1
  finalDataset$enrolled <- finalDataset$recTime <= censTime1
  finalDataset <-  finalDataset[finalDataset$enrolled==T,]
  finalDataset$survival_time <- ifelse(finalDataset$pseudo_time>censTime1, censTime1  - finalDataset$recTime, finalDataset$time)
  
  #Testing for significance
  test <- survdiff(Surv(survival_time, status)~group, data = finalDataset)
  
  #kmfit <- survfit(Surv(survival_time, status)~group, data = finalDataset)
  #plot(kmfit, col = c("blue", "red"), xlim=c(0,50))
  
  #Making sure the significance is in the correct direction
  coxmodel <- coxph(Surv(survival_time, status)~group, data = finalDataset)
  deltad <- as.numeric(exp(coef(coxmodel)))
  
  
  BPPVec[j] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
  
}

return(list(BPP = mean(BPPVec), BPPSS = nrow(dataCombined), BPPCensTime = BPPOutcome$censTime))
  
}

BPPList <- vector("list", length(ScenarioList))

for (i in 1:length(ScenarioList)){
  
  
  # Set the number of CPU cores you want to use
  num_cores <- 64 # Change this to the number of cores you want to use
  
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
    
    BPPOutcome <- "Continue"
    BPP1Outcome <- BPPFunc(dataCombined, paramsList$numEventsRequired*0.5, ScenarioList[[i]]$recTime)
    BPP1censTime <- BPP1Outcome$BPPCensTime
    BPP1SS <- BPP1Outcome$BPPSS
    if (BPP1Outcome$BPP < 0.25) BPPOutcome <- "Stop1"
    
    BPP2Outcome <- BPPFunc(dataCombined, paramsList$numEventsRequired*0.75, ScenarioList[[i]]$recTime)
    BPP2censTime <- BPP2Outcome$BPPCensTime
    BPP2SS <- BPP2Outcome$BPPSS
    if (BPP2Outcome$BPP < 0.4 & BPPOutcome=="Continue") BPPOutcome <- "Stop2"

    BPPFinalOutcome <- censFunc(dataCombined, paramsList$numEventsRequired)
    BPPFinal <- BPPFinalOutcome$dataCombined
    BPPFinalcensTime <- BPPFinalOutcome$censTime
    BPPFinalSS <- BPPFinalOutcome$sampleSize
    coxmodel <- coxph(Surv(survival_time, status) ~ group, data = BPPFinal)
    deltad <- as.numeric(exp(coef(coxmodel)))
    BPPLRT <- survdiff(Surv(time, status) ~ group, data = BPPFinal)
    BPPpower <- (BPPLRT$chisq > qchisq(0.95, 1) & deltad < 1 & BPPOutcome == "Continue")
    BPPCensTime <- ifelse(BPPOutcome == "Stop1", BPP1censTime, ifelse(BPPOutcome == "Stop2", BPP2censTime, BPPFinalcensTime))
    BPPSS <- ifelse(BPPOutcome == "Stop1", BPP1SS, ifelse(BPPOutcome == "Stop2", BPP2SS, BPPFinalSS))
    
    
    # Return results
    list(NoIApower = NoIApower, NoIAcensTime = NoIAcensTime, NoIASS = NoIASS, 
         Wpower = Wpower, WCensTime = WCensTime, WSS = WSS, 
         OBFpower = OBFpower, OBFCensTime = OBFCensTime, OBFSS = OBFSS,
         Proppower = Proppower, PropCensTime = PropCensTime, PropSS = PropSS,
         BPPpower = BPPpower, BPPCensTime = BPPCensTime, BPPSS = BPPSS,
         BPP1 = BPP1Outcome$BPP, BPP2 = BPP2Outcome$BPP)
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
  BPPpowerVec <- mean(sapply(results, function(result) result$BPPpower))
  BPPCensTimeVec <- mean(sapply(results, function(result) result$BPPCensTime))
  BPPSSVec <- mean(sapply(results, function(result) result$BPPSS))
  BPPPower <- sapply(results, function(result) result$BPPpower)
  BPP1Vec <- sapply(results, function(result) result$BPP1)
  BPP2Vec <- sapply(results, function(result) result$BPP2)
  
  BPPList[[i]] <- data.frame(power = BPPPower, BPP1 = BPP1Vec, BPP2 = BPP2Vec)
  
  # Clean up parallel resources
  stopCluster(cl)

  if (i==1){
    outcomeDF <- data.frame(NoIApowerVec = NoIApowerVec, NoIAcensTimeVec = NoIAcensTimeVec, NoIASSVec = NoIASSVec,
                            WpowerVec = WpowerVec, WCensTimeVec = WCensTimeVec, WSSVec = WSSVec,
                            OBFpowerVec = OBFpowerVec, OBFCensTimeVec = OBFCensTimeVec, OBFSSVec = OBFSSVec,
                            ProppowerVec = ProppowerVec, PropCensTimeVec = PropCensTimeVec, PropSSVec = PropSSVec,
                            BPPpowerVec = BPPpowerVec, BPPCensTimeVec = BPPCensTimeVec, BPPSSVec = BPPSSVec
                            )
    
    colnames(outcomeDF) <- c("No IA Power", "No IA Duration", "No IA SS", 
                             "Wieand Power", "Wieand Duration", "Wieand SS",
                             "OBF Power", "OBF Duration", "OBF SS",
                             "Prop Power", "Prop Duration", "Prop SS",
                             "BPP Power", "BPP Duration", "BPP SS"
                             )
  } else {
    
    tempVec <- c(NoIApowerVec, NoIAcensTimeVec, NoIASSVec,
                 WpowerVec, WCensTimeVec, WSSVec, 
                 OBFpowerVec, OBFCensTimeVec, OBFSSVec, 
                 ProppowerVec, PropCensTimeVec, PropSSVec, 
                 BPPpowerVec, BPPCensTimeVec, BPPSSVec) 
    
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




