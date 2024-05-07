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
  Outcome <- 1
  if (deltad>observedHR){ Outcome <- 0}
  return(Outcome)
}

BPPFunc <- function(dataset, numPatients, numIAEvents, numFinalEvents, recTime, targetEff){
  
  BPPOutcome <- CensFunc(dataset, numIAEvents)
  
  dataCombined <- BPPOutcome$dataCombined
  
  dataCombined <- dataCombined[order(dataCombined$group), ]
  
  
  
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

update(model, n.iter=50, progress.bar = "none")
output=coda.samples(model=model, variable.names=c("HR", "bigT", "lambda2"), n.iter = 150, progress.bar = "none")


#The number of unenrolled patients in each group
cPatientsLeft <- numPatients - sum(dataCombined$group=="Control") 
tPatientsLeft <- numPatients - sum(dataCombined$group=="Treatment") 

#Extract realisations from the MCMC
HRoutput <- as.numeric(unlist(output[,1]))
bigToutput <- as.numeric(unlist(output[,2]))
lambda2output <- as.numeric(unlist(output[,3]))

#Calculate the proportion of values lower than the target effect
propEffect <- mean(HRoutput<targetEff)

BPPVec <- rep(NA, 50)

for (j in 1:50){
  
  #Sampling the recruitment times for the unenrolled patients
  unenrolledRecTimes <- runif(cPatientsLeft+tPatientsLeft, BPPOutcome$censTime, recTime)
  
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
  
  unenrolledData$pseudoTime <- unenrolledData$time + unenrolledData$recTime
  
  
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
  cCensoredData <- cCensoredData[,c(8, 2, 3, 9)]
  tBeforeDelay1 <- tBeforeDelay1[,c(8, 2, 3, 9)]
  tBeforeDelay2 <- tBeforeDelay2[,c(9, 2, 3, 10)]
  tAfterDelay <- tAfterDelay[,c(8, 2, 3, 9)]
  
  #Keeping the column names consistent
  colnames(cCensoredData) <- c("time", "group", "recTime", "pseudoTime")
  colnames(tBeforeDelay1) <- c("time", "group", "recTime", "pseudoTime")
  colnames(tBeforeDelay2) <- c("time", "group", "recTime", "pseudoTime")
  colnames(tAfterDelay) <- c("time", "group", "recTime", "pseudoTime")
  
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
  censTime1 <- sort(finalDataset$pseudoTime)[numFinalEvents]
  finalDataset$status <- finalDataset$pseudoTime <= censTime1
  finalDataset$status <- finalDataset$status*1
  finalDataset$enrolled <- finalDataset$recTime <= censTime1
  finalDataset <-  finalDataset[finalDataset$enrolled==T,]
  finalDataset$survival_time <- ifelse(finalDataset$pseudoTime>censTime1, censTime1  - finalDataset$recTime, finalDataset$time)
  
  #Testing for significance
  test <- survdiff(Surv(survival_time, status)~group, data = finalDataset)
  
  #kmfit <- survfit(Surv(survival_time, status)~group, data = finalDataset)
  #plot(kmfit, col = c("blue", "red"), xlim=c(0,50))
  
  #Making sure the significance is in the correct direction
  coxmodel <- coxph(Surv(survival_time, status)~group, data = finalDataset)
  deltad <- as.numeric(exp(coef(coxmodel)))
  
  
  BPPVec[j] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
  
}

return(list(BPP = mean(BPPVec), BPPSS = nrow(dataCombined), BPPCensTime = BPPOutcome$censTime, propEffect = propEffect))

}

proposedRuleFunc <- function(dataCombined, numEventsRequired, monthsDelay, propEvents){
  
  
  sortedPseudoTimes <- sort(dataCombined$pseudoTime)
  propVec <- numeric(numEventsRequired)
  cutoff_index <- ceiling(numEventsRequired/2)
  
  for (k in cutoff_index:numEventsRequired){
    censTime <- sortedPseudoTimes[k]
    status <- as.integer(dataCombined$pseudoTime <= censTime)
    filteredData <- dataCombined[dataCombined$recTime <= censTime, ]
    survival_time <- ifelse(status == 1, filteredData$time, censTime - filteredData$recTime)
    propVec[k] <- mean(survival_time[status==1]>monthsDelay)
  }
  
  propVec[1:(cutoff_index-1)] <- 0
  
  Stop1 <- max(numEventsRequired*0.5, which.min(propVec<propEvents))
  Stop2 <- max(numEventsRequired*0.75, which.min(propVec<propEvents))
  
  
  Prop1Outcome <- 1
  Prop1Cens <- CensFunc(dataCombined, Stop1)
  Prop1 <- Prop1Cens$dataCombined
  Prop1censTime <- Prop1Cens$censTime
  Prop1SS <- Prop1Cens$SS
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = Prop1)
  deltad <- as.numeric(exp(coef(coxmodel)))
  if (deltad > 1) Prop1Outcome <- 0
  
  Prop2Outcome <- 1
  Prop2Cens <- CensFunc(dataCombined, Stop2)
  Prop2 <- Prop2Cens$dataCombined
  Prop2censTime <- Prop2Cens$censTime
  Prop2SS <- Prop2Cens$SS
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = Prop2)
  deltad <- as.numeric(exp(coef(coxmodel)))
  if (deltad > 1) Prop2Outcome <- 0
  
  PropFinalOutcome <- 0
  PropFinalCens <- CensFunc(dataCombined, numEventsRequired)
  PropFinal <- PropFinalCens$dataCombined
  PropFinalcensTime <- PropFinalCens$censTime
  PropFinalSS <- PropFinalCens$SS
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = PropFinal)
  deltad <- as.numeric(exp(coef(coxmodel)))
  PropLRT <- survdiff(Surv(time, status) ~ group, data = PropFinal)
  if (PropLRT$chisq > qchisq(0.95, 1) & deltad < 1) PropFinalOutcome <- 1
  
  return(list(Prop1Outcome = Prop1Outcome, Prop2Outcome = Prop2Outcome, PropFinalOutcome = PropFinalOutcome,
              Prop1SS = Prop1SS, Prop2SS = Prop2SS, PropFinalSS = PropFinalSS,
              Prop1censTime = Prop1censTime, Prop2censTime = Prop2censTime, PropFinalcensTime = PropFinalcensTime))
  
  
}

GSDFunc <- function(dataCombined, futBound, effBound){
  
  coxmodel <- coxph(formula = surv_obj ~ group, data = data)
  y <- coef(summary(coxmodel))[,4]
  
  
  
}

