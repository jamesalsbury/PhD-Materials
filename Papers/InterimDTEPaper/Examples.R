library(survival)
library(rpact)

#######################################
##First calculate assurance (no IA)
#######################################

n_c <- 250
n_t <- 250
control_rate <- log(2)/10
#control_rate <- 0.07238056
P_S <- 0.9
P_DTE <- 0.7
recruitmentTime <- 12
numEvents <- 450
assvec <- rep(NA, 5e1)

for (i in 1:length(assvec)){
  control_times <- rexp(n_c, control_rate)
  
  if (runif(1) < P_S) {
    HRStar <- rgamma(1, 29.6, 47.8)
    if (runif(1) < P_DTE) {
      bigT <- rgamma(1, 7.29, 1.76)
    } else {
      bigT <- 0
    }
  } else {
    HRStar <- 1
    bigT <- 0
  }
  
  CP <- exp(-control_rate*bigT)
  
  u <- runif(n_t)
  treatment_times <- ifelse(u > CP, -log(u)/control_rate,
                            (1/(control_rate*HRStar))*(-log(u)-control_rate*bigT+control_rate*HRStar*bigT))
  
  
  trial_data <- data.frame(time = c(control_times, treatment_times),
                           group = c(rep("Control", n_c), rep("Treatment", n_t)))                    
  
  
  trial_data$rec_time <- runif(n_c + n_t, 0, recruitmentTime)
  
  trial_data$pseudo_time <- trial_data$time + trial_data$rec_time
  
  trial_data <- trial_data[order(trial_data$pseudo_time), ]
  
  censTime <- trial_data$pseudo_time[numEvents]
  
  trial_data <- trial_data[censTime > trial_data$rec_time, ]
  
  trial_data$status <- trial_data$pseudo_time <= censTime
  
  trial_data$survival_time <- ifelse(trial_data$pseudo_time <= censTime, trial_data$time, censTime - trial_data$rec_time)
  
  test <- survdiff(Surv(survival_time, status)~group, data = trial_data)
  coxmodel <- coxph(Surv(survival_time, status)~group, data = trial_data)
  deltad <- as.numeric(exp(coef(coxmodel)))
  
  assvec[i] <- test$chisq > qchisq(0.95, 1) & deltad<1
  
}

mean(assvec)

#######################################
##Now look at one IA
#######################################

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

computeCox <- function(data, events) {
  censDF <- CensFunc(data, events)
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = censDF$dataCombined)
  SS <- censDF$SS
  Duration <- censDF$censTime
  ZScore <- -(coef(summary(coxmodel))[, 4])
  delta <- as.numeric(exp(coef(coxmodel)))
  list(SS = SS, Duration = Duration, ZScore = ZScore, delta = delta)
}

GSDOneIAFunc <- function(dataCombined, futBound, critValues, IF, numEvents) {
  
  # Compute first and final IA
  IA1 <- computeCox(dataCombined, numEvents * IF[1])
  IA2 <- computeCox(dataCombined, numEvents * IF[2])
  
  # Determine outcome based on ZScores
  Outcome <- ifelse(IA1$ZScore > critValues[1], "Efficacy", 
                    ifelse(IA1$ZScore < futBound, "Futility", 
                           ifelse(IA2$ZScore > critValues[2], "Successful", "Unsuccessful")))
  
  # Determine SS and Duration based on outcome
  SS <- ifelse(Outcome %in% c("Efficacy", "Futility"), IA1$SS, IA2$SS)
  Duration <- ifelse(Outcome %in% c("Efficacy", "Futility"), IA1$Duration, IA2$Duration)
  
  return(list(Outcome = Outcome, SS = SS, Duration = Duration, 
              IA1Time = IA1$Duration, delta1 = IA1$delta, delta2 = IA2$delta))
}

NoIAFunc <- function(dataCombined, numEvents){
  
  FinalAnalysis <- computeCox(dataCombined, numEvents)
  
  Outcome <- FinalAnalysis$ZScore > qnorm(0.975)
  
  SS <- FinalAnalysis$SS
  Duration <- FinalAnalysis$Duration
  
  
  return(list(Outcome = Outcome, SS = SS, Duration = Duration))
              
  
  
}


n_c <- 250
n_t <- 250
control_rate <- log(2)/10
#control_rate <- 0.07238056
P_S <- 0.9
P_DTE <- 0.7
recruitmentTime <- 12
numEvents <- 450
IFVec <- seq(0.2, 0.8, by=0.2)
NRep <- 1e4

assMat <- data.frame(matrix(NA, ncol = length(IFVec)+1, nrow = NRep))
durMat <- data.frame(matrix(NA, ncol = length(IFVec)+1, nrow = NRep))
SSMat <- data.frame(matrix(NA, ncol = length(IFVec)+1, nrow = NRep))
outcomeMat <- data.frame(matrix(NA, ncol = length(IFVec)+1, nrow = NRep))

for (i in 1:NRep){
  control_times <- rexp(n_c, control_rate)
  
  if (runif(1) < P_S) {
    HRStar <- rgamma(1, 29.6, 47.8)
    if (runif(1) < P_DTE) {
      bigT <- rgamma(1, 7.29, 1.76)
    } else {
      bigT <- 0
    }
  } else {
    HRStar <- 1
    bigT <- 0
  }
  
  CP <- exp(-control_rate*bigT)
  
  u <- runif(n_t)
  treatment_times <- ifelse(u > CP, -log(u)/control_rate,
                            (1/(control_rate*HRStar))*(-log(u)-control_rate*bigT+control_rate*HRStar*bigT))
  
  
  trial_data <- data.frame(time = c(control_times, treatment_times),
                           group = c(rep("Control", n_c), rep("Treatment", n_t)))                    
  
  
  trial_data$recTime <- runif(n_c + n_t, 0, recruitmentTime)
  
  trial_data$pseudoTime <- trial_data$time + trial_data$recTime
  
  trial_data <- trial_data[order(trial_data$pseudoTime), ]
  
  NoIA <- NoIAFunc(trial_data, numEvents)
  
  assMat[i,length(IFVec)+1] <- NoIA$Outcome
  durMat[i,length(IFVec)+1] <- NoIA$Duration
  SSMat[i,length(IFVec)+1] <- NoIA$SS
  outcomeMat[i,length(IFVec)+1] <- NoIA$Outcome
  
  #Do this for the first IA
  
  for (k in 1:length(IFVec)){
    
    design <- getDesignGroupSequential(typeOfDesign = "asUser",
                                       informationRates = c(IFVec[k], 1),
                                       userAlphaSpending = c(0.0125, 0.025),
                                       typeBetaSpending = "bsUser",
                                       userBetaSpending = c(0.05, 0.1))
    
    
    IAOne <- GSDOneIAFunc(trial_data, design$futilityBounds, design$criticalValues, c(IFVec[k], 1), numEvents) 
    
    assMat[i,k] <- ifelse(IAOne$Outcome %in% c("Efficacy", "Successful"), 1, 0)
    
    durMat[i,k] <- IAOne$Duration
    
    SSMat[i,k] <- IAOne$SS
    
    outcomeMat[i,k] <- IAOne$Outcome
    
  }
  
}

colMeans(durMat)
colMeans(assMat)
colMeans(SSMat)


mean(outcomeMat[outcomeMat[,1] == "Futility", ]$X5)





#######################################
##Now look at two IAs
#######################################

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

computeCox <- function(data, events) {
  censDF <- CensFunc(data, events)
  coxmodel <- coxph(Surv(survival_time, status) ~ group, data = censDF$dataCombined)
  SS <- censDF$SS
  Duration <- censDF$censTime
  ZScore <- -(coef(summary(coxmodel))[, 4])
  delta <- as.numeric(exp(coef(coxmodel)))
  list(SS = SS, Duration = Duration, ZScore = ZScore, delta = delta)
}

GSDTwoIAFunc <- function(dataCombined, futBound, critValues, IF, numEvents) {
  
  # Compute first and final IA
  IA1 <- computeCox(dataCombined, numEvents * IF[1])
  IA2 <- computeCox(dataCombined, numEvents * IF[2])
  IA3 <- computeCox(dataCombined, numEvents * IF[3])
  
  # Determine outcome based on ZScores
  Outcome <- ifelse(IA1$ZScore > critValues[1], "Efficacy1", 
                    ifelse(IA1$ZScore < futBound[1], "Futility1", 
                           ifelse(IA2$ZScore > critValues[2], "Efficacy2",
                                  ifelse(IA2$ZScore < futBound[2], "Futility2", 
                                         ifelse(IA3$ZScore > critValues[3], "Successful", "Unsuccessful")))))
  
  # Determine SS and Duration based on outcome
  SS <- ifelse(Outcome %in% c("Efficacy1", "Futility1"), IA1$SS, ifelse(Outcome %in% c("Efficacy2", "Futility2"), IA2$SS, IA3$SS))
  Duration <- ifelse(Outcome %in% c("Efficacy1", "Futility1"), IA1$Duration, ifelse(Outcome %in% c("Efficacy2", "Futility2"), IA2$Duration, IA3$Duration))
  
  return(list(Outcome = Outcome, SS = SS, Duration = Duration, IA1Time = IA1$Duration, IA2Time = IA2$Duration,
              delta1 = IA1$delta, delta2 = IA2$delta, delta3 = IA3$delta))
}

NoIAFunc <- function(dataCombined, numEvents){
  
  FinalAnalysis <- computeCox(dataCombined, numEvents)
  
  Outcome <- FinalAnalysis$ZScore > qnorm(0.975)
  
  SS <- FinalAnalysis$SS
  Duration <- FinalAnalysis$Duration
  
  
  return(list(Outcome = Outcome, SS = SS, Duration = Duration))
  
  
  
}


n_c <- 250
n_t <- 250
control_rate <- log(2)/10
#control_rate <- 0.07238056
P_S <- 0.9
P_DTE <- 0.7
recruitmentTime <- 12
numEvents <- 450
IFList <- list(A = list(IA1 = 0.2, IA2 = 0.3),
               B = list(IA1 = 0.2, IA2 = 0.5),
               C = list(IA1 = 0.2, IA2 = 0.7),
               D = list(IA1 = 0.2, IA2 = 0.9),
               E = list(IA1 = 0.4, IA2 = 0.5),
               F = list(IA1 = 0.4, IA2 = 0.7),
               G = list(IA1 = 0.4, IA2 = 0.9),
               H = list(IA1 = 0.6, IA2 = 0.7),
               I = list(IA1 = 0.6, IA2 = 0.9),
               J = list(IA1 = 0.8, IA2 = 0.9))
               
               
NRep <- 5e3

assMat <- data.frame(matrix(NA, ncol = length(IFList)+1, nrow = NRep))
durMat <- data.frame(matrix(NA, ncol = length(IFList)+1, nrow = NRep))
SSMat <- data.frame(matrix(NA, ncol = length(IFList)+1, nrow = NRep))

for (i in 1:NRep){
  control_times <- rexp(n_c, control_rate)
  
  if (runif(1) < P_S) {
    HRStar <- rgamma(1, 29.6, 47.8)
    if (runif(1) < P_DTE) {
      bigT <- rgamma(1, 7.29, 1.76)
    } else {
      bigT <- 0
    }
  } else {
    HRStar <- 1
    bigT <- 0
  }
  
  CP <- exp(-control_rate*bigT)
  
  u <- runif(n_t)
  treatment_times <- ifelse(u > CP, -log(u)/control_rate,
                            (1/(control_rate*HRStar))*(-log(u)-control_rate*bigT+control_rate*HRStar*bigT))
  
  
  trial_data <- data.frame(time = c(control_times, treatment_times),
                           group = c(rep("Control", n_c), rep("Treatment", n_t)))                    
  
  
  trial_data$recTime <- runif(n_c + n_t, 0, recruitmentTime)
  
  trial_data$pseudoTime <- trial_data$time + trial_data$recTime
  
  trial_data <- trial_data[order(trial_data$pseudoTime), ]
  
  NoIA <- NoIAFunc(trial_data, numEvents)
  
  assMat[i,length(IFList)+1] <- NoIA$Outcome
  durMat[i,length(IFList)+1] <- NoIA$Duration
  SSMat[i,length(IFList)+1] <- NoIA$SS
  
  #Do this for the first IA
  
  for (k in 1:length(IFList)){
    
    design <- getDesignGroupSequential(typeOfDesign = "asUser",
                                       informationRates = c(IFList[[k]][[1]], IFList[[k]][[2]], 1),
                                       userAlphaSpending = c(0.025*(1/3), 0.025*(2/3), 0.025),
                                       typeBetaSpending = "bsUser",
                                       userBetaSpending = c(0.1*(1/3), 0.1*(2/3), 0.1))
    
    
    IATwo <- GSDTwoIAFunc(trial_data, design$futilityBounds, design$criticalValues,  c(IFList[[k]][[1]], IFList[[k]][[2]], 1), numEvents) 
    
    assMat[i,k] <- ifelse(IATwo$Outcome %in% c("Efficacy1", "Efficacy2", "Successful"), 1, 0)
    
    durMat[i,k] <- IATwo$Duration
    
    SSMat[i,k] <- IATwo$SS
    
  }
  
  
}


#######################################
##Now we perform the BPP
#######################################


n_c <- 250
n_t <- 250
control_rate <- log(2)/10
#control_rate <- 0.07238056
P_S <- 0.9
P_DTE <- 0.7
recruitmentTime <- 12
numEvents <- 450
IF <- 0.8
Nrep <- 5e2
distParambigT <- paste0("bigT2 ~ dgamma(7.29, 1.76)")
distParamHR <- paste0("HR2 ~ dgamma(29.6, 47.8)")

# Number of cores to use (you can modify this based on your system)
num_cores <- parallel::detectCores() - 1  # Reserve 1 core for system tasks

# Create a cluster and register it
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Parallel loop using foreach
BPPHist <- foreach(i = 1:Nrep, .combine = 'c', .packages = c("rjags", "survival", "dplyr")) %dopar% {
  control_times <- rexp(n_c, control_rate)
  
  if (runif(1) < P_S) {
    HRStar <- rgamma(1, 29.6, 47.8)
    if (runif(1) < P_DTE) {
      bigT <- rgamma(1, 7.29, 1.76)
    } else {
      bigT <- 0
    }
  } else {
    HRStar <- 1
    bigT <- 0
  }
  
  CP <- exp(-control_rate*bigT)
  
  u <- runif(n_t)
  treatment_times <- ifelse(u > CP, -log(u)/control_rate,
                            (1/(control_rate*HRStar))*(-log(u)-control_rate*bigT+control_rate*HRStar*bigT))
  
  trial_data <- data.frame(time = c(control_times, treatment_times),
                           group = c(rep("Control", n_c), rep("Treatment", n_t)))                    
  
  trial_data$recTime <- runif(n_c + n_t, 0, recruitmentTime)
  trial_data$pseudoTime <- trial_data$time + trial_data$recTime
  trial_data <- trial_data[order(trial_data$pseudoTime), ]
  
  BPPOutcome <- CensFunc(trial_data, numEvents * IF)
  dataCombined <- BPPOutcome$dataCombined
  dataCombined <- dataCombined[order(dataCombined$group), ]
  
  # JAGS code for posterior distributions
  modelString <- paste0(
    "data {\n",
    "  for (j in 1:m){\n",
    "    zeros[j] <- 0\n",
    "  }\n",
    "}\n",
    "\n",
    "model {\n",
    "  C <- 10000\n",
    "  for (i in 1:n){\n",
    "    zeros[i] ~ dpois(zeros.mean[i])\n",
    "    zeros.mean[i] <-  -l[i] + C\n",
    "    l[i] <- ifelse(datEvent[i]==1, log(lambda2)-(lambda2*datTimes[i]), -(lambda2*datTimes[i]))\n",
    "  }\n",
    "  for (i in (n+1):m){\n",
    "    zeros[i] ~ dpois(zeros.mean[i])\n",
    "    zeros.mean[i] <-  -l[i] + C\n",
    "    l[i] <- ifelse(datEvent[i]==1, ifelse(datTimes[i]<bigT, log(lambda2)-(lambda2*datTimes[i]), log(lambda1)-lambda1*(datTimes[i]-bigT)-(bigT*lambda2)),\n",
    "      ifelse(datTimes[i]<bigT, -(lambda2*datTimes[i]), -(lambda2*bigT)-lambda1*(datTimes[i]-bigT)))\n",
    "  }\n",
    " \n",
    "    lambda2 ~ dbeta(1, 1)T(0,)\n",
    "   \n",
    "    mixT ~ dbern(1-P_S*P_DTE)\n",
    "    bigT <- mixT * bigT1 + (1-mixT) * bigT2\n",
    "    bigT1 ~ dnorm(0, 100)\n",
    "    ", distParambigT, "\n",
    "   \n",
    "    mixHR ~ dbern(1-P_S)\n",
    "    HR <- mixHR * HR1 + (1-mixHR) * HR2\n",
    "    HR1 ~ dnorm(1, 10000)\n",
    "    ", distParamHR, "\n",
    "   \n",
    "    lambda1 <- lambda2 * HR\n",
    "}"
  )
  
  model = jags.model(textConnection(modelString), data = list(datTimes = dataCombined$survival_time,
                                                              datEvent = dataCombined$status,
                                                              n = sum(dataCombined$group == "Control"),
                                                              m = nrow(dataCombined),
                                                              P_S = P_S,
                                                              P_DTE = P_DTE), quiet = TRUE)
  update(model, n.iter = 50, progress.bar = "none")
  output = coda.samples(model = model, variable.names = c("HR", "bigT", "lambda2"), n.iter = 1000, progress.bar = "none")
  
  #The number of unenrolled patients in each group
  cPatientsLeft <- n_c - sum(dataCombined$group=="Control")
  tPatientsLeft <- n_t - sum(dataCombined$group=="Treatment")
  
  #Extract realisations from the MCMC
  HRoutput <- as.numeric(unlist(output[,1]))
  bigToutput <- as.numeric(unlist(output[,2]))
  lambda2output <- as.numeric(unlist(output[,3]))
  
  
  # Further processing like before (sampling, testing, etc.)
  BPPVec <- rep(NA, 500)
  
  for (j in 1:500) {
    #Sampling the recruitment times for the unenrolled patients
    unenrolledRecTimes <- runif(cPatientsLeft+tPatientsLeft, BPPOutcome$censTime, recruitmentTime)
    
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
    censTime1 <- sort(finalDataset$pseudoTime)[numEvents]
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
    
    BPPVec[j] <- (test$chisq > qchisq(0.95, 1) & deltad < 1)
  }
  
  mean(BPPVec)
}

# Stop the cluster
stopCluster(cl)


vector_data <- c(1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.916, 0.232, 0.984, 1.000, 1.000, 1.000, 1.000, 0.284, 1.000, 0.984, 
                 1.000, 1.000, 0.000, 0.740, 0.000, 1.000, 1.000, 0.094, 1.000, 1.000, 0.998, 1.000, 0.034, 1.000, 1.000, 1.000, 
                 1.000, 0.462, 1.000, 1.000, 1.000, 0.998, 1.000, 1.000, 1.000, 1.000, 0.236, 1.000, 0.998, 1.000, 0.036, 1.000, 
                 0.116, 1.000, 1.000, 1.000, 1.000, 1.000, 0.000, 1.000, 0.016, 0.996, 1.000, 0.000, 1.000, 1.000, 0.000, 1.000, 
                 1.000, 1.000, 0.874, 0.354, 0.172, 1.000, 1.000, 0.956, 1.000, 1.000, 1.000, 0.994, 0.742, 1.000, 1.000, 0.030, 
                 0.000, 1.000, 1.000, 0.174, 1.000, 0.698, 0.358, 1.000, 0.984, 1.000, 1.000, 1.000, 0.000, 1.000, 0.568, 1.000, 
                 1.000, 1.000, 1.000, 0.748, 0.960, 1.000, 0.998, 1.000, 1.000, 1.000, 1.000, 0.302, 1.000, 1.000, 1.000, 1.000, 
                 1.000, 1.000, 1.000, 0.990, 1.000, 1.000, 0.006, 1.000, 1.000, 0.998, 0.100, 0.030, 1.000, 0.272, 1.000, 1.000, 
                 1.000, 0.996, 0.102, 0.024, 1.000, 1.000, 0.656, 1.000, 0.918, 1.000, 1.000, 0.000, 1.000, 0.374, 0.074, 1.000, 
                 0.000, 1.000, 0.000, 1.000, 1.000, 1.000, 0.028, 1.000, 1.000, 0.120, 1.000, 1.000, 1.000, 0.000, 0.068, 0.940, 
                 0.000, 1.000, 0.832, 0.068, 0.544, 1.000, 1.000, 1.000, 0.056, 1.000, 1.000, 1.000, 1.000, 0.000, 1.000, 1.000, 
                 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.726, 0.956, 0.990, 1.000, 1.000, 1.000, 1.000, 
                 1.000, 0.000, 0.684, 1.000, 1.000, 0.000, 0.338, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.998, 1.000, 0.332, 
                 1.000, 1.000, 1.000, 1.000, 0.960, 1.000, 1.000, 0.996, 0.998, 0.014, 1.000, 1.000, 0.000, 0.002, 1.000, 0.792, 
                 1.000, 1.000, 0.942, 0.378, 1.000, 1.000, 1.000, 0.916, 0.728, 1.000, 1.000, 1.000, 1.000, 0.736, 1.000, 0.282, 
                 0.458, 0.992, 1.000, 1.000, 1.000, 1.000, 1.000, 0.000, 1.000, 1.000, 1.000, 0.188, 1.000, 1.000, 1.000, 0.914, 
                 0.000, 1.000, 1.000, 0.006, 1.000, 1.000, 0.208, 1.000, 0.992, 0.000, 0.002, 0.948, 0.990, 1.000, 1.000, 0.734, 
                 1.000, 1.000, 1.000, 1.000, 0.992, 1.000, 1.000, 1.000, 0.268, 1.000, 1.000, 0.974, 1.000, 0.006, 1.000, 0.754, 
                 0.010, 0.524, 1.000, 0.000, 0.934, 0.000, 0.000, 1.000, 0.968, 1.000, 1.000, 1.000, 0.122, 1.000, 0.770, 1.000, 
                 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.024, 1.000, 0.864, 1.000, 0.958, 1.000, 1.000, 1.000, 1.000, 1.000, 
                 1.000, 1.000, 0.424, 1.000, 1.000, 0.984, 0.998, 1.000, 0.996, 0.648, 0.438, 1.000, 1.000, 0.440, 1.000, 1.000, 
                 0.956, 0.984, 0.812, 1.000, 0.818, 1.000, 1.000, 1.000, 0.984, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.000, 
                 1.000, 0.980, 1.000, 1.000, 0.984, 0.998, 0.998, 1.000, 1.000, 1.000, 0.876, 0.184, 0.000, 0.252, 1.000, 1.000, 
                 1.000, 0.044, 0.998, 0.038, 0.038, 0.002, 1.000, 1.000, 1.000, 0.998, 1.000, 1.000, 1.000, 0.000, 1.000, 0.148, 
                 1.000, 0.024, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.004, 1.000, 0.960, 1.000, 
                 1.000, 1.000, 1.000, 0.006, 0.000, 0.824, 0.000, 1.000, 1.000, 0.492, 0.000, 0.984, 0.998, 1.000, 1.000, 1.000, 
                 0.998, 1.000, 0.000, 1.000, 0.984, 1.000, 0.692, 0.000, 0.988, 0.998, 1.000, 1.000, 1.000, 0.812, 0.986, 0.000, 
                 1.000, 1.000, 0.206, 1.000, 1.000, 1.000, 1.000, 0.898, 0.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 
                 1.000, 0.440, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.698, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 
                 1.000, 1.000, 0.474, 1.000, 1.000, 0.374, 1.000, 1.000, 1.000, 0.986, 0.470, 0.998, 0.996, 1.000, 0.950, 0.032, 
                 0.976, 1.000, 1.000, 0.994, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.444, 1.000, 
                 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.960, 1.000, 0.980, 1.000, 1.000, 0.502, 0.420)



