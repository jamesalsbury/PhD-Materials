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
assvec <- rep(NA, 5e3)

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
NRep <- 5e3

assMat <- data.frame(matrix(NA, ncol = length(IFVec)+1, nrow = NRep))
durMat <- data.frame(matrix(NA, ncol = length(IFVec)+1, nrow = NRep))
SSMat <- data.frame(matrix(NA, ncol = length(IFVec)+1, nrow = NRep))

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
    
  }
  
  
}

colMeans(durMat)
colMeans(assMat)
colMeans(SSMat)



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

censFunc <- CensFunc(trial_data, numEvents*IF)




