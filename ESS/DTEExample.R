#***************************************************************************
#     Determining the Effective Sample Size of a Parametric Prior
#                                 by
#      Satoshi Morita(1), Peter F. Thall(2), and Peter Muller(2)
#       (1) Department of Epidemiology and Health Care Research,
#           Kyoto University Graduate School of Medicine, Kyoto, Japan
#       (2) Department of Biostatistics and Applied Mathematics,
#           University of Texas, M.D. Anderson Cancer Center,
#           1515 Holcombe Boulevard, Houston, Texas 77030, U.S.A.
#                              March XX, 2007
#***************************************************************************
#        R Program for Detemining the Effective Sample Sizes (ESSs)
#                     of the Prior in Thall and Lee (2003)
#***************************************************************************

library(purrr)
library(dplyr)

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

logistic <- function(M,d.ini,d.end)
{
  # Let p(theta_) be the prior on the parameter vector theta_ =
  #(theta.1,...,theta.d), where d denotes dim(theta_).
  # Set M being a positive integer chosen so that, initially, it is reasonable
  #to assume m, an ESS of p(theta_), <= M.
  # If M is not sufficiently large, 'NA' returns as a result of the
  # computations.
  # For Case 1, set both d.ini and d.end at 1.
  # For Case 2, set d.ini and d.end at 1 and d, respectively.
  # For Case 3, where theta = (theta_1,...,theta_K) is partitioned into K
  #  subvectors,
  # run the program to determine each mk, a prior ESS of each subvector
  #theta_k for k=1,...,K.
  # Set M, assuming mk being less than or equal to M. Set d.ini and d.end at
  # k.ini and k.end, respectively,
  # where (k.ini, ..., k.end) denoting the set of indices of the elements of
  # theta_k.
  
  Mrep <- M+1
  T    <- 1000
  c    <- 10000
  DqYMrep.out <- numeric(Mrep)
  trialTime <- 20
  
  # Set-ups for the study design of the example.
  d      <- 3  # d = dim(theta_)
  
  # Modification is required for the specifications in STEPs 1 and 2.
  # STEP 1
  # 1-1. Specify the prior p(theta_)
  # theta.1 ~ N(a1,b1) for intercept, theta.2 ~ N(a2,b2) for dose effect
  a1 <- 45.5
  b1 <- 780
  a2 <- 102
  b2 <- 136
  a3 <- 0.976
  b3 <- 0.25
  
  # 1-2. Specify theta_bar, the prior mean under p(theta_)
  # Refer to suitable textbooks such as Gelman et al. (2004)
  Etheta1 <- a1/b1
  Etheta2 <- a2/b2
  Etheta3 <- a3/b3
  
  # 1-3. Specify the hyperparamters of the epsilon-information prior
  #q0(theta_)
  # Refer to Table 1 of the paper */
  a1.0 <- a1/c
  b1.0 <- b1/c
  a2.0 <- a2/c
  b2.0 <- b2/c
  a3.0 <- a3/c
  b3.0 <- b3/c
  
  # STEP 2
  # 2-1. Compute the infomation matrix of p(theta_)
  Dp      <- c((a1-1)/Etheta1^2, (a2-1)/Etheta2^2,
               (a3-1)/Etheta3^2) # Specify Dp,j(theta_) for
  #j=1,...,d.
  Dp.plus <- sum(Dp[d.ini:d.end])
  # 2-2. Compute the expected information matrix of qm(theta_|data)
  Dq0      <- c((a1.0-1)/Etheta1^2, (a2.0-1)/Etheta2^2,
                (a3.0-1)/Etheta3^2)   # Specify - 2nd derivative of the
  #log {q0(theta_)}
  # regarding theta.j for j=1,...,d..
  Dq0.plus <- sum(Dq0[d.ini:d.end])
  for (t in 1:T) { # Simulate Monte Carlo samples
    DqYm.out <- numeric(M)
    DqY <- numeric(d)
    
    
    dataCombined <- generateData(Etheta1, 1, Etheta3, Etheta2, 340, 12)
    
    dataCombined <- censFunc(dataCombined, 512)$dataCombined
    
    dataCombined <- dataCombined[sample(1:nrow(dataCombined)), ]
    
    
    #We also need to make sure that we loop over the correct things here
    #Compare with the logistic regression case
    
    for (i in 1:M) {
      
      Dq.1 <- dataCombined[i,]$status/Etheta1^2
      
      if (dataCombined[i,]$group=="Treatment"&dataCombined[i,]$survival_time>Etheta3){
        Dq.2 <- dataCombined[i,]$status/Etheta2^2
      } else {
        Dq.2 <- 0
      }
      
      Dq.3 <- 0
      
      Dq   <- c(Dq.1, Dq.2, Dq.3)              
      
      DqY  <- DqY + Dq
      Dqm.plus    <- sum(DqY[d.ini:d.end])
      DqYm.out[i] <- Dqm.plus + Dq0.plus
    }
    DqYm.out    <- c(Dq0.plus, DqYm.out)
    DqYMrep.out <- rbind(DqYMrep.out,DqYm.out)
  }
  T1  <- T+1
  DqYMrep.out <- DqYMrep.out[c(2:T1),]
  Dqm.out     <- numeric(Mrep)
  Dqm <- 0
  for (i in 1:Mrep) {
    Dqm.out[i] <- mean(DqYMrep.out[,i])
  }
  # STEP 3
  # Compute delta(m,theta_bar,p,q0) to determine the ESS.
  D.m     <- Dqm.out - Dp.plus
  D.min.n <- which(abs(D.m) == min(abs(D.m)))
  D.min.v <- D.m[which(abs(D.m) == min(abs(D.m)))]
  {
    if (D.min.v < 0)       {
      D.min.v.nxt <- D.m[D.min.n+1]
      ESS <- D.min.n - 1 + (-D.min.v / (-D.min.v + D.min.v.nxt))
    }
    else if (D.min.v > 0)  {
      D.min.v.prv <- D.m[D.min.n-1]
      ESS <- D.min.n - 1 - (D.min.v / (D.min.v - D.min.v.prv))
    }
    else if (D.min.v == 0) {
      ESS <- D.min.n -1
    }
  }
  ESS
}

# For determining m
logistic(200,1,3)

