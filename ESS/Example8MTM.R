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


#Etheta1
Etheta1Func <- function(Etheta1, Etheta2, Etheta3, Etheta4, Etheta5, Etheta6, X1, X2, Y){
  result <- -(X1^(2*Etheta2)*((X1^(2*Etheta2)*Y-X1^(2*Etheta2))*Etheta1^2+((2*X1^Etheta2*(X1^Etheta2*X2^Etheta4)^Etheta6*Y-2*X1^Etheta2*(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta5+(2*X1^Etheta2*X2^Etheta4*Y-2*X1^Etheta2*X2^Etheta4)*Etheta3+2*X1^Etheta2*Y)*Etheta1+((X1^Etheta2*X2^Etheta4)^(2*Etheta6)*Y-(X1^Etheta2*X2^Etheta4)^(2*Etheta6))*Etheta5^2+((2*X2^Etheta4*(X1^Etheta2*X2^Etheta4)^Etheta6*Y-2*X2^Etheta4*(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta3+2*(X1^Etheta2*X2^Etheta4)^Etheta6*Y)*Etheta5+(X2^(2*Etheta4)*Y-X2^(2*Etheta4))*Etheta3^2+2*X2^Etheta4*Y*Etheta3+Y))/((X1^Etheta2*Etheta1+(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3)^2*(X1^Etheta2*Etheta1+(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+1)^2)
  return(-result)
}


#Etheta2
Etheta2Func <- function(Etheta1, Etheta2, Etheta3, Etheta4, Etheta5, Etheta6, X1, X2, Y){
  result <- (log(X1)*(((X1^Etheta2*X2^Etheta4)^Etheta6*Y-(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta5+(X2^Etheta4*Y-X2^Etheta4)*Etheta3+(X1^Etheta2*Y-X1^Etheta2)*Etheta1+Y)*(log(X1)*(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*Etheta6^2+X1^Etheta2*log(X1)*Etheta1))/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1))+(log(X1)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*Etheta6+X1^Etheta2*Etheta1)*(Etheta5*(log(X1)*(X1^Etheta2*X2^Etheta4)^Etheta6*Y*Etheta6-log(X1)*(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta6)+(X1^Etheta2*log(X1)*Y-X1^Etheta2*log(X1))*Etheta1))/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1))-(log(X1)*(((X1^Etheta2*X2^Etheta4)^Etheta6*Y-(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta5+(X2^Etheta4*Y-X2^Etheta4)*Etheta3+(X1^Etheta2*Y-X1^Etheta2)*Etheta1+Y)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*Etheta6+X1^Etheta2*Etheta1)*(log(X1)*(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*Etheta6+X1^Etheta2*log(X1)*Etheta1))/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)^2*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1))-(log(X1)*(((X1^Etheta2*X2^Etheta4)^Etheta6*Y-(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta5+(X2^Etheta4*Y-X2^Etheta4)*Etheta3+(X1^Etheta2*Y-X1^Etheta2)*Etheta1+Y)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*Etheta6+X1^Etheta2*Etheta1)*(log(X1)*(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*Etheta6+X1^Etheta2*log(X1)*Etheta1))/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1)^2)
  return(-result)
}


#Etheta3
Etheta3Func <- function(Etheta1, Etheta2, Etheta3, Etheta4, Etheta5, Etheta6, X1, X2, Y){
  result <- -(X2^(2*Etheta4)*((X2^(2*Etheta4)*Y-X2^(2*Etheta4))*Etheta3^2+((2*X2^Etheta4*(X1^Etheta2*X2^Etheta4)^Etheta6*Y-2*X2^Etheta4*(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta5+(2*X1^Etheta2*X2^Etheta4*Y-2*X1^Etheta2*X2^Etheta4)*Etheta1+2*X2^Etheta4*Y)*Etheta3+((X1^Etheta2*X2^Etheta4)^(2*Etheta6)*Y-(X1^Etheta2*X2^Etheta4)^(2*Etheta6))*Etheta5^2+((2*X1^Etheta2*(X1^Etheta2*X2^Etheta4)^Etheta6*Y-2*X1^Etheta2*(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta1+2*(X1^Etheta2*X2^Etheta4)^Etheta6*Y)*Etheta5+(X1^(2*Etheta2)*Y-X1^(2*Etheta2))*Etheta1^2+2*X1^Etheta2*Y*Etheta1+Y))/((X2^Etheta4*Etheta3+(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X1^Etheta2*Etheta1)^2*(X2^Etheta4*Etheta3+(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X1^Etheta2*Etheta1+1)^2)
  return(-result)
}


#Etheta4
Etheta4Func <- function(Etheta1, Etheta2, Etheta3, Etheta4, Etheta5, Etheta6, X1, X2, Y){
  result <- (log(X2)*(((X1^Etheta2*X2^Etheta4)^Etheta6*Y-(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta5+(X2^Etheta4*Y-X2^Etheta4)*Etheta3+(X1^Etheta2*Y-X1^Etheta2)*Etheta1+Y)*(log(X2)*(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*Etheta6^2+X2^Etheta4*log(X2)*Etheta3))/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1))+(log(X2)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*Etheta6+X2^Etheta4*Etheta3)*(Etheta5*(log(X2)*(X1^Etheta2*X2^Etheta4)^Etheta6*Y*Etheta6-log(X2)*(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta6)+(X2^Etheta4*log(X2)*Y-X2^Etheta4*log(X2))*Etheta3))/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1))-(log(X2)*(((X1^Etheta2*X2^Etheta4)^Etheta6*Y-(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta5+(X2^Etheta4*Y-X2^Etheta4)*Etheta3+(X1^Etheta2*Y-X1^Etheta2)*Etheta1+Y)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*Etheta6+X2^Etheta4*Etheta3)*(log(X2)*(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*Etheta6+X2^Etheta4*log(X2)*Etheta3))/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)^2*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1))-(log(X2)*(((X1^Etheta2*X2^Etheta4)^Etheta6*Y-(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta5+(X2^Etheta4*Y-X2^Etheta4)*Etheta3+(X1^Etheta2*Y-X1^Etheta2)*Etheta1+Y)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*Etheta6+X2^Etheta4*Etheta3)*(log(X2)*(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*Etheta6+X2^Etheta4*log(X2)*Etheta3))/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1)^2)
  return(-result)
}


#Etheta5
Etheta5Func <- function(Etheta1, Etheta2, Etheta3, Etheta4, Etheta5, Etheta6, X1, X2, Y){
  result <- -((X1^Etheta2*X2^Etheta4)^(2*Etheta6)*(((X1^Etheta2*X2^Etheta4)^(2*Etheta6)*Y-(X1^Etheta2*X2^Etheta4)^(2*Etheta6))*Etheta5^2+((2*X2^Etheta4*(X1^Etheta2*X2^Etheta4)^Etheta6*Y-2*X2^Etheta4*(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta3+(2*X1^Etheta2*(X1^Etheta2*X2^Etheta4)^Etheta6*Y-2*X1^Etheta2*(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta1+2*(X1^Etheta2*X2^Etheta4)^Etheta6*Y)*Etheta5+(X2^(2*Etheta4)*Y-X2^(2*Etheta4))*Etheta3^2+((2*X1^Etheta2*X2^Etheta4*Y-2*X1^Etheta2*X2^Etheta4)*Etheta1+2*X2^Etheta4*Y)*Etheta3+(X1^(2*Etheta2)*Y-X1^(2*Etheta2))*Etheta1^2+2*X1^Etheta2*Y*Etheta1+Y))/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)^2*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1)^2)
  return(-result)
}


#Etheta6
Etheta6Func <- function(Etheta1, Etheta2, Etheta3, Etheta4, Etheta5, Etheta6, X1, X2, Y){
  result <- (log(X1^Etheta2*X2^Etheta4)*(X1^Etheta2*X2^Etheta4)^Etheta6*(log(X1^Etheta2*X2^Etheta4)*(X1^Etheta2*X2^Etheta4)^Etheta6*Y-log(X1^Etheta2*X2^Etheta4)*(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta5^2)/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1))+(log(X1^Etheta2*X2^Etheta4)^2*(X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5*(((X1^Etheta2*X2^Etheta4)^Etheta6*Y-(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta5+(X2^Etheta4*Y-X2^Etheta4)*Etheta3+(X1^Etheta2*Y-X1^Etheta2)*Etheta1+Y))/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1))-(log(X1^Etheta2*X2^Etheta4)^2*(X1^Etheta2*X2^Etheta4)^(2*Etheta6)*Etheta5^2*(((X1^Etheta2*X2^Etheta4)^Etheta6*Y-(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta5+(X2^Etheta4*Y-X2^Etheta4)*Etheta3+(X1^Etheta2*Y-X1^Etheta2)*Etheta1+Y))/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)^2*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1))-(log(X1^Etheta2*X2^Etheta4)^2*(X1^Etheta2*X2^Etheta4)^(2*Etheta6)*Etheta5^2*(((X1^Etheta2*X2^Etheta4)^Etheta6*Y-(X1^Etheta2*X2^Etheta4)^Etheta6)*Etheta5+(X2^Etheta4*Y-X2^Etheta4)*Etheta3+(X1^Etheta2*Y-X1^Etheta2)*Etheta1+Y))/(((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1)*((X1^Etheta2*X2^Etheta4)^Etheta6*Etheta5+X2^Etheta4*Etheta3+X1^Etheta2*Etheta1+1)^2)
  return(-result)
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
  T    <- 100
  c    <- 10000
  DqYMrep.out <- numeric(Mrep)
  d <- 6
  
  
  # Modification is required for the specifications in STEPs 1 and 2.
  # STEP 1
  # 1-1. Specify the prior p(theta_)
  # theta.i ~ Ga(ai,bi)
  a1 <- 1.74
  b1 <- 4.07
  a2 <- 10.24
  b2 <- 1.34
  a3 <- 2.32
  b3 <- 5.42
  a4 <- 15.24
  b4 <- 1.95
  a5 <- 0.33
  b5 <- 0.33
  a6 <- 0.0008
  b6 <- 0.0167
  # 1-2. Specify theta_bar, the prior mean under p(theta_)
  # Refer to suitable textbooks such as Gelman et al. (2004)
  #theta_bar in this Gamma case is a/b
  Etheta1 <- a1/b1
  Etheta2 <- a2/b2
  Etheta3 <- a3/b3
  Etheta4 <- a4/b4
  Etheta5 <- a5/b5
  Etheta6 <- a6/b6
  # 1-3. Specify the hyperparamters of the epsilon-information prior
  #q0(theta_)
  # Refer to Table 1 of the paper */
  a1.0 <- a1/c
  b1.0 <- b1/c
  a2.0 <- a2/c
  b2.0 <- b2/c
  a3.0 <- a3/c
  b3.0 <- b3/c
  a4.0 <- a4/c
  b4.0 <- b4/c
  a5.0 <- a5/c
  b5.0 <- b5/c
  a6.0 <- a6/c
  b6.0 <- b6/c
  
  # STEP 2
  # 2-1. Compute the infomation matrix of p(theta_)
  Dp      <- c((a1-1)/Etheta1^2, (a2-1)/Etheta2^2,
               (a3-1)/Etheta3^2, (a4-1)/Etheta4^2,
               (a5-1)/Etheta5^2, (a6-1)/Etheta6^2) # Specify Dp,j(theta_) for
  #j=1,...,d.
  Dp.plus <- sum(Dp[d.ini:d.end])
  # 2-2. Compute the expected information matrix of qm(theta_|data)
  Dq0      <- c((a1.0-1)/Etheta1^2, (a2.0-1)/Etheta2^2,
                (a3.0-1)/Etheta3^2, (a4.0-1)/Etheta4^2,
                (a5.0-1)/Etheta5^2, (a6.0-1)/Etheta6^2)   # Specify - 2nd derivative of the
  #log {q0(theta_)}
  # regarding theta.j for j=1,...,d..
  Dq0.plus <- sum(Dq0[d.ini:d.end])
  for (t in 1:T) { # Simulate Monte Carlo samples
    DqYm.out <- numeric(M)
    DqY <- numeric(d)
    for (i in 1:M) {
      u <- runif(1)
      if (u < 1/3){
        X1 <- 0.001
        X2 <- runif(1, 0, 1)
      } else if (u > 2/3){
        X1 <- runif(1, 0, 1)
        X2 <- 0.001
      } else {
        X1 <- runif(1, 0, 1)
        X2 <- runif(1, 0, 1)
      }
      
      
      probY <- (Etheta1*X1^Etheta2+Etheta3*X2^Etheta4+Etheta5*(X1^Etheta2*X2^Etheta4)^Etheta6)/
        (1+Etheta1*X1^Etheta2+Etheta3*X2^Etheta4+Etheta5*(X1^Etheta2*X2^Etheta4)^Etheta6)
      
      Y <- rbinom(n = 1, size = 1, prob = probY)
      
      # Specify the model employed
      # in the example.
      
      
      Dq.1 <- Etheta1Func(Etheta1, Etheta2, Etheta3, Etheta4, Etheta5, Etheta6, X1, X2, Y)          
      
      Dq.2 <- Etheta2Func(Etheta1, Etheta2, Etheta3, Etheta4, Etheta5, Etheta6, X1, X2, Y)          
      
      Dq.3 <- Etheta3Func(Etheta1, Etheta2, Etheta3, Etheta4, Etheta5, Etheta6, X1, X2, Y)          
      
      Dq.4 <- Etheta4Func(Etheta1, Etheta2, Etheta3, Etheta4, Etheta5, Etheta6, X1, X2, Y)          
      
      Dq.5 <- Etheta5Func(Etheta1, Etheta2, Etheta3, Etheta4, Etheta5, Etheta6, X1, X2, Y)  
      
      Dq.6 <- Etheta6Func(Etheta1, Etheta2, Etheta3, Etheta4, Etheta5, Etheta6, X1, X2, Y)          
      
      
      # for j=1,...,d.
      Dq   <- c(Dq.1, Dq.2, Dq.3, Dq.4, Dq.5, Dq.6)               # Set Dq,j(m,theta_bar) for
      # j=1,...,d, as a vector.
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
logistic(1000,3,4)
  # For determining m1
#logistic(10,1,1)
# For determining m2
#logistic(10,2,2)