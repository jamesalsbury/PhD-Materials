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

linearRegression <- function(M,d.ini,d.end)
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
  
  # Set-ups for the study design of the example.
  d      <- 3  # d = dim(theta_)
  
  
  
  # Modification is required for the specifications in STEPs 1 and 2.
  # STEP 1
  # 1-1. Specify the prior p(theta_)
  # theta.1 ~ N(a1,b1) 
  # theta.2 ~ N(a2,b2) 
  # theta.3 ~ Ga(a3, b3)
  
  a1 <- 0
  b1 <- 100
  a2 <- 0
  b2 <- 10
  a3 <- 1
  b3 <- 1
  # 1-2. Specify theta_bar, the prior mean under p(theta_)
  # Refer to suitable textbooks such as Gelman et al. (2004)
  Etheta1 <- a1
  Etheta2 <- a2
  Etheta3 <- a3/b3
  # 1-3. Specify the hyperparamters of the epsilon-information prior
  #q0(theta_)
  # Refer to Table 1 of the paper */
  a1.0 <- a1
  b1.0 <- b1*c
  a2.0 <- a2
  b2.0 <- b2*c
  a3.0 <- a3/c
  b3.0 <- b3/c
  # STEP 2
  # 2-1. Compute the infomation matrix of p(theta_)
  Dp      <- c(1/b1, 1/b2, (a3-1)/(Etheta3)^2)        # Specify Dp,j(theta_) for
  #j=1,...,d.
  Dp.plus <- sum(Dp[d.ini:d.end])
  # 2-2. Compute the expected information matrix of qm(theta_|data)
  Dq0      <- c(1/b1.0, 1/b2.0, (a3.0-1)/(Etheta3)^2)     # Specify - 2nd derivative of the
  #log {q0(theta_)}
  # regarding theta.j for j=1,...,d..
  Dq0.plus <- sum(Dq0[d.ini:d.end])
  for (t in 1:T) { # Simulate Monte Carlo samples
    DqYm.out <- numeric(M)
    DqY <- numeric(d)
    
    
    
    Dq.3 <- ((a3/c)-1)/Etheta3^2
    
    
    for (i in 1:M) {
      x <- rnorm(1)
     
      # Specify the model employed
      # in the example.
      Dq.1 <- 1/(c*b1)+Etheta3*i              # Specify Dq,j(m,theta_bar)
      #for j=1,...,d.
      Dq.2 <- 1/(c*b2)+Etheta3*x^2            # Specify Dq,j(m,theta_bar)
      # for j=1,...,d.
      Dq.3 <- i/(2*Etheta3^2)        # Specify Dq,j(m,theta_bar)
      # for j=1,...,d.
      
   # print(Dq.3)
      
      Dq   <- c(Dq.1, Dq.2, Dq.3)               # Set Dq,j(m,theta_bar) for
      # j=1,...,d, as a vector.
      DqY  <- DqY + Dq
      
     # print(DqY)
      Dqm.plus    <- sum(DqY[d.ini:d.end])
      
      #print(Dqm.plus)
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

# For determining m1
linearRegression(10,1,2)

# For determining m1
linearRegression(10,3,3)

