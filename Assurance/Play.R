AssFunc <- function(n1, n2, mu, var){
  
  covs <- data.frame(id = 1:n1)
  ass <- rep(NA, 100)
  expectedevents <- 0.7*n1+0.5*n2
  for (i in 1:100){
    S1 <- rbeta(1, 30, 70)
    rho <- rtruncnorm(1, a = -S1, b = 1-S1, mean = mu, sd = sqrt(var))
    S2 <- S1 + rho
    lambda1 <- -log(S2)/12
    lambda2 <- -log(S1)/12
    Control <- data.frame(simsurv(dist = "exponential", x = covs, lambdas = lambda2, maxt = 12), group = rep("control", n1))
    Treatment <- data.frame(simsurv(dist = "exponential", x = covs, lambdas = lambda1, maxt = 12), group = rep("treatment", n2))
    
    
    
    Combined <- rbind(Control, Treatment)
    
    SortedHalfEvents <- Combined[Combined$status==1,][order(Combined[Combined$status==1,]$eventtime),][1:(expectedevents*0.5),]
    
    if (sum(SortedHalfEvents$group=="treatment")==0){
      
    } else if (sum(SortedHalfEvents$group=="treatment")==expectedevents){
      
    } else {
      test1 <- survdiff(Surv(eventtime, status)~group, data = SortedHalfEvents)
      
      if (test1$chisq > qchisq(1-0.0054, 1)){
        ass[i] <- 1
      } else {
        test2 <- survdiff(Surv(eventtime, status)~group, data = Combined)
        if (test2$chisq >qchisq(1-0.0492, 1)){
          ass[i] <- 1
        } else {
          ass[i] <- 0
        }
      }
    }
    
  }
  mean(ass)
}

N1Vec <- seq(40, 300, by=10)
N2Vec <- N1Vec
assvec2 <- rep(NA, length(N1Vec))

for (i in 1:length(N1Vec)){
  assvec2[i] <- AssFunc(N1Vec[i], N2Vec[i], mu = 0.2, var = 0.05)
}

lo3 <- loess(assvec2~N1Vec)













expectedevents <- 0.7*n1+0.5*n2
i <- 1
n2 <- 40
mu <- 0.2
var <- 0.05
n1 <- 40
covs <- data.frame(id = 1:n1)
ass <- rep(NA, 100)
S1 <- rbeta(1, 30, 70)
rho <- rtruncnorm(1, a = -S1, b = 1-S1, mean = mu, sd = sqrt(var))
S2 <- S1 + rho
lambda1 <- -log(S2)/12
lambda2 <- -log(S1)/12
Control <- data.frame(simsurv(dist = "exponential", x = covs, lambdas = lambda2, maxt = 12), group = rep("control", n1))
Treatment <- data.frame(simsurv(dist = "exponential", x = covs, lambdas = lambda1, maxt = 12), group = rep("treatment", n2))

Combined <- rbind(Control, Treatment)

SortedHalfEvents <- Combined[Combined$status==1,][order(Combined[Combined$status==1,]$eventtime),][1:(expectedevents*0.5),]

if (sum(SortedHalfEvents$group=="treatment")==0){
  
} else if (sum(SortedHalfEvents$group=="treatment")==expectedevents){
  
} else {
  test1 <- survdiff(Surv(eventtime, status)~group, data = SortedHalfEvents)
  
  if (test1$chisq > qchisq(1-0.0054, 1)){
    ass[i] <- 1
  } else {
    test2 <- survdiff(Surv(eventtime, status)~group, data = Combined)
    if (test2$chisq >qchisq(1-0.0492, 1)){
      ass[i] <- 1
    } else {
      ass[i] <- 0
    }
  }
}
