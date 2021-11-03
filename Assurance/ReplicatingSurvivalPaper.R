
#Exponential example


###Power calculation
N1vec <- seq(1:500)
N2vec <- seq(1:500)
alpha <- 0.05


#Function to calculate the power for two user-specified sample sizes
powerfunc <- function(N1, N2){
  lambda1 <- -log(0.6)/5
  lambda2 <- -log(0.8)/5
  theta <- log(lambda2/lambda1)
  T <- 5
  R <- 3
  P1e <- 1-(exp(-lambda1*(T-R))-exp(-lambda1*T))/(lambda1*R)
  P2e <- 1-(exp(-lambda2*(T-R))-exp(-lambda2*T))/(lambda2*R)
  power <- vector(length = 500)
  pnorm(-theta/(sqrt(1/(N1*P1e)+1/(N2*P2e)))-qnorm(alpha/2), lower.tail = F)+pnorm(theta/(sqrt(1/(N1*P1e)+1/(N2*P2e)))-qnorm(alpha/2), lower.tail = F)
}

#Calculating power for different values of the sample sizes
power <- vector(length = 500)
for (i in 1:500){
  power[i] <- powerfunc(N1vec[i], N2vec[i])  
}

#Plotting the power calculations
plot(power, type = "l", ylim = c(0,1), xaxt = "n", xlab = "Total sample size", ylab = "Power/Assurance", lty=2)
axis(1, at = seq(0, 500,by = 100), labels = seq(0, 1000, by = 200))

###Assurance methods now

#function to calculate assurance for user-specified sample sizes (m is mean of rho, v is variance of rho)
assurancefunc <- function(N1, N2, m, v){
  M <- 1000
  T <- 5
  R <- 3
  assurance <- vector(length = M)
  for (j in 1:M){
    S1 <- rbeta(1, 60, 40)
    rho <- rnorm(1, m, sqrt(v))
    S2 <- S1+rho
    lambda1 <- -log(S1)/5
    lambda2 <- -log(S2)/5
    theta <- log(lambda2/lambda1)
    P1e <- 1-(exp(-lambda1*(T-R))-exp(-lambda1*T))/(lambda1*R)
    P2e <- 1-(exp(-lambda2*(T-R))-exp(-lambda2*T))/(lambda2*R)
    assurance[j] <- pnorm(theta/(sqrt(1/(N1*P1e)+1/(N2*P2e)))-qnorm(alpha/2), lower.tail = F)
  }
  assurance <- na.omit(assurance)
  mean(assurance)  
}


#Need to calculate assurance for the three different scenarios

#Scenario 1

ass1 <- vector(length = 500)
for (i in 1:500){
  ass1[i] <- assurancefunc(N1vec[i], N2vec[i], m = 0.2, v = 0.001)  
}
lines(ass1, lty = 3, col = "blue")

#Scenario 2

ass2 <- vector(length = 500)
for (i in 1:500){
  ass2[i] <- assurancefunc(N1vec[i], N2vec[i], m = 0.2, v = 0.05)  
}
lines(ass2, lty = 4, col = "purple")


#Scenario 3

ass3 <- vector(length = 500)
for (i in 1:500){
  ass3[i] <- assurancefunc(N1vec[i], N2vec[i], m = 0.3, v = 0.005)  
}
lines(ass3, lty = 5, col = "red")

legend("bottomright", legend = c("Power", "Assurance: scenario 1", "Assurance: scenario 2", "Assurance: scenario 3"),
       col = c("black", "blue", "purple", "red"), lty = 2:5, cex = 0.8)


###Weibull example
N1vec <- seq(1:1000)
N2vec <- seq(1:1000)
alpha <- 0.05

#Function to calculate the power for two user-specified sample sizes
powerfunc <- function(N1, N2){
  kappa1 <- log(log(0.1)/log(0.2))/(log(2/1))
  kappa2 <- log(log(0.2)/log(0.3))/(log(2/1))
  lambda1 <- -log(0.2)/(1^kappa1)
  lambda2 <- -log(0.3)/(1^kappa1)
  mu1 <- gamma(1+1/kappa1)/(lambda1^(1/kappa1))
  mu2 <- gamma(1+1/kappa2)/(lambda2^(1/kappa2))
  sigma1 <- (gamma(1+2/kappa1)-(gamma(1+1/kappa1))^2)/lambda1^(2/kappa1)
  sigma2 <- (gamma(1+2/kappa2)-(gamma(1+1/kappa2))^2)/lambda2^(2/kappa2)
  alpha <- 0.05
  pnorm((mu2-mu1)/(sqrt(sigma1/N1+sigma2/N2))-qnorm(alpha/2), lower.tail = F)+pnorm(-(mu2-mu1)/(sqrt(sigma1/N1+sigma2/N2))-qnorm(alpha/2), lower.tail = F)
}

#Calculating power for different values of the sample sizes
power <- vector(length = 1000)
for (i in 1:1000){
  power[i] <- powerfunc(N1vec[i], N2vec[i])  
}

#Plotting the power calculations
plot(power, type = "l", ylim = c(0,1), xaxt = "n", xlab = "Total sample size", ylab = "Power/Assurance", lty=2)
axis(1, at = seq(0, 1000,by = 250), labels = seq(0, 2000, by = 500))


#Assurance calculations now


assurancefunc <- function(N1, N2){
  M = 500
  assurance <- vector(length = M)
  for (j in 1:M){
    S1t0 <- rbeta(1, 8.10, 32.81)
    Delta11 <- rbeta(1, 4.36, 32.61)
    Delta12 <- rnorm(1, 0.097, sqrt(0.004))
    Delta22 <- rbeta(1, 2.23, 20.10)
    S1t0prime <- S1t0 - Delta11
    S2t0 <- S1t0 + Delta12
    S2t0prime <- S2t0 - Delta22
    kappa1 <- log(log(S1t0prime)/log(S1t0))/(log(2/1))
    kappa2 <- log(log(S2t0prime)/log(S2t0))/(log(2/1))
    lambda1 <- -log(S1t0)
    lambda2 <- -log(S2t0)
    if (is.na(kappa1)|is.na(kappa2)){
      assurance[j] = NA
    } else{
      Survtimes1 <- rweibull(N1, scale = lambda1, shape = kappa1)
      Survtimes2 <- rweibull(N2, scale = lambda2, shape = kappa2)
      X2bar = mean(Survtimes2)
      X1bar = mean(Survtimes1)
      sigma1 = var(Survtimes1)
      sigma2 = var(Survtimes2)
      assurance[j]   =  (X2bar - X1bar)/(sqrt(sigma1/N1+sigma2/N2))>qt(0.975 , df = N1+N2-1)
    }
  }
  assurance <- na.omit(assurance)
  mean(assurance)  
}


N1vec <- seq(1,1000)
N2vec <- seq(1,1000)

ass <- vector(length = 1000)
for (i in 1:1000){
  ass[i] <- assurancefunc(N1vec[i], N2vec[i])  
}

lines(ass, lty=3, col="blue")

##Nonparametric example

Standard = data.frame(Time = c(0.8, 1.0, 2.7, 3.1, 5.4, 7.0, 9.2, 12.1), Cens = c(1, 0, 0, 1, 1, 0, 1, 0))
fit = survfit(Surv(Time, Cens)~1, data = Standard)
plot(fit, conf.int = F)
summary(fit)

