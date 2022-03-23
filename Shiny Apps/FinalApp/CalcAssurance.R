

#Need to look at what we actually need to compute assurance??
#We are sampling from the control and treatment curves
#But maybe we only need number of events that occur (or do not occur) by the time-point that the assurance is going to take place>

#We could compute assurance in the way that we have done previously
#Then find a 'quicker' way by just estimating the number of events


#We assume the control is lambda2 = 0.06, gamma2 = 0.8

lambda2 <- 0.06
gamma2 <- 0.8

controltime <- seq(0, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)
plot(controltime, controlsurv, type="l", col="blue", ylab="Survival", xlab = "Time")

#We assume T comes from a Normal(3, 0.01) dist

x <- seq(0, 6, by=0.01)
y <- dnorm(x, mean = 3, sd = sqrt(0.01))
plot(x, y, type="l")

#We assume HR comes from a Beta(6, 4) dist

x <- seq(0, 1, by=0.01)
y <- dbeta(x, 6, 4)
plot(x, y, type="l")

#Plot a sampled treatment curve

bigT <- rnorm(1, 3, sd = sqrt(0.01))
HR <- rbeta(1, 6, 4)
lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
gamma1 <- gamma2

treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))

timecombined <- c(controltime, treatmenttime)
survcombined <- c(controlsurv, treatmentsurv)

lines(treatmenttime, treatmentsurv, col="red")

#Calculate assurance the way we have been doing
#Firstly we assume n1 = 100, n2 = 100
#We assume the trial will be stopped after 40 months
#0 = alive
#1 = dead

n1 <- 100
n2 <- 100
assvec <- rep(NA, 10000)
gamma1 <- gamma2

for (i in 1:10000){
  bigT <- rnorm(1, 3, sd = sqrt(0.01))
  HR <- rbeta(1, 6, 4)
  lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
  CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
  u <- runif(n2)
  
  suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
  
  dataCombined <- data.frame(Time = c(rweibull(n1, gamma2, 1/lambda2), z), Group = c(rep("Control", n1), rep("Treatment", n2)), 
                             Cens = rep(1, n1+n2))
  
  dataCombined$Cens <-  dataCombined$Time < 40
  dataCombined$Cens <- dataCombined$Cens*1
  test <- survdiff(Surv(Time, Cens)~Group, data = dataCombined)
  assvec[i] <-test$chisq > qchisq(0.95, 1)
}

mean(assvec)

# 'True' assurance is 0.591

#Can we replicate this by just looking at the survival probabilities?

#We can calculate how many patients we expect to be alive after 40 months in the control group fairly easily

n1 <- 100
n2 <- 100
lambda2 <- 0.06
gamma2 <- 0.8
gamma1 <- gamma2
bigT <- rnorm(1, 3, sd = sqrt(0.01))
HR <- rbeta(1, 6, 4)
lambda1 <- exp((log(HR)/gamma2)+log(lambda2))

controlprob <- exp(-(lambda2*40)^gamma2)*n1
treatmentprob <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(40^gamma1-bigT^gamma1))


treatmentprob <- rep(NA, 10000)

for (i in 1:10000){
  bigT <- rnorm(1, 3, sd = sqrt(0.01))
  HR <- rbeta(1, 6, 4)
  lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
  treatmentprob[i] <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(40^gamma1-bigT^gamma1))
}
hist(treatmentprob, breaks = 100, xlim=c(0,0.8))








