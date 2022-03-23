

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

#Now we calculate assurance for these chosen inputs - the way we have previously been doing

bigT <- rnorm(1, 3, sd = sqrt(0.01))
HR <- rbeta(1, 6, 4)
lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
gamma1 <- gamma2

treatmenttime <- seq(bigT, exp((1.527/gamma2)-log(lambda2))*1.1, by=0.01)
treatmentsurv <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))

timecombined <- c(controltime, treatmenttime)
survcombined <- c(controlsurv, treatmentsurv)

lines(treatmenttime, treatmentsurv, col="red")











