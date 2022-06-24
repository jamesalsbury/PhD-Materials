
####
#Looking at stopping for futility in the case we have
####


#Start off with traditional power


#We assume the following parameters for our model
lambda2 <- 0.04
gamma2 <- 0.8

#Draw the control curve
controltime <- seq(0, 100, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)
plot(controltime, controlsurv, type="l", col="blue")

#We assume that the delay is 3 months
T <- 3