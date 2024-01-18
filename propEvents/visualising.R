lambdac <- log(2)/24
T1 <- 3
HR1 <- 5
T2 <- 7
HR2 <- 2.5
recTime <- 0

numPatients <- 340

#Control survival lines
controlTime1 <- seq(0, T1, by=0.01)
controlSurv1 <- exp(-lambdac*controlTime1)

controlTime2 <- seq(T1, 100, by=0.01)
controlSurv2 <- exp(-lambdac*T1 - lambdac*HR1*(controlTime2-T1))


#Treatment survival lines
treatmentTime1 <- seq(0, T1, by=0.01)
treatmentSurv1 <- exp(-lambdac*treatmentTime1)

treatmentTime2 <- seq(T1, T2, by=0.01)
treatmentSurv2 <- exp(-lambdac*T1 - lambdac*HR1*(treatmentTime2-T1))

treatmentTime3 <- seq(T2, 100, by=0.01)
treatmentSurv3 <- exp(-lambdac*T1 - lambdac*HR1*(T2-T1) - lambdac*HR2*(treatmentTime3-T2))

# Plotting
plot(controlTime1, controlSurv1, type = "l", ylim = c(0,1), xlim = c(0,60), xlab = "Time", ylab = "Survival", col = "blue")
lines(controlTime2, controlSurv2, col = "blue")
lines(treatmentTime1, treatmentSurv1, col = "red", lty = 2)
lines(treatmentTime2, treatmentSurv2, col = "red", lty = 2)
lines(treatmentTime3, treatmentSurv3, col = "red", lty  = 2)
legend("topright", legend = c("Control", "Treatment"), lty = 1:2, col = c("blue", "red"))

abline(v = 7)


# 
# #Simulate from these survival curves
# CP1 <- exp(-lambdac*T1)
# CP2 <- exp(-lambdac*T1 - lambdac*HR1*(T2-T1))
# u <- runif(numPatients)
# controlData <- ifelse(u>CP1, -log(u)/lambdac, (1/(lambdac*HR1))*(lambdac*HR1*T1-log(u)-lambdac*T1))
# u <- runif(numPatients)
# treatmentData <- ifelse(u>CP1, -log(u)/lambdac, ifelse(u<CP2, (1/(lambdac*HR2))*(lambdac*HR2*T2-log(u)-lambdac*T1-lambdac*HR1*(T2-T1)), (1/(lambdac*HR1))*(lambdac*HR1*T1-log(u)-lambdac*T1)))
# 
# controlData <- data.frame(time = controlData, status = 1)
# kmcontrol <- survfit(Surv(time, status)~1, data = controlData)
# #lines(kmcontrol, conf.int = F, col = "red", lty = 2)
# 
# treatmentData <- data.frame(time = treatmentData, status = 1)
# kmtreatment <- survfit(Surv(time, status)~1, data = treatmentData)
# 
# lines(kmtreatment, conf.int = F, col = "red", lty = 2)
