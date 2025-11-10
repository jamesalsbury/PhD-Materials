# Efficient data generator
generateData <- function(lambdac, HR1, T1, HR2, numPatients, recTime) {
  # control
  u_c <- runif(numPatients)
  control_time <- -log(u_c)/lambdac
  
  # treatment
  u_t <- runif(numPatients)
  CP <- exp(-lambdac * HR1 * T1)
  treatment_time <- numeric(numPatients)
  idx <- u_t > CP
  treatment_time[idx] <- -log(u_t[idx]) / (lambdac * HR1)
  treatment_time[!idx] <- (1/(lambdac * HR2)) * (lambdac * HR2 * T1 - log(u_t[!idx]) - lambdac * HR1 * T1)
  
  time <- c(control_time, treatment_time)
  group <- rep(c("Control", "Treatment"), each = numPatients)
  rec_time <- runif(numPatients * 2, 0, recTime)
  pseudo_time <- time + rec_time
  
  data.frame(time, group, rec_time, pseudo_time)
}

# Vectorized censoring calculation
censFunc <- function(time, rec_time, pseudo_time, censTime) {
  enrolled <- rec_time <= censTime
  status <- as.integer(pseudo_time <= censTime)
  surv_time <- ifelse(status == 1, time, censTime - rec_time)
  list(status = status[enrolled], surv_time = surv_time[enrolled])
}

# Vectorized proportion of events > 3 months
propEventsGivenCensor <- function(calTimes, time, rec_time, pseudo_time, threshold = 3) {
  sapply(calTimes, function(ct) {
    tmp <- censFunc(time, rec_time, pseudo_time, ct)
    mean(tmp$surv_time[tmp$status == 1] > threshold)
  })
}

# Main function
propEventFunc <- function(lambdac, HR1, T1, HR2, numPatients, recTime, numSimulations = 100) {
  # Preallocate list for results
  eventProps <- vector("list", numSimulations)
  
  for (k in seq_len(numSimulations)) {
    dat <- generateData(lambdac, HR1, T1, HR2, numPatients, recTime)
    dat <- dat[order(dat$pseudo_time), ]
    eventProps[[k]] <- propEventsGivenCensor(
      calTimes = dat$pseudo_time,
      time = dat$time,
      rec_time = dat$rec_time,
      pseudo_time = dat$pseudo_time
    )
  }
  
  # Combine results into a matrix
  propMatrix <- do.call(rbind, eventProps)
  
  # Compute summary statistics
  eventMean <- colMeans(propMatrix)
  eventSD   <- apply(propMatrix, 2, sd)
  eventSE   <- eventSD / sqrt(numSimulations)
  
  list(
    eventProp = eventMean,
    eventVec  = seq_len(numPatients * 2),
    calTime   = dat$pseudo_time,
    eventMean = eventMean,
    eventSD   = eventSD,
    eventSE   = eventSE
  )
}




scen_list <- list(A = list(HR1 = 0.75,
                           HR2 = 0.75,
                           delay = 0),
                  B = list(HR1 = 1,
                           HR2 = 1,
                           delay = 0),
                  C = list(HR1 = 1.3,
                           HR2 = 1.3,
                           delay = 0),
                  D = list(HR1 = 1,
                           HR2 = 0.693,
                           delay = 3),
                  E = list(HR1 = 1,
                           HR2 = 0.620,
                           delay = 6),
                  F = list(HR1 = 1.3,
                           HR2 = 0.628,
                           delay = 3)
)

lambda_c <- log(2)/12
n_c <- 340
n_t <- 340
num_events <- 512





events_v_prop_func <- function(scenario, rec_duration){
  
  # Example run
  x1 <- propEventFunc(lambda_c,
                      scen_list[[scenario]]$HR1, 
                      scen_list[[scenario]]$delay, 
                      scen_list[[scenario]]$HR2,
                      n_c, 
                      rec_duration)
  
  #plotTime <- seq(1, max(x1$calTime), by = 6)
  #plotEvents <- sapply(plotTime, function(pt) sum(x1$calTime < pt))
  
  plot(x1$eventVec, x1$eventProp, type = "l", ylim = c(0, 1),
       xlab = "Number of events", ylab = "Proportion of events > 3 months",
       col = "red", 
       main = paste0("Scenario ", scenario, ", Recruitment: ", rec_duration), 
       cex.axis=1.5, cex.lab=1.5, cex.main=2)
  lines(x1$eventVec, eventMean + 1.96 * eventSE, lty = 2, col = "gray")
  lines(x1$eventVec, eventMean - 1.96 * eventSE, lty = 2, col = "gray")
  #points(x1$eventVec[plotEvents], x1$eventProp[plotEvents])
  abline(h = 2/3, lty = 2)
  abline(v = 512 * 0.5, lty = 3)
  abline(v = 512 * 0.75, lty = 3)
  abline(v = 512, lty = 3)
  
  #print(which(x1$eventProp > (2/3))[1])
}
# 
# for (i in 1:6){
#   events_v_prop_func(i, 34)
# }
# 
# for (i in 1:6){
#   events_v_prop_func(i, 12)
# }

png("Events_Prop_Scen1.png", units="in", width=12, height=5, res=700)
par(mfrow=c(1,2))
events_v_prop_func(1, 34)
events_v_prop_func(1, 12)
dev.off()

png("Events_Prop_Scen6.png", units="in", width=12, height=5, res=700)
par(mfrow=c(1,2))
events_v_prop_func(6, 34)
events_v_prop_func(6, 12)
dev.off()






#############
#Finding parameters that 'break' it
##############


lambda_c <- log(2)/18
HR1 <- 1.3
HR2 <- 0.6
delay <- 5 
n_c <- 340
rec_duration <- 0

x1 <- propEventFunc(lambda_c,
                    HR1, 
                    delay, 
                    HR2,
                    n_c, 
                    rec_duration)

plotTime <- seq(1, max(x1$calTime), by = 6)
plotEvents <- sapply(plotTime, function(pt) sum(x1$calTime < pt))

plot(x1$eventVec, x1$eventProp, type = "l", ylim = c(0, 1),
     xlab = "Number of events", ylab = "Proportion of events > 3 months",
     col = "red", 
     #main = paste0("Scenario ", scenario, ", Recruitment: ", rec_duration), 
     cex.axis=1.5, cex.lab=1.5, cex.main=2)
lines(x1$eventVec, x1$eventMean + 1.96 * x1$eventSE, lty = 2, col = "gray")
lines(x1$eventVec, x1$eventMean - 1.96 * x1$eventSE, lty = 2, col = "gray")
points(x1$eventVec[plotEvents], x1$eventProp[plotEvents])
abline(h = 2/3, lty = 2)
abline(v = 512 * 0.5, lty = 3)
abline(v = 512 * 0.75, lty = 3)
abline(v = 512, lty = 3)


print(which(x1$eventProp > (2/3))[1])
IA_Time <- max(x1$calTime[which(x1$eventProp > (2/3))[1]], x1$calTime[256])
print(IA_Time)
Final_Time <- x1$calTime[512]
print(Final_Time)

#Control survival lines
controlTime1 <- seq(0, 100, by=0.01)
controlSurv1 <- exp(-lambda_c*controlTime1)

#Treatment survival lines
treatmentTime1 <- seq(0, delay, by=0.01)
treatmentSurv1 <- exp(-lambda_c*treatmentTime1*HR1)

treatmentTime2 <- seq(delay, 100, by=0.01)
treatmentSurv2 <- exp(-lambda_c*delay*HR1 - lambda_c*HR2*(treatmentTime2-delay))


# Plotting
plot(controlTime1, controlSurv1, type = "l", ylim = c(0,1), xlim = c(0,60), xlab = "Time", ylab = "Survival", col = "blue")
lines(treatmentTime1, treatmentSurv1, col = "red", lty = 2)
lines(treatmentTime2, treatmentSurv2, col = "red", lty = 2)
legend("topright", legend = c("Control", "Treatment"), lty = 1:2, col = c("blue", "red"))
abline(v = IA_Time)
abline(v = Final_Time)










