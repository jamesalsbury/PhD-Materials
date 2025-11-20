BPPFunc <- function(IAData, censTime){
  
  n <- sum(IAData$group=="Control")
  m <- nrow(IAData)
  
  IAData <- IAData[order(IAData$group),]
  
  IAkmfit <- survfit(Surv(survival_time, status)~group, data = IAData) 
  
  #plot(IAkmfit, col = c("blue", "red"), xlim=c(0,50))
  
  
  
  #We need to update the priors using this interim data
  modelstring="

data {
  for (j in 1:m){
    zeros[j] <- 0
  }
}

model {
  C <- 10000
  for (i in 1:n){
    zeros[i] ~ dpois(zeros.mean[i])
    zeros.mean[i] <-  -l[i] + C
    l[i] <- ifelse(datEvent[i]==1, log(lambda2)-(lambda2*datTimes[i]), -(lambda2*datTimes[i]))
  }
  for (i in (n+1):m){                                                                                                             
    zeros[i] ~ dpois(zeros.mean[i])
    zeros.mean[i] <-  -l[i] + C
    l[i] <- ifelse(datEvent[i]==1, ifelse(datTimes[i]<bigT, log(lambda2)-(lambda2*datTimes[i]), log(lambda1)-lambda1*(datTimes[i]-bigT)-(bigT*lambda2)), 
      ifelse(datTimes[i]<bigT, -(lambda2*datTimes[i]), -(lambda2*bigT)-lambda1*(datTimes[i]-bigT)))
  }
  
    lambda2 ~ dbeta(1,1)
    HR ~ dnorm(0.75, 1/HRsd)T(0,)
    bigT ~ dnorm(3, 1/bigTsd^2)T(0,)
    lambda1 <- lambda2*HR
    
    }
"

model = jags.model(textConnection(modelstring), data = list(datTimes = IAData$survival_time, datEvent = IAData$status, bigTsd = 1, HRsd = 0.1, n= n, m=m), quiet = T) 


update(model, n.iter=100)
output=coda.samples(model=model, variable.names=c("HR", "bigT", "lambda2"), n.iter = 500)

#plot(output)

#The number of unenrolled patients in each group
cPatientsLeft <- (nPatients/2) - sum(IAData$group=="Control") 
tPatientsLeft <- (nPatients/2) - sum(IAData$group=="Treatment") 

BPPVec <- rep(NA, 500)

for (i in 1:500){
  
  #Sampling the recruitment times for the unenrolled patients
  unenrolledRecTimes <- runif(cPatientsLeft+tPatientsLeft, censTime, recTime)
  
  #Extract realisations from the MCMC
  HRoutput <- as.numeric(unlist(output[,1]))
  bigToutput <- as.numeric(unlist(output[,2]))
  lambda2output <- as.numeric(unlist(output[,3]))
  
  #Sample values from the MCMC output
  sampledHR <- sample(HRoutput, 1)
  sampledbigT <- sample(bigToutput, 1)
  sampledlambda2 <- sample(lambda2output, 1)
  sampledlambda1 <- sampledlambda2*sampledHR
  
  
  #For the unenrolled data, we can sample the remaining data according to the updated (sampled) parameters
  CP <- exp(-(sampledlambda2*sampledbigT))
  u <- runif(tPatientsLeft)
  
  unenrolledData <- data.frame(time = c(rexp(cPatientsLeft, rate = sampledlambda2), ifelse(u>CP, (-log(u))/sampledlambda2, (1/sampledlambda1)*(sampledbigT*sampledlambda1-log(u)-sampledbigT*sampledlambda2))), group = c(rep("Control", cPatientsLeft),
                                                                                                                                                                                                                          rep("Treatment", tPatientsLeft)), recTime = unenrolledRecTimes)
  
  unenrolledData$psuedoTime <- unenrolledData$time + unenrolledData$recTime
  
  
  #Extracting the observations that were censored at the IA  
  censoredData <- IAData[IAData$status==0,]
  
  #Number of censored observations in each group
  cCensored <- sum(censoredData$group=="Control")
  tCensored <- sum(censoredData$group=="Treatment")
  
  #Extracting the censored observations in the control group
  cCensoredData <- censoredData %>%
    filter(group=="Control")
  
  #Adding a exp(sampledlambda2) value to the censored value
  cCensoredData$finalsurvTime <- cCensoredData$survival_time + rexp(cCensored, rate = sampledlambda2)
  
  #Calculating the psuedo time
  cCensoredData$finalPsuedoTime <- cCensoredData$recTime + cCensoredData$finalsurvTime
  
  #Extacting the observations in the treatment group which may still be influenced by the delay (their observation time is smaller than the sampled delay time)
  tBeforeDelay <- censoredData %>%
    filter(group=="Treatment") %>%
    filter(survival_time < sampledbigT)
  
  #Extracting the observations in the treatment group which will not be influenced by the delay (their observation time is bigger than the sampled delay time)
  tAfterDelay <- censoredData %>%
    filter(group=="Treatment") %>%
    filter(survival_time > sampledbigT)
  
  #As these observations are still subject to a delay, we add on a Exp(lambda2) (lambdac) time
  tBeforeDelay$IASurv <- tBeforeDelay$survival_time + rexp(nrow(tBeforeDelay), rate = sampledlambda2)
  
  #Extracting the observations in which the survival time is smaller than the sampled delay time
  tBeforeDelay1 <- tBeforeDelay %>%
    filter(IASurv < sampledbigT)
  
  #Extracting the observations in which the survival time is bigger than the sampled delay time
  tBeforeDelay2 <- tBeforeDelay %>%
    filter(IASurv > sampledbigT)
  
  #For the observations in which the survival time is bigger, we sample a Exp(lambda1) and add it to the sampled delay time
  tBeforeDelay2$IASurv2 <- sampledbigT + rexp(nrow(tBeforeDelay2), rate = sampledlambda1)
  
  #For the observations not influenced by the delay, we sample a Exp(lambda1) time and add it to the current survival time
  tAfterDelay$IASurv <- tAfterDelay$survival_time + rexp(nrow(tAfterDelay), rate = sampledlambda1)
  
  #Calculate the pseudo time for all the data frames
  tBeforeDelay1$IApsuedoTime <- tBeforeDelay1$IASurv + tBeforeDelay1$recTime
  tBeforeDelay2$IApsuedoTime <- tBeforeDelay2$IASurv2 + tBeforeDelay2$recTime
  tAfterDelay$IApsuedoTime <- tAfterDelay$IASurv + tAfterDelay$recTime
  
  #Only keeping the columns of interest
  cCensoredData <- cCensoredData[,c(2:3, 8:9)]
  tBeforeDelay1 <- tBeforeDelay1[,c(2:3, 8:9)]
  tBeforeDelay2 <- tBeforeDelay2[,c(2:3, 9:10)]
  tAfterDelay <- tAfterDelay[,c(2:3, 8:9)]
  
  #Keeping the column names consistent
  colnames(cCensoredData) <- c("group", "recTime", "time", "psuedoTime")
  colnames(tBeforeDelay1) <- c("group", "recTime", "time", "psuedoTime")
  colnames(tBeforeDelay2) <- c("group", "recTime", "time", "psuedoTime")
  colnames(tAfterDelay) <- c("group", "recTime", "time", "psuedoTime")
  
  #Only keeping observations which are dead
  finalDataset <- IAData %>%
    filter(status==1)
  
  finalDataset <- finalDataset[,1:4]
  
  #Combining all the above data sets 
  finalDataset <- rbind(finalDataset, tBeforeDelay1)
  finalDataset <- rbind(finalDataset, tBeforeDelay2)
  finalDataset <- rbind(finalDataset, tAfterDelay)
  finalDataset <- rbind(finalDataset, unenrolledData)
  finalDataset <- rbind(finalDataset, cCensoredData)
  
  #Making sure the final data set is correct
  censTime1 <- sort(finalDataset$psuedoTime)[nEvents]
  finalDataset$status <- finalDataset$psuedoTime <= censTime1
  finalDataset$status <- finalDataset$status*1
  finalDataset$enrolled <- finalDataset$recTime <= censTime1
  finalDataset <-  finalDataset[finalDataset$enrolled==T,]
  finalDataset$survival_time <- ifelse(finalDataset$psuedoTime>censTime1, censTime1  - finalDataset$recTime, finalDataset$time)
  
  #Testing for significance
  test <- survdiff(Surv(survival_time, status)~group, data = finalDataset)
  
  #kmfit <- survfit(Surv(survival_time, status)~group, data = finalDataset)
  #plot(kmfit, col = c("blue", "red"), xlim=c(0,50))
  
  #Making sure the significance is in the correct direction
  coxmodel <- coxph(Surv(survival_time, status)~group, data = finalDataset)
  deltad <- as.numeric(exp(coef(coxmodel)))
  
  
  BPPVec[i] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
  
}


return(BPP = mean(BPPVec))
}

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


no_IA_func <- function(df){
  
  df <- df[order(df$psuedo_time), ]
  final_time <- df$psuedo_time[num_events]
  
  df_final <- df[df$rec_time <= final_time, ]
  df_final$event <- df_final$psuedo_time <= final_time
  df_final$obs_time <- pmin(df_final$time, final_time - df_final$rec_time)
  
  fit <- coxph(Surv(obs_time, event) ~ group, data = df_final)
  test_statistic <- summary(fit)$sctest["test"]
  observed_HR <- exp(coef(fit))
  
  power <- (observed_HR < 1) & (test_statistic > qchisq(0.95, 1))
  duration <- final_time
  ss <- nrow(df_final)
  
  return(list(power = power,
              duration = duration,
              ss = ss))
  
}


wieand_func <- function(df){
  
  df <- df[order(df$psuedo_time), ]
  interim1_time <- df$psuedo_time[num_events*0.5]
  interim2_time <- df$psuedo_time[num_events*0.75]
  final_time <- df$psuedo_time[num_events]
  
  # Evaluate at first interim
  df1 <- df[df$rec_time <= interim1_time, ]
  df1$event <- df1$psuedo_time <= interim1_time
  df1$obs_time <- pmin(df1$time, interim1_time - df1$rec_time)
  
  fit1 <- coxph(Surv(obs_time, event) ~ group, data = df1)
  hr1 <- exp(coef(fit1))
  if (hr1 > 1){
    return(list(power = 0, duration = interim1_time, ss = nrow(df1),
                interim = 1))
    
  }
  
  # Evaluate at second interim if not stopped
  
  df2 <- df[df$rec_time <= interim2_time, ]
  df2$event <- df2$psuedo_time <= interim2_time
  df2$obs_time <- pmin(df2$time, interim2_time - df2$rec_time)
  
  fit2 <- coxph(Surv(obs_time, event) ~ group, data = df2)
  hr2 <- exp(coef(fit2))
  if (hr2 > 1){
    return(list(power = 0,  duration = interim2_time, ss = nrow(df2),
                interim = 2))
  }
  
  
  
  # Evaluate at final analysis if not stopped
  
  
  df_final <- df[df$rec_time <= final_time, ]
  df_final$event <- df_final$psuedo_time <= final_time
  df_final$obs_time <- pmin(df_final$time, final_time - df_final$rec_time)
  
  fit <- coxph(Surv(obs_time, event) ~ group, data = df_final)
  test_statistic <- summary(fit)$sctest["test"]
  observed_HR <- exp(coef(fit))
  
  power <- (observed_HR < 1) & (test_statistic > qchisq(0.95, 1))
  duration <- final_time
  ss <- nrow(df_final)
  
  
  
  return(list(power = power,
              duration = duration,
              ss = ss, interim = 3))
  
}


OBF_func <- function(df){
  
  df <- df[order(df$psuedo_time), ]
  interim1_time <- df$psuedo_time[num_events*(1/3)]
  interim2_time <- df$psuedo_time[num_events*(2/3)]
  final_time <- df$psuedo_time[num_events]
  
  # Evaluate at first interim
  df1 <- df[df$rec_time <= interim1_time, ]
  df1$event <- df1$psuedo_time <= interim1_time
  df1$obs_time <- pmin(df1$time, interim1_time - df1$rec_time)
  
  fit1 <- coxph(Surv(obs_time, event) ~ group, data = df1)
  hr1 <- exp(coef(fit1))
  
  
  
  if (hr1 > 0.998){
    return(list(power = 0, duration = interim1_time, ss = nrow(df1),
                interim = 1))
    
  }
  
  # Evaluate at second interim if not stopped
  
  df2 <- df[df$rec_time <= interim2_time, ]
  df2$event <- df2$psuedo_time <= interim2_time
  df2$obs_time <- pmin(df2$time, interim2_time - df2$rec_time)
  
  fit2 <- coxph(Surv(obs_time, event) ~ group, data = df2)
  hr2 <- exp(coef(fit2))
  if (hr2 > 0.913){
    return(list(power = 0, duration = interim2_time, ss = nrow(df2),
                interim = 2))
  }
  
  
  
  # Evaluate at final analysis if not stopped
  
  
  df_final <- df[df$rec_time <= final_time, ]
  df_final$event <- df_final$psuedo_time <= final_time
  df_final$obs_time <- pmin(df_final$time, final_time - df_final$rec_time)
  
  fit <- coxph(Surv(obs_time, event) ~ group, data = df_final)
  test_statistic <- summary(fit)$sctest["test"]
  observed_HR <- exp(coef(fit))
  
  power <- (observed_HR < 1) & (test_statistic > qchisq(0.95, 1))
  duration <- final_time
  ss <- nrow(df_final)
  
  
  
  return(list(power = power,
              duration = duration,
              ss = ss, interim = 3))
}


proposed_func <- function(df) {
  df <- df[order(df$psuedo_time), ]
  
  # Event thresholds
  target1 <- as.integer(ceiling(num_events * 0.50))
  target2 <- as.integer(ceiling(num_events * 0.75))
  final_time <- df$psuedo_time[num_events]
  
  # Function to find earliest time meeting both event and delay conditions
  find_earliest_time <- function(df, target_events) {
    for (i in seq_len(nrow(df))) {
      t_i <- df$psuedo_time[i]
      df_i <- df[df$rec_time <= t_i, ]
      df_i$event <- df_i$psuedo_time <= t_i
      
      n_events <- sum(df_i$event)
      if (n_events >= target_events && n_events > 0) {
        prop_delayed <- mean((df_i$psuedo_time - df_i$rec_time)[df_i$event] >= 3)
        if (!is.na(prop_delayed) && prop_delayed >= 2/3) {
          return(list(time = t_i, prop_delayed = prop_delayed))
        }
      }
    }
    return(NULL)
  }
  
  # ---- Interim 1 ----
  ia1 <- find_earliest_time(df, target1)
  if (!is.null(ia1) && ia1$time <= final_time) {
    interim1_time <- ia1$time
    df1 <- df[df$rec_time <= interim1_time, ]
    df1$event <- df1$psuedo_time <= interim1_time
    df1$obs_time <- pmin(df1$time, interim1_time - df1$rec_time)
    
    fit1 <- coxph(Surv(obs_time, event) ~ group, data = df1)
    hr1 <- exp(coef(fit1))
    
    if (hr1 > 1) {
      return(list(power = 0, duration = interim1_time, ss = nrow(df1),
                  interim = 1, hr = hr1, prop_delayed = ia1$prop_delayed))
    }
  }
  
  # ---- Interim 2 ----
  ia2 <- find_earliest_time(df, target2)
  if (!is.null(ia2) && ia2$time <= final_time) {
    interim2_time <- ia2$time
    df2 <- df[df$rec_time <= interim2_time, ]
    df2$event <- df2$psuedo_time <= interim2_time
    df2$obs_time <- pmin(df2$time, interim2_time - df2$rec_time)
    
    fit2 <- coxph(Surv(obs_time, event) ~ group, data = df2)
    hr2 <- exp(coef(fit2))
    
    if (hr2 > 1) {
      return(list(power = 0, duration = interim2_time, ss = nrow(df2),
                  interim = 2))
    }
  }
  
  # ---- Final analysis ----
  df_final <- df[df$rec_time <= final_time, ]
  df_final$event <- df_final$psuedo_time <= final_time
  df_final$obs_time <- pmin(df_final$time, final_time - df_final$rec_time)
  
  fit <- coxph(Surv(obs_time, event) ~ group, data = df_final)
  test_statistic <- summary(fit)$sctest["test"]
  observed_HR <- exp(coef(fit))
  
  power <- (observed_HR < 1) & (test_statistic > qchisq(0.95, 1))
  duration <- final_time
  ss <- nrow(df_final)
  
  return(list(power = power, duration = duration, ss = ss, interim = 3))
}

reproduce_func <- function(scenario, rec_duration, n_sims) {
  
  wieand_list <- vector("list", length = n_sims)
  OBF_list <- vector("list", length = n_sims)
  no_IA_list <- vector("list", length = n_sims)
  proposed_list <- vector("list", length = n_sims)
  
  for (i in 1:n_sims) {
    u <- runif(n_c)
    control_times <- -log(u) / scen_list[[scenario]]$lambda_c
    
    u <- runif(n_t)
    CP <- exp(-scen_list[[scenario]]$lambda_c * scen_list[[scenario]]$HR1* scen_list[[scenario]]$delay)
    treatment_times <- ifelse(
      u > CP,
      -log(u) / (scen_list[[scenario]]$lambda_c * scen_list[[scenario]]$HR1),
      (1 / (scen_list[[scenario]]$lambda_c * scen_list[[scenario]]$HR2)) * (scen_list[[scenario]]$lambda_c*scen_list[[scenario]]$HR2*scen_list[[scenario]]$delay -log(u) - scen_list[[scenario]]$lambda_c * scen_list[[scenario]]$HR1 * scen_list[[scenario]]$delay)
    )
    
    df <- data.frame(
      time = c(control_times, treatment_times),
      group = c(rep("Control", n_c), rep("Treatment", n_t))
    )
    
    df$rec_time <- runif(n_c + n_t, 0, rec_duration)
    df$psuedo_time <- df$time + df$rec_time
    
    wieand_list[[i]] <- wieand_func(df)
    OBF_list[[i]] <- OBF_func(df)
    no_IA_list[[i]] <- no_IA_func(df)
    proposed_list[[i]] <- proposed_func(df)
    
    
    
  }
  
  return(list(no_IA_list = no_IA_list,
              wieand_list = wieand_list,
              OBF_list = OBF_list,
              proposed_list = proposed_list))
  
}


plot_events_v_prop_func <- function(scenario, rec_duration){
  
  # Example run
  x1 <- propEventFunc(lambda_c,
                      scen_list[[scenario]]$HR1, 
                      scen_list[[scenario]]$delay, 
                      scen_list[[scenario]]$HR2,
                      n_c, 
                      rec_duration)
  
  plot(x1$eventVec, x1$eventProp, type = "l", ylim = c(0, 1),
       xlab = "Number of events", ylab = "Proportion of events > 3 months",
       col = "red", 
       main = paste0("Scenario ", scenario, ", Recruitment: ", rec_duration), 
       cex.axis=1.5, cex.lab=1.5, cex.main=2)
  abline(h = 2/3, lty = 2)
  abline(v = 512 * 0.5, lty = 3)
  abline(v = 512 * 0.75, lty = 3)
  abline(v = 512, lty = 3)
}

plot_break_scenarios <- function(scenario, rec_duration){
  
  x1 <- propEventFunc(scen_list[[scenario]]$lambda_c,
                      scen_list[[scenario]]$HR1, 
                      scen_list[[scenario]]$delay, 
                      scen_list[[scenario]]$HR2,
                      n_c, 
                      rec_duration)
  
  
  par(mfrow = c(1,2))
  plot(x1$eventVec, x1$eventProp, type = "l", ylim = c(0, 1),
       xlab = "Number of events", ylab = "Proportion of events > 3 months",
       col = "red", 
       #main = paste0("Scenario ", scenario, ", Recruitment: ", rec_duration), 
       cex.axis=1.5, cex.lab=1.5, cex.main=2)
  abline(h = 2/3, lty = 2)
  abline(v = 512 * 0.5, lty = 3)
  abline(v = 512 * 0.75, lty = 3)
  abline(v = 512, lty = 3)
  
  
  IA_Time1 <- min(max(x1$calTime[which(x1$eventProp > (2/3))[1]], x1$calTime[256]), x1$calTime[512])
  IA_Time2 <- min(max(x1$calTime[which(x1$eventProp > (2/3))[1]], x1$calTime[384]), x1$calTime[512])
  Final_Time <- x1$calTime[512]
  
  #Control survival lines
  controlTime1 <- seq(0, 100, by=0.01)
  controlSurv1 <- exp(-scen_list[[scenario]]$lambda_c*controlTime1)
  
  #Treatment survival lines
  treatmentTime1 <- seq(0, scen_list[[scenario]]$delay, by=0.01)
  treatmentSurv1 <- exp(-scen_list[[scenario]]$lambda_c*treatmentTime1*scen_list[[scenario]]$HR1)
  
  treatmentTime2 <- seq(scen_list[[scenario]]$delay, 100, by=0.01)
  treatmentSurv2 <- exp(-scen_list[[scenario]]$lambda_c*scen_list[[scenario]]$delay*scen_list[[scenario]]$HR1 - scen_list[[scenario]]$lambda_c*scen_list[[scenario]]$HR2*(treatmentTime2-scen_list[[scenario]]$delay))
  
  
  # Plotting
  plot(controlTime1, controlSurv1, type = "l", ylim = c(0,1),
       xlim = c(0,60), xlab = "Time",
       ylab = "Survival", col = "blue",
       cex.axis=1.5, cex.lab=1.5, cex.main=2)
  lines(treatmentTime1, treatmentSurv1, col = "red", lty = 1)
  lines(treatmentTime2, treatmentSurv2, col = "red", lty = 1)
  legend("topright", legend = c("Control", "Treatment", "Interim Analysis 1", "Interim Analysis 2", "Final Analysis"),
         lty = c(1, 1, 2, 3, 4), 
         col = c("blue", "red", "black", "black", "black"))
  abline(v = IA_Time1, lty = 2)
  abline(v = IA_Time2, lty = 3)
  abline(v = Final_Time, lty = 4)
  
}

