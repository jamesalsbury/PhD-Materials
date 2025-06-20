library(survival)
library(rjags)

IPD_INTEREST_OS <- read.csv(file = "Thesis/Chapter 3/Data/IPD_INTEREST_OS.csv")
IPD_ZODIAC_OS <- read.csv(file = "Thesis/Chapter 3/Data/IPD_ZODIAC_OS.csv")
IPD_REVEL_OS <- read.csv(file = "Thesis/Chapter 3/Data/IPD_REVEL_OS.csv")

#png("Doce_Surv.png", units="in", width=10, height=8, res=700)
INTEREST_kmfit <- survfit(Surv(Survival.time, Status)~1, data = IPD_INTEREST_OS)
plot(INTEREST_kmfit, xlim = c(0, 21), col = "blue", conf.int = F, ylab = "Overall Survival", xlab = "Time (months)")
INTEREST_myfit <- survreg(Surv(Survival.time, Status)~1, data = IPD_INTEREST_OS, dist = "exponential")
INTEREST_lambda <- exp(-INTEREST_myfit$coefficients)
time_seq <- seq(0, 21, by=0.01)
surv_probs <- exp(-INTEREST_lambda*time_seq)
#lines(time_seq, surv_probs, col = "blue")


ZODIAC_kmfit <- survfit(Surv(Survival.time, Status)~1, data = IPD_ZODIAC_OS)
lines(ZODIAC_kmfit, xlim = c(0, 21), col = "red", conf.int = F)
ZODIAC_myfit <- survreg(Surv(Survival.time, Status)~1, data = IPD_ZODIAC_OS, dist = "exponential")
ZODIAC_lambda <- exp(-ZODIAC_myfit$coefficients)
surv_probs <- exp(-ZODIAC_lambda*time_seq)
#lines(time_seq, surv_probs, col = "red")

REVEL_kmfit <- survfit(Surv(Survival.time, Status)~1, data = IPD_REVEL_OS)
lines(REVEL_kmfit, xlim = c(0, 21), col = "green", conf.int = F)
REVEL_myfit <- survreg(Surv(Survival.time, Status)~1, data = IPD_REVEL_OS, dist = "exponential")
REVEL_lambda <- exp(-REVEL_myfit$coefficients)
surv_probs <- exp(-REVEL_lambda*time_seq)
#lines(time_seq, surv_probs, col= "green")

legend("topright", legend = c("INTEREST", "ZODIAC", "REVEL"), col = c("blue", "red", "green"), lty = 1)

med_INTEREST <- summary(INTEREST_kmfit)$table["median"]
med_ZODIAC <- summary(ZODIAC_kmfit)$table["median"]
med_REVEL <- summary(REVEL_kmfit)$table["median"]



lam <- c(log(2)/med_INTEREST, log(2)/med_ZODIAC, log(2)/med_REVEL) 
dc <- c(sum(IPD_INTEREST_OS$Status==1), sum(IPD_ZODIAC_OS$Status==1), sum(IPD_REVEL_OS$Status==1))
hist_data <- list(H = length(lam),                  
                dc = dc,                
                Tc = dc/lam,            
                prior_prec_tau = 1)     

# function to generate inital values for MCMC
hist_init <- function() {
  list(log_lam = rnorm(hist_data$H, mean = 0, sd = 0.5),
       mu = rnorm(1, mean = 1, sd = 0.5),
       tau = rexp(1, rate = 3)  
  )  
}

# Call WinBUGS
model <- jags.model("Thesis/Chapter 3/jags_model.txt",
                    data = hist_data,
                    inits = hist_init,
                    n.chains = 5,
                    n.adapt = 1000)

update(model, 1000) 

params <- c("lam_pred")


samples <- coda.samples(model,
                        variable.names = params,
                        n.iter = 400000,
                        thin = 20)

library(coda)
combined_samples <- do.call(rbind, samples)
lam_pred_all <- as.numeric(combined_samples[,1])


# Plot histogram
hist(lam_pred_all,
     breaks = 5000,
     main = "",
     xlim = c(0, 0.2),
     freq = F,
     ylim = c(0, 25),
     xlab = expression("Predicted "*lambda[c])
)

plot(INTEREST_kmfit, xlim = c(0, 21), col = "blue", conf.int = F, ylab = "Overall Survival", xlab = "Time (months)")
lines(ZODIAC_kmfit, xlim = c(0, 21), col = "red", conf.int = F)
lines(REVEL_kmfit, xlim = c(0, 21), col = "green", conf.int = F)
Med <- exp(-quantile(lam_pred_all, 0.5)*time_seq)
LQ <- exp(-quantile(lam_pred_all, 0.05)*time_seq)
UQ <- exp(-quantile(lam_pred_all, 0.95)*time_seq)
lines(time_seq, Med)
lines(time_seq, LQ, lty = 2)
lines(time_seq, UQ, lty = 2)

legend("topright", legend = c("INTEREST", "ZODIAC", "REVEL", "MAP Prior", "90% Interval for MAP Prior"), 
       col = c("blue", "red", "green", "black", "black"), lty = c(1, 1, 1, 1, 2))


MCMC_sample <- read.csv(file = "Thesis/Chapter 3/Data/MCMC_sample.csv")
Assurance_Calc <- calc_dte_assurance(n_c = c(5, 60, 115, 170, 225, 280, 335, 390, 445, 500),
                             n_t = c(5, 60, 115, 170, 225, 280, 335, 390, 445, 500),
                             control_dist = "Exponential",
                             control_parameters = "Distribution",
                             control_distribution = "MCMC",
                             MCMC_sample = MCMC_sample,
                             delay_time_SHELF = SHELF::fitdist(c(3, 4, 5), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                             delay_time_dist = "gamma",
                             post_delay_HR_SHELF = SHELF::fitdist(c(0.55, 0.6, 0.7), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                             post_delay_HR_dist = "gamma",  
                             P_S = 0.9,
                             P_DTE = 0.8,
                             cens_method = "Events",
                             cens_IF = 0.8,
                             rec_method = "power",
                             rec_period = 12,
                             rec_power = 1,
                             analysis_method = "LRT",
                             alternative = "one.sided",
                             alpha = 0.025,
                             nSims = 25)

Power_Calc <- calc_dte_assurance(n_c = c(5, 60, 115, 170, 225, 280, 335, 390, 445, 500),
                              n_t = c(5, 60, 115, 170, 225, 280, 335, 390, 445, 500),
                              control_dist = "Exponential",
                              control_parameters = "Fixed",
                              fixed_parameters_type = "Parameter",
                              lambda_c = 0.077,
                              delay_time_SHELF = SHELF::fitdist(c(3.999, 4, 4.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                              delay_time_dist = "gamma",
                              post_delay_HR_SHELF = SHELF::fitdist(c(0.599, 0.6, 0.601), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                              post_delay_HR_dist = "gamma",  
                              P_S = 1,
                              P_DTE = 1,
                              cens_method = "Events",
                              cens_IF = 0.8,
                              rec_method = "power",
                              rec_period = 12,
                              rec_power = 1,
                              analysis_method = "LRT",
                              alternative = "one.sided",
                              alpha = 0.025,
                              nSims = 25)


Power_ND_Calc <- calc_dte_assurance(n_c = c(5, 60, 115, 170, 225, 280, 335, 390, 445, 500),
                              n_t = c(5, 60, 115, 170, 225, 280, 335, 390, 445, 500),
                              control_dist = "Exponential",
                              control_parameters = "Fixed",
                              fixed_parameters_type = "Parameter",
                              lambda_c = 0.077,
                              delay_time_SHELF = SHELF::fitdist(c(3.999, 4, 4.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                              delay_time_dist = "gamma",
                              post_delay_HR_SHELF = SHELF::fitdist(c(0.599, 0.6, 0.601), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                              post_delay_HR_dist = "gamma",  
                              P_S = 1,
                              P_DTE = 0,
                              cens_method = "Events",
                              cens_IF = 0.8,
                              rec_method = "power",
                              rec_period = 12,
                              rec_power = 1,
                              analysis_method = "LRT",
                              alternative = "one.sided",
                              alpha = 0.025,
                              nSims = 25)

outcome_DF <- data.frame(N = seq(10, 1000, length = 10),
                          Ass = sapply(Assurance_Calc, `[[`, 1),
                          Power = sapply(Power_Calc, `[[`, 1),
                          Power_ND = sapply(Power_ND_Calc, `[[`, 1))

png("PowerAss_Exp_Example.png", units="in", width=10, height=8, res=700)
plot(outcome_DF$N, outcome_DF$Ass, type= "l", ylim = c(0,1),
     xlab = "Number of Patients", ylab = "Assurance/Power")
lines(outcome_DF$N, outcome_DF$Power, type= "l", col = "blue", lty = 2)
lines(outcome_DF$N, outcome_DF$Power_ND, type= "l", col = "red", lty = 3)

legend("bottomright", legend = c("Assurance", "Power", "Power assuming no delay"),
       col = c("black", "blue", "red"), lty = 1:3)
dev.off()


############
#Simplifying invesitgation
###########
DTEDataSetsFunc <- function(author){
  
  if (author=="Brahmer"){
    trialLength <- 15
    THat <- 3
    ylabel <- "Progression free survival (% of patients)"
  } else if (author=="Yen"){
    trialLength <- 30
    THat <- 3.5
    ylabel <- "Overall survival (%)"
  } else if (author=="Borghaei"){
    trialLength <- 60
    THat <- 6
    ylabel <- "Overall survival (%)"
  }
  
  controldata <- read.csv(file = paste0("Thesis/Chapter 3/Data/", author, "/IPD-control.csv"))
  treatmentdata <- read.csv(file = paste0("Thesis/Chapter 3/Data/", author, "/IPD-treatment.csv"))
  
  combinedData <- data.frame(time = c(controldata$Survival.time, treatmentdata$Survival.time), 
                             status = c(controldata$Status, treatmentdata$Status), 
                             group = c(rep("Control", nrow(controldata)), rep("Treatment", nrow(treatmentdata))))
  
  kmfit <- survfit(Surv(time, status)~group, data = combinedData)
  plot(kmfit, conf.int = F, col=c("blue", "red"), xlab = "Time (months)", ylab = ylabel, yaxt = "n",
       cex.axis=1.5, cex.lab=1.5, cex.main=2)
  axis(2, at=seq(0, 1, by=0.2), labels=seq(0, 100, by=20))
  
  
  #Finding the Weibull parameters and plotting the line
  weibfit <- survreg(Surv(time, status)~1, data = combinedData[combinedData$group=="Control",], dist = "weibull")
  gammac <- as.numeric(exp(-weibfit$icoef[2]))
  lambdac <- as.numeric(1/(exp(weibfit$icoef[1])))
  controltime <- seq(0, trialLength, by=0.01)
  controlsurv <- exp(-(lambdac*controltime)^gammac)
  #lines(controltime, controlsurv, col="blue")
  
  #Now we look at the treatment curve
  #Need to find a least squares estimate for this Weibull parameterisation
  
  kmfit <- survfit(Surv(time, status)~1, data = combinedData[combinedData$group=="Treatment",])
  treatmenttime <- seq(0, THat, by=0.01)
  treatmentsurv <- controlsurv <- exp(-(lambdac*treatmenttime)^gammac)
  #lines(treatmenttime, treatmentsurv, col="red")
  
  treatmenttime1 <- seq(THat, trialLength, by=0.01)
  
  
  # #Scenario 1
  # 
  # optimfunc1 <- function(par){
  #   diff <- 0
  #   gammat <- gammac
  #   for (i in 1:length(treatmenttime1)){
  #     y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
  #     treatmentsurv <- exp(-(lambdac*THat)^gammac - par[1]^gammat*(treatmenttime1[i]^gammat-THat^gammat))
  #     diff <- diff + (y-treatmentsurv)^2 
  #   }
  #   return(diff) 
  # }
  # 
  # s1gammat <- gammac
  # s1lambdat <- optimize(optimfunc1, c(0, 2))$minimum
  # treatmentsurv1 <- exp(-(lambdac*THat)^gammac - s1lambdat^s1gammat*(treatmenttime1^s1gammat-THat^s1gammat))
  # lines(treatmenttime1, treatmentsurv1, col="red", lty=1)
  # 
  # 
  # #Scenario 2
  # 
  # 
  # optimfunc2 <- function(par){
  #   diff <- 0
  #   for (i in 1:length(treatmenttime1)){
  #     y <-  kmfit$surv[sum(kmfit$time<treatmenttime1[i])]
  #     treatmentsurv <- exp(-(lambdac*THat)^gammac - par[1]^par[2]*(treatmenttime1[i]^par[2]-THat^par[2]))
  #     diff <- diff + (y-treatmentsurv)^2 
  #   }
  #   return(diff) 
  # }
  # 
  # optimoutput <-  optim(par = c(0.2, 1), fn = optimfunc2)
  # s2lambdat <- optimoutput$par[1]
  # s2gammat <- optimoutput$par[2]
  # treatmentsurv1 <- exp(-(lambdac*THat)^gammac - s2lambdat^s2gammat*(treatmenttime1^s2gammat-THat^s2gammat))
  # lines(treatmenttime1, treatmentsurv1, col="red", lty=2)
  
  #Plotting horizontal lines to indicate the delay
  # y <- seq(0, 1, by=0.01)
  # x <- rep(THat, length(y))
  # lines(x,y, lty = 2)
  # text(x = THat, y = 0.9, labels = paste0("Delay = ", THat, "months"), pos = 4, cex = 1.5)
  
  # legend("topright", legend = c("Control", "Method A", "Method B"), col=c("blue", "red", "red"), lty=c(1, 1,2),
  #        cex = 1.5)
  
  legend("topright", legend = c("Control", "Treatment"), col=c("blue", "red"), lty=c(1, 1),
         cex = 1.5)
  
  #Now look at calculating power under the two different scenarios
  
  # powerfunc <- function(n, lambdat, gammat){
  #   powervec <- rep(NA, 1000)
  #   for (i in 1:length(powervec)){
  #     #Simulating control times
  #     u <- runif(n)
  #     controltimes <- (1/lambdac)*(-log(u))^(1/gammac)
  # 
  #     #Simulating treatment times
  #     u <- runif(n)
  #     CP <- exp(-(lambdac*THat)^gammac)
  #     treatmenttimes <- ifelse(u>=CP, (1/lambdac)*(-log(u))^(1/gammac), (1/(lambdat^gammat)*((lambdat*THat)^gammat-log(u)-(lambdac*THat)^gammac))^(1/gammat))
  # 
  #     #combining control and treatment
  #     combinedDataPower <- data.frame(time = c(controltimes, treatmenttimes),
  #                                     status = rep(1, 2*n), group = c(rep("Control", n), rep("Treatment", n)))
  # 
  #     test <- survdiff(Surv(time, status)~group, data = combinedDataPower)
  #     
  #     
  #     #If the p-value of the test is less than 0.05 then assvec = 1, 0 otherwise
  #     powervec[i] <- test$chisq > qchisq(0.95, 1)
  #   }
  # 
  #   return(mean(powervec))
  # }
  
  
  # nvec <- round(seq(10, 500, length = 30))
  # powervecs1 <- mapply(powerfunc, nvec, s1lambdat, s1gammat)
  # powervecs2 <- mapply(powerfunc, nvec, s2lambdat, s2gammat)
  # 
  # powers1smooth <- loess(powervecs1~nvec)
  #  powers2smooth <- loess(powervecs2~nvec)
  # 
  #  plot(nvec*2, predict(powers1smooth), type="l", col="red", ylim=c(0, 1), xlab = "Total sample size", ylab = "Power")
  #  lines(nvec*2, predict(powers2smooth), col="blue", lty=2)
  # 
  # legend("bottomright", legend = c("Method A", "Method B"), col = c("red", "blue"), lty=1:2)
}

png("KM_All_3.png", units="in", width=15, height=5, res=700)
par(mfrow=c(1,3))
DTEDataSetsFunc("Brahmer")
DTEDataSetsFunc("Yen")
DTEDataSetsFunc("Borghaei")
dev.off()








