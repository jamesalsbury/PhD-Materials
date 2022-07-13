#Libraries
library(plyr)
library(survival)

#We have the following parameters for the control
lambda2 <- 0.08
gamma2 <- 0.8

#We draw the control group
controltime <- seq(0, 100, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)
plot(controltime, controlsurv, type="l", col="blue", ylim=c(0,1), xlim=c(0, 20))

#We have the following parameters for the treatment
gamma1 <- 0.8

#We have elicited the following distributions
#bigT <- normal(6, 1)
#HR <- beta(10, 6)

#Therefore, using the mean values from the above distribtuions, we have the following treatment curve
bigTMean <- 6
HRMean <- 10/16

lambda1 <- exp((log(HRMean)/gamma2)+log(lambda2))
treatmenttime <- seq(bigTMean, 100 , by=0.01)
treatmentsurv <- exp(-(lambda2*bigTMean)^gamma2-lambda1^gamma1*(treatmenttime^gamma1-bigTMean^gamma1))
lines(treatmenttime, treatmentsurv, type="l", col="red")


#So this is what we are planning our trial on

#We can use this to calculate the sample size needed for 80% assurance (power)

#We assume the maximum length of trial will be 60 months
trialLength <- 60

assfunc <- function(n1, n2){
  assvec <- rep(NA, 250)
  eventvec <- rep(NA, 250)
  for (i in 1:length(assvec)){
    bigT <- rnorm(1, 6, 1)
    HR <- rbeta(1, 10, 6)
    lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
    
    controldata <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))
    
    CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
    u <- runif(n2)
    
    suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
    
    DataCombined <- data.frame(time = c(controldata$time, z),
                               group = c(rep("Control", n1), rep("Treatment", n2)))
    
    DataCombined$cens <- DataCombined$time < trialLength
    
    DataCombined$cens <- DataCombined$cens*1
    
    
    
    if (sum(DataCombined$cens)==(n1+n2)){
      
    } else {
      DataCombined[DataCombined$cens==0,]$time <- trialLength
    }
    
    
    eventvec <- sum(DataCombined$cens==1)
    
    test <- survdiff(Surv(time, cens)~group, data = DataCombined)
    assvec[i] <- test$chisq > qchisq(0.95, 1)
  }
  return(list(assvec = mean(assvec), eventvec=mean(na.omit(eventvec))))
}


n1vec <- seq(10, 500, by=20)
n2vec <- seq(10, 500, by=20)
assvec <- rep(NA, length(n1vec))
eventvec <- rep(NA, length(n1vec))
for (j in 1:length(n1vec)){
  assvec[j] <- assfunc(n1vec[j], n2vec[j])$assvec
  eventvec[j] <- assfunc(n1vec[j], n2vec[j])$eventvec
}

nvec <- n1vec+n2vec
asssmooth <- loess(assvec~nvec)
par(mar = c(10, 10, 10, 10))
plot(nvec, predict(asssmooth), ylim=c(0,1), ylab = "Assurance", xlab="Total sample size", type="l")

eventsmooth <- loess(eventvec~nvec)

#How many patients do we need for 80% assurance?


for (i in 20:1000){
  if (predict(asssmooth, newdata = i)>0.8){
    break
  }
}

npatients <- round_any(i, 2)

#How many events do we expect to see with this number of patients?

events <- round(predict(eventsmooth, newdata = npatients))

#We need to determine when is best to perform any IA

#We require 458 events for 80% assurance, 496 patients

#Lets simulate some data, we look at 50% through the information fraction, so we look when 229 events have happened
#We simulate the delay to be 8 months, but keep the HR the same
#We can combine the data with the prior to compute BPP


n1 <- npatients/2
n2 <- npatients/2

lambda2 <- 0.08
gamma2 <- 0.8
gamma1 <- 0.8
HR <- 0.625
lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
bigT <- 3

controldata <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))
CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
u <- runif(n2)

suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))

DataCombined <- data.frame(time = c(controldata$time, z),
                           group = c(rep("Control", n1), rep("Treatment", n2)))


IATime <- DataCombined[order(DataCombined$time),][events*0.75,]$time

#We are performing tha IA at 75% through thr IF

#What will the data set look like here?

DataCombined$cens <- DataCombined$time < IATime

OriginalData <- DataCombined

for (i in 1:nrow(DataCombined)){
  if (DataCombined$time[i]>IATime){
    DataCombined$time[i]=IATime
  }
}

#Need to work out what direction to go down
if (bigT<bigTMean){
  
}



#How to estimate the HR

#We need to use the control data to esimtate lambda2 and gamma2
controlSample <- read_excel(chosenFile$datapath, sheet=1)
weibfit <- survreg(Surv(time, cens)~1, data = controlSample, dist = "weibull")
updateTextInput(session, "lambda2", value = round(as.numeric(1/(exp(weibfit$icoef[1]))), 3))
updateTextInput(session, "gamma2", value = round(as.numeric(exp(-weibfit$icoef[2])), 3))


















