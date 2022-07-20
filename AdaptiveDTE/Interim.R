#Libraries
library(plyr)
library(survival)
library(tidyverse)


#We change the underlying true values of T and HR and see how that affects our judgements

#We will keep HR fixed for now, at 0.625  
#We will fix T at 4.04

lambda2 <- 0.08
gamma2 <- 0.8
bigT <- 4.04
gamma1 <- 0.8
HR <- 0.625
lambda1 <- exp((log(HR)/gamma2)+log(lambda2))
totalLength <- 60

n1 <- 490/2
n2 <- 490/2

IAFunc <- function(){
  IAVec <- rep(NA, 50)
  for (i in 1:50){
    
    #Generate the initial data
    CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
    u <- runif(n2)
    
    suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
    
    DataCombined <- data.frame(time = c(rweibull(n1, gamma2, 1/lambda2),z),
                               group = c(rep("Control", n1), rep("Treatment", n2)))
    
    

    DataCombined <- DataCombined[order(DataCombined$time),]
    
    IATime <- DataCombined[345,]$time
    
    DataCombined$event <- DataCombined$time<IATime
    
    DataCombined[DataCombined$event==0,]$time <- IATime
    
    
    #Extract the control group
    controlSample <- DataCombined %>%
      filter(group=="Control")
    weibfit <- survreg(Surv(time, event)~1, data = controlSample, dist = "weibull")
    #Estimate gamma2 and lambda2 from the control group
    lambda2est <- as.numeric(1/(exp(weibfit$icoef[1])))
    gamma2est <- as.numeric(exp(-weibfit$icoef[2]))
    
    
    #We can use these estimated parameters to simulate the remaining control data
    
    controleventsremaining <- n1 - sum(controlSample$time<IATime)
    
    u <- runif(controleventsremaining, 0, controleventsremaining/n1)
    
    remainingcontroltimes <- exp((1/gamma2est)*log(-log(u))-log(lambda2est))
    
    finalcontrolsample <- data.frame(time = c(controlSample$time[1:(sum(controlSample$time<IATime))], remainingcontroltimes), group = rep("Control", n1), event = rep(0, n1))
    
    finalcontrolsample$event <- finalcontrolsample$time<totalLength
  
    #We now need to estimate the remaining treatment data 
    gamma1est <- gamma2est
    
    treatmentSample <- DataCombined %>%
      filter(group=="Treatment")
    
    kmtreatment <- survfit(Surv(time, event)~1, data = treatmentSample)
    
    estT <- seq(bigT, IATime, by = 0.5)
    
    findlambda1 <- function(par){
      diff <- 0
      for (i in 1:length(estT)){
        diff <- diff + (summary(kmtreatment,time=estT[i])$surv - exp(-(lambda2est*bigT)^gamma2est - par[1]^gamma1est*(estT[i]^gamma1est-bigT^gamma1est)))^2
      }
      return(diff)
    }
    
    lambda1est <-  optimise(findlambda1, c(0,0.2))$minimum
    
    treatmenteventsremaining <- n2 - sum(treatmentSample$time<IATime)
    
    u <- runif(treatmenteventsremaining, 0, treatmenteventsremaining/n2)
    
    remainingtreatmenttimes <- exp((1/gamma1est)*log(1/(lambda1est^gamma1est)*(-log(u)-(lambda2est*bigT)^gamma2est+lambda1est^gamma1est*bigT*gamma1est)))
    
    finaltreatmentsample <- data.frame(time = c(treatmentSample$time[1:(sum(treatmentSample$time<IATime))], remainingtreatmenttimes), group = rep("Treatment", n2), event = rep(0, n2))
    
    finaltreatmentsample$event <- finaltreatmentsample$time<totalLength
    

    #Combine the two groups
    finalDataCombined <- rbind(finalcontrolsample, finaltreatmentsample)
    test <- survdiff(Surv(time, event)~group, data = finalDataCombined)
    
    IAVec[i] <- test$chisq > qchisq(0.95, 1)
    
  }
  return(mean(IAVec))
}

y <- IAFunc()





















