library(simsurv)
library(survival)
library(truncnorm)

##############################################
#Power with point-estimates of the treatment effects
##############################################
lambda1 <- -log(0.8)/5
lambda2 <- -log(0.6)/5

PowerFunc <- function(n1, n2){
  covs <- data.frame(id = 1:n1)
  power <- rep(NA, 100)
  for (i in 1:100){
    Control <- data.frame(simsurv(dist = "exponential", x = covs, lambdas = lambda2, maxt = 5), group = rep("control", n1))  
    Treatment <- data.frame(simsurv(dist = "exponential", x = covs, lambdas = lambda1, maxt = 5), group = rep("treatment", n2)) 
    
    Combined <- rbind(Control, Treatment)
    
    test <- survdiff(Surv(eventtime, status)~group, data = Combined)
    power[i] <- test$chisq > qchisq(0.95, 1)
  }
  mean(power)
}

N1Vec <- N2Vec <- seq(10, 300, by=10)
powervec <- rep(NA, 30)

for (i in 1:30){
  powervec[i] <- PowerFunc(N1Vec[i], N2Vec[i])
}

lo <- loess(powervec~N1Vec)
plot(N1Vec,predict(lo), type="l", xlab="Total sample size", ylab="Power", xaxt="n", ylim=c(0,1))
axis(1, at = seq(0, 300,by = 50), labels = seq(0, 600, by = 100))

##############################################
#Power with 1 interim analysis - when half the expected events have occurred
##############################################


lambda1 <- -log(0.8)/5
lambda2 <- -log(0.6)/5


AssFunc <- function(n1, n2, percent1){
  #We expect 0.4n1+0.2n2 events to occur
  expect <- 0.4*n1+0.2*n2
  covs <- data.frame(id = 1:n1)
  ass1 <- rep(NA, 100)
  ass2 <- rep(NA, 100)
  ass3 <- rep(NA, 100)
  
  for (i in 1:100){
    Control <- data.frame(simsurv(dist = "exponential", x = covs, lambdas = lambda2, maxt = 5), group = rep("control", n1))  
    Treatment <- data.frame(simsurv(dist = "exponential", x = covs, lambdas = lambda1, maxt = 5), group = rep("treatment", n2)) 
    
    Combined <- rbind(Control, Treatment)
    
    #Get the dataframe of the first percentage of total events only
    Percent1SortedEvents <- Combined[Combined$status==1,][order(Combined[Combined$status==1,]$eventtime),][1:(expect*percent1),]

    #Test the first percentage
    if (sum(Percent1SortedEvents$group=="treatment")==0){
      
    } else if (sum(Percent1SortedEvents$group=="treatment")==expect){
      
    } else {
      test1 <- survdiff(Surv(eventtime, status)~group, data = Percent1SortedEvents)
      
      if (test1$chisq > qchisq(1-0.0054, 1)){
        ass1[i] <- 1
      } else {
        test2 <- survdiff(Surv(eventtime, status)~group, data = Combined)
        if (test2$chisq >qchisq(1-0.0492, 1)){
          ass1[i] <- 1
        } else {
          ass1[i] <- 0
        }
      }
    }
       
  }
  return(mean(na.omit(ass1)))
}


N1Vec <- N2Vec <- seq(20, 300, by=10)
assvec <- rep(NA, length(N1Vec))


for (i in 1:length(N1Vec)){
  assvec[i] <- AssFunc(N1Vec[i], N2Vec[i], 0.5)
}

lo1 <- loess(assvec~N1Vec)
lines(N1Vec,predict(lo1), col="blue")

##Can see that it doesn't really make a difference, having 1 interim analysis here
#Need to write down what I have done and what assumptions I have made; then can look at assurance in this setup




rtruncnorm()





for (i in 1:10000){
  mu <-  0.2
  var <- 0.05
  S1 <- rbeta(1, 30, 70)
  rho <- rtruncnorm(1, a = -S1, b = 1-S1, mean = mu, sd = sqrt(var))
  S2 <- S1 + rho
  lambda1 <- -log(S2)/12
  lambda2 <- -log(S1)/12
}








