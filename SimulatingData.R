lambda2 <- 0.08
gamma2 <- 1.5
n1 <- 150
bigT <- 3
n2 <- 150

controldata <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))


lambda1 <- 0.08
gamma1 <- 0.6


CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
u <- runif(n2)

suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))

treatmentdata <- data.frame(time = z)


combinedData <- data.frame(time = c(controldata$time, treatmentdata$time), cens = rep(1, n1+n2), group = c(rep("control", n1), rep("treatment", n2)))

saveRDS(combinedData, file = "simData5.rds")


controldata <- simData5[simData5$group=="control",]

weibcontrolfit <- survreg(Surv(time, cens)~1, data = controldata, dist = "weibull")
kmcontrol <- survfit(Surv(time, cens)~1, data = controldata)

gamma2 <- as.numeric(exp(-weibcontrolfit$icoef[2]))
lambda2 <- as.numeric(1/(exp(weibcontrolfit$icoef[1])))

plot(kmcontrol, conf.int = F)

controltime <- seq(0, 50, by=0.01)
controlsurv <- exp(-(lambda2*controltime)^gamma2)
lines(controltime, controlsurv)


treatmentdata <- simData5[simData5$group=="treatment",]

kmtreatment <- survfit(Surv(time, cens)~1, data = treatmentdata)

lines(kmtreatment, conf.int = F)


optimfunc1 <- function(par){
  
  gamma1 <- gamma2
  se <- 0
  
  for (i in 1:150){
    if (kmtreatment$time[i]<3){
      y <- exp(-(lambda2*kmtreatment$time[i])^gamma2)
      se <- se + (y-kmtreatment$surv[i])^2
    } else {
      y <- exp(-(lambda2*3)^gamma2-par^gamma1*(kmtreatment$time[i]^gamma1-3^gamma1))
      se <- se + (y-kmtreatment$surv[i])^2
    }
  }
  return(se)
  
}

 optimoutput1 <-  optimize(f = optimfunc1, interval = c(0,1))


lambda1 <- optimoutput$minimum



optimfunc2 <- function(par){
  
  se <- 0
  
  for (i in 1:150){
    if (kmtreatment$time[i]<3){
      y <- exp(-(lambda2*kmtreatment$time[i])^gamma2)
      se <- se + (y-kmtreatment$surv[i])^2
    } else {
      y <- exp(-(lambda2*3)^gamma2-par[1]^par[2]*(kmtreatment$time[i]^par[2]-3^par[2]))
      se <- se + (y-kmtreatment$surv[i])^2
    }
  }
  return(se)
  
}


optimoutput2 <-  optim(par = c(0, 1), fn = optimfunc2)
optimoutput2





