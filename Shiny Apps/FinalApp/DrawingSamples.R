library(dplyr)
library(survival)
bigT <- rgamma(1, 4.09, 1.28)
HR <- rbeta(1, 5.82, 5.82)
lambda2 <- 0.2
gamma1 <- gamma2 <- 0.7
lambda1 <- exp((log(HR)/gamma2)+log(lambda2))

controltime <- seq(0, 50, by=0.1)
controlsurv <- exp(-(lambda2*controltime)^gamma2)

plot(controltime, controlsurv, type="l", col="blue")

treatmenttime <- seq(0, bigT, by=0.1)
treatmentsurv <- exp(-(lambda2*treatmenttime)^gamma2)

treatmenttime1 <- seq(bigT, 50, by=0.1)
treatmentsurv1 <- exp(-(lambda2*bigT)^gamma2 - lambda1^gamma1*(treatmenttime1^gamma1-bigT^gamma1))

lines(treatmenttime, treatmentsurv)
lines(treatmenttime1, treatmentsurv1, col="red")

#Generate data according to this scheme

controldata <- data.frame(time = rweibull(10000, gamma2, 1/lambda2), cens = rep(1, 10000))

controlfit <- survfit(Surv(time, cens)~1, data = controldata)

lines(controlfit, conf.int = F, col="green")

CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
u <- runif(10000)
suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))

treatmentdata <- data.frame(time = z, cens = rep(1, 10000))

treatmentfit <- survfit(Surv(time, cens)~1, data = treatmentdata)

lines(treatmentfit, conf.int = F)


DataCombined <- data.frame(time = c(controldata$time, treatmentdata$time),
                           group = c(rep("Control", 10000), rep("Treatment", 10000)))

df <- DataCombined %>%
  filter(time<24)


combinedfit <- survfit(Surv(time)~group, data = df)

#plot(combinedfit)

survdiff(Surv(time)~group, data = df)
