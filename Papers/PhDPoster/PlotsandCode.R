png("PowerAss.png", units = "in", width = 8, height = 5, res = 700)

muc <- 87
mut <- 85
sigma <- 5

powerFunc <- function(nc, nt){
  powervec <- rep(NA, 1000)
  for (i in 1:1000){
    controlData <- rnorm(nc, muc, sd = sigma)
    treatmentData <- rnorm(nt, mut, sd = sigma)
    test <- t.test(controlData, treatmentData)
    powervec[i] <- test$p.value < 0.05
  }
  return(mean(powervec))
}



ncvec <- seq(2, 200, by=5)
ntvec <- seq(2, 200, by=5)

powervec <- mapply(FUN = powerFunc, nc = ncvec, nt = ntvec)

smoothpower <- loess(powervec~ncvec)

plot(2*ncvec, predict(smoothpower), type="l", xlab = "Total sample size",
     ylab = "Power/Assurance", ylim=c(0,1))

assFunc <- function(nc, nt){
  assvec <- rep(NA, 1000)
  for (i in 1:1000){
    controlData <- rnorm(nc, muc, sd = sigma)
    theta <- rnorm(1, 1.82, sd = 0.577)
    treatmentData <- rnorm(nt, muc+theta, sd = sigma)
    test <- t.test(controlData, treatmentData)
    assvec[i] <- test$p.value < 0.05
  }
  return(mean(assvec))
}

assvec <- mapply(FUN = assFunc, nc = ncvec, nt = ntvec)

smoothass <- loess(assvec~ncvec)

lines(2*ncvec, predict(smoothass), col = "red", lty = 2)


legend("bottomright", legend = c("Power", "Assurance"), col = c("black", "red"), lty = 1:2)
dev.off()

predict(smoothass, newdata = 98)


predict(smoothpower, newdata = 98)





png("DTEKM.png", units="in", width=8, height=5, res=700)

controldata <- read.csv(file = "Papers/DTEPaper/data/Brahmer/IPD-control.csv")
treatmentdata <- read.csv(file = "Papers/DTEPaper/data/Brahmer/IPD-treatment.csv")

combinedData <- data.frame(time = c(controldata$Survival.time, treatmentdata$Survival.time), 
                           status = c(controldata$Status, treatmentdata$Status), 
                           group = c(rep("Control", nrow(controldata)), rep("Treatment", nrow(treatmentdata))))

kmfit <- survfit(Surv(time, status)~group, data = combinedData)
plot(kmfit, conf.int = F, col=c("blue", "red"), xlab = "Time (months)", 
     ylab = "Progression free survival (% of patients)", xaxt = "n", yaxt = "n")
axis(1, at=seq(0, 21, by=3), labels=seq(0, 21, by=3))
axis(2, at=seq(0, 1, by=0.1), labels=seq(0, 100, by=10))
legend("topright", legend = c("Docetaxel", "Nivolumab"), col = c("blue", "red"), lty=1)

dev.off()









