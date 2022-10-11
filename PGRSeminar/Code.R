controldata <- rnorm(24, mean = 87, sd = 4.75)
treatmentdata <- rnorm(25, mean = 85, sd = 5.2)

t.test(controldata, treatmentdata)


hist(controldata, col=rgb(0,0,1,1/4), breaks = 10, freq = F)
hist(treatmentdata, add=T, col=rgb(1,0,0,1/4), breaks=10, freq=F)

x <- seq(70, 100, by=0.01)
y <- dnorm(x, mean = 87, sd = 5)
plot(x,y, type="l")



ncontrol <- 98
ntreatment <- 98


power.t.test(n = 98, delta = 1.5, sd = 5)

x <- seq(0.5, 3, by=0.01)
y <- dnorm(x, 1.81, sd = 0.413)
plot(x, y, type = "l", xlab="delta", ylab="density")

powervec <- rep(NA, 100000)
for (i in 1:100000){
  controldata <- rnorm(98, mean = 87, sd = 5)
  treatmentdata <- rnorm(98, mean = 85, sd = 5)
  test <- t.test(controldata, treatmentdata)
  powervec[i] <- test$p.value<0.05
}
mean(powervec)

assfunc <- function(nsamples){
  assvec <- rep(NA, 10000)
  for (i in 1:10000){
    controldata <- rnorm(nsamples, mean = 87, sd = 5)
    delta <- rnorm(1, mean = 1.81, sd = 0.413)
    treatmentdata <- rnorm(nsamples, mean = 87 - delta, sd = 5)
    test <- t.test(controldata, treatmentdata)
    assvec[i] <- test$p.value<0.05
  }
  mean(assvec)
}


nsamples <- seq(100, 200, by = 5)
outcomes <-sapply(nsamples, assfunc)
outcomes[8]
nsamples[9]

plot(nsamples, outcomes)
