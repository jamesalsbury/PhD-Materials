library(pwr)

###Plot for the what is power slide?
nvec <- seq(5, 200, by=5)
powervec1 <- rep(NA, length = length(nvec))
powervec2 <- rep(NA, length = length(nvec))
powervec3 <- rep(NA, length = length(nvec))
for (i in 1:length(nvec)){
  powervec1[i] <- pwr.t.test(d = 1, n = nvec[i], sig.level=0.05, type = "two.sample")$power
  powervec2[i] <- pwr.t.test(d = 0.5, n = nvec[i], sig.level=0.05, type = "two.sample")$power
  powervec3[i] <- pwr.t.test(d = 0.2, n = nvec[i], sig.level=0.05, type = "two.sample")$power
}


plot(nvec*2, powervec1, type="l", ylim=c(0,1), xlab="Total sample size", ylab= "Power", font=2, cex= 1.5)
lines(nvec*2, powervec2, col="blue")
lines(nvec*2, powervec3, col="red")


legend("bottomright", legend = c(expression(paste(delta, " = 1.0")), expression(paste(delta, " = 0.5")),
                                 expression(paste(delta, " = 0.2")) ), col=c("black", "blue", "red"), lty=1, cex=1.5)

##Plots for the uncertainty about variances plots
#diastolic blood pressure


x <- seq(30, 150, by = 0.01)
controly <- dnorm(x, 90, sd = 10)
treatmenty <- dnorm(x, 88, sd = sqrt(15.6))
plot(x, controly, type="l", ylim=c(0, 0.25), col="blue", xlab="Diastolic blood pressure", ylab = "Density", font=2, cex=1.5)
lines(x, treatmenty, col="red")
legend("topright", legend = c("Control", "Treatment"), col = c("blue", "red"), lty=1)

notbenefitx <- seq(80, 140, by=0.1)
nottreatmenty <- dnorm(notbenefitx, 88, sd = sqrt(15.6))

polygon(c(notbenefitx[notbenefitx>=89], max(notbenefitx), 89), c(nottreatmenty[notbenefitx>=89], 0, 0), col="green")

axis(1, at = 89,
     labels = 89, col="green", col.axis = "red")




1 - pnorm(89, 88, sd = 15.65)

sdvec <- seq(1, 20, by=0.01)

areaundercurve <- rep(NA, length = length(sdvec))

areaundercurve  <- (1 - pnorm(89, 88, sd = sdvec))



plot(sdvec, areaundercurve, type="l", xlab = expression(paste(sigma[t])), ylab= expression(paste("Area [80, ", infinity, "]")))




diff = abs(areaundercurve-0.3)

which.min(diff)                                                                                          

sdvec[92]




##How much difference does the different variation make in terms of power?

powerfunc <- function(n){
  N <- 1000
  powervec1 <- rep(NA, N)
  powervec2 <- rep(NA, N)
  powervec3 <- rep(NA, N)
  for (i in 1:N){
    controlgroup <- rnorm(n, 90, sd = 10)
    treatmentgroup1 <- rnorm(n, 88, sd = 4)
    treatmentgroup2 <- rnorm(n, 88, sd = 10)
    treatmentgroup3 <- rnorm(n, 88, sd = 17)
    powervec1[i] <- t.test(controlgroup, treatmentgroup1)$p.value<0.05
    powervec2[i] <- t.test(controlgroup, treatmentgroup2)$p.value<0.05
    powervec3[i] <- t.test(controlgroup, treatmentgroup3)$p.value<0.05
  }
  return(list(powervec1 = mean(powervec1), powervec2 = mean(powervec2), powervec3 = mean(powervec3)))
}

ssvec <- seq(10, 500, by=10)
powervec1 <- rep(NA, 50)
powervec2 <- rep(NA, 50)
powervec3 <- rep(NA, 50)

for (i in 1:50){
  powervec1[i] <- powerfunc(ssvec[i])$powervec1
  powervec2[i] <- powerfunc(ssvec[i])$powervec2
  powervec3[i] <- powerfunc(ssvec[i])$powervec3
}

power1smooth <- loess(powervec1~ssvec)
power2smooth <- loess(powervec2~ssvec)
power3smooth <- loess(powervec3~ssvec)

plot(ssvec*2, predict(power1smooth), type="l", ylab = "", font=2, xlab = "")
lines(ssvec*2, predict(power2smooth), col="blue", lty=2)
lines(ssvec*2, predict(power3smooth), col="red", lty=3)


legend("bottomright", legend = c(expression(paste(sigma, " = 4")), expression(paste(sigma, " = 10")),
                                 expression(paste(sigma, " = 17")) ), col=c("black", "blue", "red"), lty=1:3, cex=1)


mtext(side=2, line=3, "Power", font=2, cex=1.2)
mtext(side=1, line=3, "Total sample size", font=2, cex=1.2)




