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


##Plots for the uncertainty about variances plots
#diastolic blood pressure


x <- seq(30, 150, by = 0.1)
controly <- dnorm(x, 90, sd = 10)
treatmenty <- dnorm(x, 75, sd = 9.53)
plot(x, controly, type="l", ylim=c(0, 0.1), col="blue", xlab="Diastolic blood pressure", ylab = "Density")
lines(x, treatmenty, col="red")
legend("topright", legend = c("Control", "Treatment"), col = c("blue", "red"), lty=1)

notbenefitx <- seq(80, 140, by=0.1)
nottreatmenty <- dnorm(notbenefitx, 75, sd = 9.53)

polygon(c(notbenefitx[notbenefitx>=80], max(notbenefitx), 80), c(nottreatmenty[notbenefitx>=80], 0, 0), col="green")

1 - pnorm(80, 75, sd = 17)

sdvec <- seq(1, 20, by=0.01)

areaundercurve <- rep(NA, length = length(sdvec))

areaundercurve  <- (1 - pnorm(80, 75, sd = sdvec))


plot(sdvec, areaundercurve, type="l", xlab = expression(paste(sigma[t])), ylab= expression(paste("Area [80, ", infinity, "]")))


areaundercurve==0.2

diff = abs(areaundercurve-0.3)

which.min(diff)                                                                                          

sdvec[854]

legend("bottomright", legend = c(expression(paste(delta, " = 1.0")), expression(paste(delta, " = 0.5")),
                                 expression(paste(delta, " = 0.2")) ), col=c("black", "blue", "red"), lty=1, cex=1.5)

