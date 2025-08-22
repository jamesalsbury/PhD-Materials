probs <- c(0.05, 0.05,
           0.15, 0.15, 0.15, 0.15,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.35, 0.35, 0.35,
           0.45, 0.45,
           0.55, 0.55, 0.65)
png("Trial_Roulette_Dist.png", units="in", width=10, height=6, res=700)
par(mfrow = c(1,2))
hist(probs, ylim = c(0, 4), xlim = c(0,1), col = "orange",
     xlab = expression("Probability ("*theta*")"),
     main = "", freq=F, ylab = "", yaxt="n", cex.axis=1.5, cex.lab=1.5, cex.main=2)
x <- seq(0, 1, by=0.01)
lines(x, dbeta(x, 2.04, 4.89), col = "red", lwd = 2)
legend("topright", legend = "Beta(2.04, 4.89)", lty = 1, col = "red", lwd = 2, cex = 1.5)

xAxis <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
CDF <- c(0.1, 0.3, 0.6, 0.75, 0.85, 0.95)
plot(xAxis, CDF, xlab = expression(theta), xlim = c(0,1), ylim = c(0,1),
     cex.axis=1.5, cex.lab=1.5, cex.main=2, pch = 16, cex = 2)
points(c(0, 1), c(0,1), pch = 1, cex = 2)
lines(x, pbeta(x, 2.04, 4.89))
dev.off()


png("Math_Pooling.png", units="in", width=12, height=6, res=700)
par(mfrow = c(1,2))
x <- seq(0, 1.2, by=0.01)
f1 <- dnorm(x, 0.4, 0.08)
f2 <- dnorm(x, 0.6, 0.08)
w1 <- 0.5
w2 <- 1 - w1

#Plotting the linear pool
plot(x, f1, type = "l",
     xlab = expression("Probability ("*theta*")"), col = "red",
     ylab = "Density", xlim = c(0, 1.2),
     cex.axis=1.5, cex.lab=1.5, cex.main=2)
lines(x, f2, col = "blue")
lines(x, w1*f1 + w2*f2, lty = 2, lwd = 2)
legend("topright", legend = c("Expert 1", "Expert 2", "Linear Pool"), col = c("red", "blue", "black"),
       lty = c(1, 1, 2), lwd = c(1, 1, 2))

#Plotting the log pool
log_pool_raw <- f1^w1 * f2^w2
log_pool <- log_pool_raw / sum(log_pool_raw * diff(x)[1])
plot(x, f1, type = "l",
     xlab = expression("Probability ("*theta*")"), col = "red",
     ylab = "Density", xlim = c(0, 1.2),
     cex.axis=1.5, cex.lab=1.5, cex.main=2)
lines(x, f2, col = "blue")
lines(x, log_pool, lty = 2, lwd = 2)
legend("topright", legend = c("Expert 1", "Expert 2", "Log Pool"), col = c("red", "blue", "black"),
       lty = c(1, 1, 2), lwd = c(1, 1, 2))
dev.off()

png("Elicit_Rio.png", units="in", width=10, height=6, res=700)
x <- seq(0, 7, by=0.01)
ex1 <- dgamma(x, 11.4, 4.45)
ex2 <- dgamma(x, 10.7, 2.87)
ex3 <- dgamma(x, 5.68, 1.9)
ex4 <- dgamma(x, 4.09, 1.28)
w1 <- 0.25
w2 <- 0.25
w3 <- 0.25
w4 <- 0.25
RIO <- dlnorm(x, 1.16, 0.437)
plot(x, ex1, type = "l", ylab = "Density", xlab = expression(theta),
     col = "red", lty = 2,
     cex.axis=1.5, cex.lab=1.5, cex.main=2)
lines(x, ex2, col = "blue", lty = 2)
lines(x, ex3, col = "green", lty = 2)
lines(x, ex4, col = "purple", lty = 2)
lines(x, RIO, lwd = 2)
lines(x, w1*ex1 + w2*ex2 + w3*ex3 + w4*ex4, lty = 2, lwd = 2)
legend("topright", legend = c("Expert 1", "Expert 2", "Expert 3", "Expert 4", "RIO"),
       lty = c(2, 2, 2, 2, 1), col = c("red", "blue", "green", "purple", "black"),
       lwd = c(1, 1, 1, 1, 2))
dev.off()
