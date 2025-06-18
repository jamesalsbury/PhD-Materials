sigma <- 10
delta_0 <- 5
alpha <- 0.05
beta <- 0.2


n1 <- 2*sigma^2*(qnorm(alpha, lower.tail = F) + qnorm(beta, lower.tail = F))^2
n1 <- n1/(delta_0)^2

n_0 <- 4
n <- 1:500
power <- pnorm(qnorm(alpha) + (delta_0*sqrt(n))/sqrt(2*sigma^2))
assurance <- pnorm(sqrt(n_0/(n+n_0))*(qnorm(alpha)+(delta_0*sqrt(n))/sqrt(2*sigma^2)))
assurance_upper_bound <- pnorm(sqrt(n_0)*delta_0/sqrt(2*sigma^2))
Norm_Ass <- assurance/assurance_upper_bound



png("PowerAssurance.png", units="in", width=10, height=6, res=700)
plot(n, power, type = "l", col = "blue",
     xlab = "Number of patients (in each group)", ylab = "Power/Assurance", 
     cex.axis=1.5, cex.lab=1.5, cex.main=2)
lines(n, assurance, col = "red", lty = 2)
abline(h = assurance_upper_bound, lty = 3)
legend("bottomright", legend = c("Power", "Assurance", "Assurance bound"), 
       lty = 1:3, col = c("blue", "red", "black"), cex = 1.25)
dev.off()

png("Norm_Assurance.png", units="in", width=10, height=6, res=700)
plot(n, Norm_Ass, type = "l", ylim = c(0,1), xlab = "Number of patients (in each group)",
     ylab = "Normalised Assurance",
     cex.axis=1.5, cex.lab=1.5, cex.main=2)
dev.off()

png("Pred_delta_hat.png", units="in", width=10, height=8, res=1000)
par(mfrow = c(2,2))

produce_plot <- function(n){
  critValue <- qnorm(0.95) * sqrt(2*sigma^2) / sqrt(n)
  delta_hat <- seq(-20, 30, by = 0.1)
  
  plot(delta_hat,
       dnorm(delta_hat, mean = delta_0, sd = sqrt(2*sigma^2 * (1 / n))),
       type = "l", col = "blue",
       main = bquote("Predictive " * hat(delta) * " for " * n * " = " * .(n)),
       ylab = "Density",
       xlab = expression(hat(delta)),
       ylim = c(0, 0.35),
       cex.axis=1.5, cex.lab=1.5, cex.main=2)
  
  lines(delta_hat,
        dnorm(delta_hat, mean = delta_0, sd = sqrt(2*sigma^2 * (1 / n_0 + 1 / n))),
        col = "red")
  
  abline(v = critValue, lwd = 1, col = "green")
  
  legend("topright",
         legend = c("Predictive (Power)", "Predictive (Assurance)", "Critical Value"),
         col = c("blue", "red", "green"), lty = 1)
}

produce_plot(5)
produce_plot(25)
produce_plot(50)
produce_plot(75)
dev.off()


sim_data_power <- function(n1){
  placebo <- rbinom(1, n1, prob = 0.45)
  moxo <- rbinom(1, n1, prob = 0.3)
  M <- as.table(rbind(c(moxo, n1-moxo), c(placebo, n1-placebo)))
  dimnames(M) <- list(group = c("M", "P"),
                      Outcome = c("Raised","Not_Raised"))
  test <- prop.test(M, correct = F)
  test$p.value < 0.05
}

sim_data_ass <- function(n1, m, nu){
  theta_c <- rbeta(1, 10.7, 13.1)
  rho <- truncnorm::rtruncnorm(1, mean = m, sd = sqrt(nu), a = (theta_c-1), b = theta_c)
  theta_t <- theta_c - rho #Treatment
  placebo <- rbinom(1, n1, prob = theta_c)
  moxo <- rbinom(1, n1, prob = theta_t)
  M <- as.table(rbind(c(moxo, n1-moxo), c(placebo, n1-placebo)))
  dimnames(M) <- list(group = c("M", "P"),
                      Outcome = c("Raised","Not_Raised"))
  test <- prop.test(M, correct = F)
  test$p.value < 0.05
}

png("Power_Assurance_Moxo.png", units="in", width=10, height=6, res=700)
n1 <- 10:500
power_estimates <- sapply(n1, function(n) mean(replicate(500, sim_data_power(n))))
smooth_power <- loess(power_estimates~n1)
plot(n1, predict(smooth_power, newdata = n1), type = "l",
     ylim = c(0,1), xlim = c(0, 500),
     ylab = "Power/Assurance",
     xlab = "Number of Patients (in each group)")


ass_estimates1 <- sapply(n1, function(n) mean(replicate(500, sim_data_ass(n, 0.15, 0.0001))))
smooth_ass1 <- loess(ass_estimates1~n1)
lines(n1, predict(smooth_ass1, newdata = n1), col = "blue", lty = 2)


ass_estimates2 <- sapply(n1, function(n) mean(replicate(500, sim_data_ass(n, 0.15, 0.01))))
smooth_ass2 <- loess(ass_estimates2~n1)
lines(n1, predict(smooth_ass2, newdata = n1), col = "red", lty = 3)


ass_estimates3 <- sapply(n1, function(n) mean(replicate(500, sim_data_ass(n, 0.1, 0.01))))
smooth_ass3 <- loess(ass_estimates3~n1)
lines(n1, predict(smooth_ass3, newdata = n1), col = "green", lty = 4)

legend("bottomright", legend = c("Power", "Assurance: Scenario 1", "Assurance: Scenario 2", "Assurance: Scenario 3"),
       lty = 1:4, col = c("black", "blue", "red", "green"))
dev.off()



probs <- c(0.05, 0.05,
           0.15, 0.15, 0.15, 0.15,
           0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
           0.35, 0.35, 0.35,
           0.45, 0.45,
           0.55, 0.55, 0.65)
png("Trial_Roulette_Dist.png", units="in", width=10, height=6, res=700)
hist(probs, ylim = c(0, 4), xlim = c(0,1), col = "orange",
     xlab = expression("Probability ("*theta*")"),
     main = "", freq=F, ylab = "", yaxt="n")
x <- seq(0, 1, by=0.01)
lines(x, dbeta(x, 2.04, 4.89), col = "red", lwd = 2)
legend("topright", legend = "Be(2.04, 4.89)", lty = 1, col = "red", lwd = 2)
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
     ylab = "Density", xlim = c(0, 1.2))
lines(x, f2, col = "blue")
lines(x, w1*f1 + w2*f2, lty = 2, lwd = 2)
legend("topright", legend = c("Expert 1", "Expert 2", "Linear Pool"), col = c("red", "blue", "black"),
       lty = c(1, 1, 2), lwd = c(1, 1, 2))

#Plotting the log pool
log_pool_raw <- f1^w1 * f2^w2
log_pool <- log_pool_raw / sum(log_pool_raw * diff(x)[1])
plot(x, f1, type = "l",
     xlab = expression("Probability ("*theta*")"), col = "red",
     ylab = "Density", xlim = c(0, 1.2))
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
     col = "red", lty = 2)
lines(x, ex2, col = "blue", lty = 2)
lines(x, ex3, col = "green", lty = 2)
lines(x, ex4, col = "purple", lty = 2)
lines(x, RIO, lwd = 2)
lines(x, w1*ex1 + w2*ex2 + w3*ex3 + w4*ex4, lty = 2, lwd = 2)
legend("topright", legend = c("Expert 1", "Expert 2", "Expert 3", "Expert 4", "RIO"),
       lty = c(2, 2, 2, 2, 1), col = c("red", "blue", "green", "purple", "black"),
       lwd = c(1, 1, 1, 1, 2))
dev.off()