sigma <- 10
delta_0 <- 5
alpha <- 0.025
beta <- 0.2


n1 <- 2*sigma^2*(qnorm(alpha) + qnorm(beta))^2
n1 <- n1/(delta_0)^2

n_0 <- 4
n <- 1:500
power <- pnorm(qnorm(alpha) + (delta_0*sqrt(n))/sqrt(2*sigma^2))
assurance <- pnorm(sqrt(n_0/(n+n_0))*(qnorm(alpha)+(delta_0*sqrt(n))/sqrt(2*sigma^2)))
assurance_upper_bound <- pnorm(sqrt(n_0)*delta_0/sqrt(2*sigma^2))
Norm_Ass <- assurance/assurance_upper_bound



png("PowerAssurance.png", units="in", width=10, height=6, res=700)
plot(n, power, type = "l", col = "blue",
     xlab = "Number of patients (in each group)", ylab = "Power or Assurance", 
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




simulate_data_power <- function(n) {
    r_c <- rbinom(1, n, 0.45)
    r_t <- rbinom(1, n, 0.3)
    p_hat_c <- r_c / n
    p_hat_t <- r_t / n
    p_pool <- (r_c + r_t) / (2 * n)
    se <- sqrt(2 * p_pool * (1 - p_pool) / n)
    z <- (p_hat_c - p_hat_t) / se
    z_crit <- qnorm(1 - 0.025)
    z > z_crit
}



sim_data_ass <- function(n, m, nu){
  
  theta_c <- rbeta(1, 10.7, 13.1)
  rho <- truncnorm::rtruncnorm(1, mean = m, sd = sqrt(nu), a = (theta_c-1), b = theta_c)
  theta_t <- theta_c - rho #Treatment
  r_c <- rbinom(1, n, prob = theta_c)
  r_t <- rbinom(1, n, prob = theta_t)
  
  p_hat_c <- r_c / n
  p_hat_t <- r_t / n
  p_pool <- (r_c + r_t) / (2 * n)
  se <- sqrt(2 * p_pool * (1 - p_pool) / n)
  z <- (p_hat_c - p_hat_t) / se
  z_crit <- qnorm(1 - 0.025)
  z > z_crit
}

png("Power_Assurance_Moxo.png", units="in", width=10, height=6, res=700)
n1 <- 10:500
power_estimates <- sapply(n1, function(n) mean(replicate(500, simulate_data_power(n))))
smooth_power <- loess(power_estimates~n1)
plot(n1, predict(smooth_power, newdata = n1), type = "l",
     ylim = c(0,1), xlim = c(0, 500),
     ylab = "Power/Assurance",
     xlab = "Number of Patients (in each group)",
     cex.axis=1.5, cex.lab=1.5, cex.main=2)


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



