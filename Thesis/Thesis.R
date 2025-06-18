library(pwr)
png("PowerPlot.png", units="in", width=8, height=5, res=700)
outcome <- sapply(1:200, function(x) pwr.t.test(n = x, d = 0.5)$power)
plot(1:200, outcome, type = "l", xlab = "Number of Patients (in each group)", 
     ylab = "Power", ylim = c(0,1), col = "blue")
legend("topleft", legend = "Power", col = "blue", lty = 1)
dev.off()