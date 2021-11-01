

#####Analsyis for this clinical trial: https://pubmed.ncbi.nlm.nih.gov/17968979/

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

estBetaParams(0.45, 0.1^2)

x = seq(0,1,by=0.01)
theta1 = dbeta(x, 10.7, 13)
plot(x,theta1, type="l", ylab="density", col="blue")

estBetaParams(0.5985, 0.2^2)
theta2 = dbeta(x, 3, 2)

lines(x, theta2, type="l", lty=2, col="red")
legend(0, 3.5, legend=c("theta1", "theta2"), lty=1:2, col=c("blue", "red"))

library(rjags)
data = list(n1=180, n2=180, zsig=1.96)
modelstring="
model{
  theta1~dbeta(10.7, 13)
  theta2~dbeta(3, 2)
  r1~dbin(theta1, n1)
  r2~dbin(theta2, n2)
  p1 <-  r1/n1
  p2 <- r2/n2
  se <-  sqrt((p1*(1-p1)/n1)+(p2*(1-p2)/n2))
  z <-  (p2-p1)/se
  ass <- step(z-zsig) # The mean of ass is the assurance
}
"
model = jags.model(textConnection(modelstring), data = data)
update(model, n.iter = 1000)
output = coda.samples(model=model, variable.names = c("ass"), n.iter = 1000000)
mean(output[[1]])




####Analysis for the binary data section in O'Hagan (2005) - Example 4
library(rjags)
data = list(n1=200, n2=400, zsig=1.96)
modelstring="
model{
  theta1~dbeta(5,20)
  theta2 <- (eff * theff)+((1-eff)*thneff) # These lines set out
  eff~dbern(0.85) # the prior as in Example 4.
  theff~dbeta(3,4.5) # They should be modified
  thneff~dbeta(2, 23) # as appropriate for any real example.
  r1~dbin(theta1, n1)
  r2~dbin(theta2, n2)
  p1 <-  r1/n1
  p2 <- r2/n2
  se <-  sqrt((p1*(1-p1)/n1)+(p2*(1-p2)/n2))
  z <-  (p2-p1)/se
  ass <- step(z-zsig) # The mean of ass is the assurance
}
"
model = jags.model(textConnection(modelstring), data = data)
update(model, n.iter = 1000)
output = coda.samples(model=model, variable.names = c("ass"), n.iter = 10000000)
mean(output[[1]])

