

#####Analsyis for this clinical trial: https://pubmed.ncbi.nlm.nih.gov/17968979/



##Chi square test for the data

MoxonidineData = data.frame(RaisedcTni = c(23, 31), NotRaisedcTni = c(40, 47), row.names = c("Control", "Treatment"))
chisq.test(moxonidine, correct = F)

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
output = coda.samples(model=model, variable.names = c("ass"), n.iter = 1000000)
mean(output[[1]])



#Same example but not using RJags
m = 10e5
zvec = vector(length=m)
for (i in 1:m){
  n1 <- 200 #Control
  n2 <- 400 #Treatment
  theta1 <- rbeta(1, 5, 20) #Control
  eff <- rbinom(1,1, 0.85)
  theff <- rbeta(1, 3, 4.5)
  thneff <- rbeta(1, 2, 23)
  theta2 <- (eff * theff)+((1-eff)*thneff) #Treatment
  control <- rbinom(1, n1, theta1)
  treatment <- rbinom(1, n2, theta2)
  control/n1
  p1 <-  control/n1
  p2 <-treatment/n2
  se <-  sqrt((p1*(1-p1)/n1)+(p2*(1-p2)/n2))
  zvec[i] <-  (p2-p1)/se
}

mean(zvec>1.96)



#Using the same method for our moxonidine trial

m = 10e3
zvec = vector(length=m)
for (i in 1:m){
  n1 <- 180 #Control
  n2 <- 180 #Treatment
  theta1 <- rbeta(1, 10.7, 13.1) #Control
  rho <- rnorm(1, 0.15, sd = sqrt(0.01))
  theta2 <- theta1 - rho #Treatment
  if (theta2<0){
    zvec[i] = NA
  } else{
    control <- rbinom(1, n1, theta1)
    treatment <- rbinom(1, n2, theta2)
    MoxonidineData = data.frame(RaisedcTni = c(control, treatment), NotRaisedcTni = c(n1 - control, n2 - treatment), row.names = c("Control", "Treatment"))
    test = chisq.test(MoxonidineData, correct = F)
    zvec[i] = test$p.value<0.05
  }
}

zvec = na.omit(zvec)
mean(zvec)


power = power.prop.test(n=300,p1=0.45,p2=0.3)
power$power

