

#####Analsyis for this clinical trial: https://pubmed.ncbi.nlm.nih.gov/17968979/


#Chi square test for the data

MoxonidineData <- data.frame(RaisedcTni = c(23, 31), NotRaisedcTni = c(40, 47), row.names = c("Control", "Treatment"))
chisq.test(MoxonidineData, correct = F)


#Eatimating beta parameters when we have an idea of mu and variance - appropriate for data [0,1]
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

#We have mean = 0.45 and var = 0.1
estBetaParams(0.45, 0.1^2)

#Plotting the prior beliefs for theta1
x = seq(0,1,by=0.01)
theta1 = dbeta(x, 10.7, 13)
plot(x,theta1, type="l", ylab="density", col="blue")


#Plotting the power curve for the moxonidine trial
N1vec <- seq(1,500)
N2vec <- seq(1,500)
powervec <- vector(length = 500)
for (i in 1:500){
  power = power.prop.test(n=N1vec[i],p1=0.45,p2=0.3)
  powervec[i] = power$power
  
}
plot(powervec, type = "l", ylim = c(0,1), xaxt = "n", xlab = "Total sample size", ylab = "Power/Assurance", lty=2)
axis(1, at = seq(0, 500,by = 125), labels = seq(0, 1000, by = 250))


#Assurance function for the moxonidine trial
assurancefunc <- function(n1, n2, m, v){
n = 10e1
zvec = vector(length=n)
for (i in 1:n){
  theta1 <- rbeta(1, 10.7, 13.1) #Control
  rho <- rnorm(1, m, sd = sqrt(v))
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
}

#Scenario 1
ass1 <- vector(length = 500)
for (i in 1:500){
  ass1[i] = assurancefunc(N1vec[i], N2vec[i], m = 0.15, v = 0.0001)
}
lo1 = loess(ass1~N1vec)
lines(predict(lo1), lty=3, col="blue")


#Scenario 2
ass2 <- vector(length = 500)
for (i in 1:500){
  ass2[i] = assurancefunc(N1vec[i], N2vec[i], m = 0.15, v = 0.01)
}

lo2a = loess(ass2~N1vec)
lines(predict(lo2a), lty=4, col="red")


#Scenario 3
ass3 <- vector(length = 500)
for (i in 1:500){
  ass3[i] = assurancefunc(N1vec[i], N2vec[i], m = 0.1, v  = 0.01)
}

lo3 = loess(ass3~N1vec)
lines(predict(lo3), lty=5, col="green")


legend("bottomright", legend = c("Power", "Scenario 1", "Scenario 2", "Scenario 3"), 
       col=c("black", "blue", "red", "green"), lty=2:5)



#Assurance at the interim analysis



MoxonidineData <- data.frame(cE = 23, tE = 31, cN = 63, tN = 78)


library(rjags)
data = list(cE = MoxonidineData$cE, cN = MoxonidineData$cN, 
            tE = MoxonidineData$tE, tN = MoxonidineData$tN)
modelstring="
model{
  cE~dbin(theta1,cN)
  theta1~dbeta(10.7, 13.1)
  tE~dbin(theta2, tN)
  rho~dnorm(0.15,1/(0.01))
  theta2 <- theta1 - rho
}
"
model = jags.model(textConnection(modelstring), data = data)
update(model, n.iter = 1000)
output = coda.samples(model=model, variable.names = c("theta1", "theta2"), n.iter = 1000)

interimassfunc = function(n1, n2){
  zvec = vector(length=1000)
  for (i in 1:1000){
    control <- rbinom(1, n1, output[[1]][i,1])
    treatment <- rbinom(1, n2, output[[1]][i,2])
    MoxonidineData <- data.frame(RaisedcTni = c(control, treatment),
                                 NotRaisedcTni = c(n1 - control, n2 - treatment), row.names = c("Control", "Treatment"))
    test <- chisq.test(MoxonidineData, correct = F)
    zvec[i] <- test$p.value < 0.05
  }
  mean(zvec)
}

#Scenario 1
N1vec <- seq(1,500)
N2vec <- seq(1,500)
interimass1 <- vector(length = 500)
for (i in 1:500){
  interimass1[i] = interimassfunc(N1vec[i], N2vec[i])
}
lo1 = loess(interimass1~N1vec)
lines(predict(lo1), lty=3, col="orange")

#Scenario 2
N1vec <- seq(1,500)
N2vec <- seq(1,500)
interimass2 <- vector(length = 500)
for (i in 1:500){
  interimass2[i] = interimassfunc(N1vec[i], N2vec[i])
}
lo2 = loess(interimass2~N1vec)
lines(predict(lo2), lty=3, col="blue")

