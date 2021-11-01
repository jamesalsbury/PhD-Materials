

lambda1 = -log(0.6)/5
lambda2 = -log(0.8)/5
theta = log(lambda2/lambda1)
alpha = 0.05
N1 = 200
N2 = 200
T = 5
R = 3
P1e = 1-(exp(-lambda1*(T-R))-exp(-lambda1*T))/(lambda1*R)
P2e = 1-(exp(-lambda2*(T-R))-exp(-lambda2*T))/(lambda2*R)
N1vec = seq(1:500)
N2vec = seq(1:500)
power = vector(length=500)

powerfunc = function(N1, N2){
  pnorm(-theta/(sqrt(1/(N1*P1e)+1/(N2*P2e)))-qnorm(alpha/2), lower.tail = F)+pnorm(theta/(sqrt(1/(N1*P1e)+1/(N2*P2e)))-qnorm(alpha/2), lower.tail = F)
}

for (i in 1:500){
power[i] = powerfunc(N1vec[i], N2vec[i])  
}
plot(power, type="l", ylim=c(0,1), xaxt="n", xlab="Total sample size", ylab="Power")
axis(1, at=seq(0,500,by=100), labels = seq(0,1000, by=200))

#Assurance methods now
assurancefunc1 = function(N1, N2){
  M = 1000
  T = 5
  R = 3
  assurance = vector(length=M)
  for (j in 1:M){
    S1 = rbeta(1, 60, 40)
    rho = rnorm(1, 0.2, sqrt(0.001))
    S2 = S1+rho
    lambda1 = -log(S1)/5
    lambda2 = -log(S2)/5
    theta = log(lambda2/lambda1)
    P1e = 1-(exp(-lambda1*(T-R))-exp(-lambda1*T))/(lambda1*R)
    P2e = 1-(exp(-lambda2*(T-R))-exp(-lambda2*T))/(lambda2*R)
    assurance[j] = pnorm(theta/(sqrt(1/(N1*P1e)+1/(N2*P2e)))-qnorm(alpha/2), lower.tail = F)
  }
  assurance = na.omit(assurance)
  mean(assurance)  
}


ass1 = vector(length=500)
for (i in 1:500){
  ass1[i] = assurancefunc1(N1vec[i], N2vec[i])  
}
lines(ass1, lty=3, col="blue")

###

#Assurance methods now
assurancefunc2 = function(N1, N2){
  M = 1000
  T = 5
  R = 3
  assurance = vector(length=M)
  for (j in 1:M){
    S1 = rbeta(1, 60, 40)
    rho = rnorm(1, 0.2, sqrt(0.05))
    S2 = S1+rho
    lambda1 = -log(S1)/5
    lambda2 = -log(S2)/5
    theta = log(lambda2/lambda1)
    P1e = 1-(exp(-lambda1*(T-R))-exp(-lambda1*T))/(lambda1*R)
    P2e = 1-(exp(-lambda2*(T-R))-exp(-lambda2*T))/(lambda2*R)
    assurance[j] = pnorm(theta/(sqrt(1/(N1*P1e)+1/(N2*P2e)))-qnorm(alpha/2), lower.tail = F)
  }
  assurance = na.omit(assurance)
  mean(assurance)  
}


ass2 = vector(length=500)
for (i in 1:500){
  ass2[i] = assurancefunc2(N1vec[i], N2vec[i])  
}
lines(ass2, lty=4, col="purple")
ass2
###
#Assurance methods now
assurancefunc3 = function(N1, N2){
  M = 1000
  T = 5
  R = 3
  assurance = vector(length=M)
  for (j in 1:M){
    S1 = rbeta(1, 60, 40)
    rho = rnorm(1, 0.3, sqrt(0.005))
    S2 = S1+rho
    lambda1 = -log(S1)/5
    lambda2 = -log(S2)/5
    theta = log(lambda2/lambda1)
    P1e = 1-(exp(-lambda1*(T-R))-exp(-lambda1*T))/(lambda1*R)
    P2e = 1-(exp(-lambda2*(T-R))-exp(-lambda2*T))/(lambda2*R)
    assurance[j] = pnorm(theta/(sqrt(1/(N1*P1e)+1/(N2*P2e)))-qnorm(alpha/2), lower.tail = F)
  }
  assurance = na.omit(assurance)
  mean(assurance)  
}


ass3 = vector(length=500)
for (i in 1:500){
  ass3[i] = assurancefunc1(N1vec[i], N2vec[i])  
}
lines(ass3, lty=5, col="purple")








