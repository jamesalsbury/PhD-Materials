

#####Analsyis for this clinical trial: https://pubmed.ncbi.nlm.nih.gov/17968979/


#Chi square test for the data

MoxonidineData <- data.frame(RaisedcTni = c(23, 31), NotRaisedcTni = c(40, 47), row.names = c("Control", "Treatment"))
chisq.test(MoxonidineData, correct = F)


#Eatimating beta parameters when we have an idea of mu and variance - appropriate for data [0,1]
estBetaParams <- function(mu, var) {
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

#We have mean = 0.45 and var = 0.1
estBetaParams(0.45, 0.1^2)

#Plotting the prior beliefs for theta1
x = seq(0,1,by=0.01)
theta1 = dbeta(x, 10.7, 13)
plot(x,theta1, type="l", ylab="density", col="blue")


#Plotting the power curve for the moxonidine trial - at the design stage

powerfunc <- function(){
  N1vec <- seq(1,500)
  N2vec <- seq(1,500)
  powervec <- vector(length = 500)
  for (i in 1:500){
    power = power.prop.test(n=N1vec[i],p1=0.45,p2=0.3)
    powervec[i] = power$power
  }
  return(powervec)
}

#Assurance function for the moxonidine trial - at the design stage

#We input mu and nu (the parameters for the rho distribution)
DesignAssuranceFunc <- function(mu, nu){
  
  #We calculate the assurance for user-specified n1 and n2
  assurancefunc <- function(n1, n2, mu, nu){
    n = 10e1
    zvec = vector(length=n)
    for (i in 1:n){
      theta1 <- rbeta(1, 10.7, 13.1) #Control
      rho <- rnorm(1, mu, sd = sqrt(nu))
      theta2 <- theta1 - rho #Treatment
      if (theta2<0){
        zvec[i] <- NA
      } else{
        control <- rbinom(1, n1, theta1)
        treatment <- rbinom(1, n2, theta2)
        MoxonidineData <- data.frame(RaisedcTni = c(control, treatment), NotRaisedcTni = c(n1 - control, n2 - treatment), row.names = c("Control", "Treatment"))
        test <- chisq.test(MoxonidineData, correct = F)
        zvec[i] <- test$p.value<0.05
      }
    }
    zvec <- na.omit(zvec)
    mean(zvec)
  }
  
  #We compute the assurance for a range of different n1, n2 values
  assurance <- vector(length = 500)
  N1vec = 1:500
  N2vec = 1:500
  for (i in 1:500){
    assurance[i] <- assurancefunc(N1vec[i], N2vec[i], mu = mu, nu = nu)
  }
  
  #Smooth the output
  lo <- loess(assurance~N1vec)
  return(lo)
}


#Assurance at the interim analysis


InterimAssuranceFunc <- function(mu, nu){

  #Find the posterior dist. for theta1 and theta2
  library(rjags)
  MoxonidineData <- data.frame(cE = 23, tE = 31, cN = 63, tN = 78)
  data = list(cE = MoxonidineData$cE, cN = MoxonidineData$cN, 
              tE = MoxonidineData$tE, tN = MoxonidineData$tN, mu=mu, nu=nu)
  modelstring="
  model{
    cE~dbin(theta1,cN)
    theta1~dbeta(10.7, 13.1)
    tE~dbin(theta2, tN)
    rho~dnorm(mu,1/(nu))
    theta2 <- theta1 - rho
}
"
model <- jags.model(textConnection(modelstring), data = data)
update(model, n.iter = 1000)
output <- coda.samples(model=model, variable.names = c("theta1", "theta2"), n.iter = 1000)


#For a given n1 and n2, compute the assurance
interimassfunc <- function(n1, n2, output){
  zvec <- vector(length=1000)
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

#Compute the assurance for a range of n1 and n2 values
N1vec <- seq(6,501, by=5)
N2vec <- seq(6,501, by=5)
interimass <- vector(length = 100)
for (i in 1:100){
  interimass[i] <- interimassfunc(N1vec[i], N2vec[i], output = output)
}
lo <- loess(interimass~N1vec)
return(lo)
}


#Plotting power, assurance before the trial and assurance at the interim analysis for different scenarios
MoxonidineTrialPowerAndAssurance <- function(mu, nu, scenario){
  plot(powerfunc(), type = "l", ylim = c(0,1), xaxt = "n", xlab = "Total sample size", ylab = "Power/Assurance", lty=2,
       main = bquote(phantom(0)["Scenario "*.(scenario)* ":" ~ rho  ~ "~ N(" *.(mu)*","*.(nu)*")"~.]))
  axis(1, at = seq(0, 500,by = 125), labels = seq(0, 1000, by = 250))
  lines(predict(DesignAssuranceFunc(mu, nu)), col="blue", lty=3)
  N1vec <- seq(6,496, by=5)
  lines(N1vec, predict(InterimAssuranceFunc(mu, nu)), col='red', lty=4)
  legend("bottomright", legend = c("Power", "Design assurance", "Interim assurance"), 
         col=c("black", "blue", "red"), lty=2:4)
}


#Scenario 1
MoxonidineTrialPowerAndAssurance(mu = 0.15, nu = 0.0001, scenario = 1)

#Scenario 2
MoxonidineTrialPowerAndAssurance(mu = 0.15, nu = 0.01, scenario = 2)

#Scenario 3
MoxonidineTrialPowerAndAssurance(mu = 0.10, nu = 0.01, scenario = 3)




#Try and do it using STAN instead
#This is what we are trying to emulate

#Find the posterior dist. for theta1 and theta2
library(rjags)
mu = 0.15
nu = 0.01
MoxonidineData <- data.frame(cE = 23, tE = 31, cN = 63, tN = 78)
data = list(cE = MoxonidineData$cE, cN = MoxonidineData$cN, 
            tE = MoxonidineData$tE, tN = MoxonidineData$tN, mu=mu, nu=nu)
modelstring="
  model{
    cE~dbin(theta1,cN)
    theta1~dbeta(10.7, 13.1)
    tE~dbin(theta2, tN)
    rho~dnorm(mu,1/(nu))
    theta2 <- theta1 - rho
}
"
model <- jags.model(textConnection(modelstring), data = data, quiet=T)
update(model, n.iter = 1000)
output <- coda.samples(model=model, variable.names = c("theta1", "theta2"), n.iter = 1000)



#For a given n1 and n2, compute the assurance
interimassfunc <- function(n1, n2, output){
  zvec <- vector(length=1000)
  for (i in 1:1000){
    control <- rbinom(1, n1, output[[1]][i,1])
    treatment <- rbinom(1, n2, output[[1]][i,2])
    MoxonidineData <- data.frame(RaisedcTni = c(control, treatment),
                                 NotRaisedcTni = c(n1 - control, n2 - treatment), 
                                 row.names = c("Control", "Treatment"))
    test <- chisq.test(MoxonidineData, correct = F)
    zvec[i] <- test$p.value < 0.05
  }
  mean(zvec)
}

interimassfunc(180, 180, output)


#STAN code
library(rstan)

bern.stan =
  "
data {
  int<lower=0> cE;               // control events
  int<lower=0> cN;               // number of patients in control group
  int<lower=0> tE;               // treatment events
  int<lower=0> tN;               // number of patients in treatment group
  real mu;                       // mean of rho
  real sigma;                    // sd of rho
}
parameters {
  real<lower=0, upper=1> theta1; // chance of success in control group
  real rho;                      // difference between treatment and control (theta1-theta2)
}
transformed parameters{
  real<lower=0, upper=1> theta2; // chance of success in treatment group
  theta2 = theta1 - rho;
}
model {
  theta1 ~ beta(10.7, 13.1);       
  rho ~ normal(mu, sigma);
  cE ~ binomial(cN, theta1);        
  tE ~ binomial(tN, theta2);
}
"
MoxonidineData <- data.frame(cE = 23, tE = 31, cN = 63, tN = 78)
data = list(cE = MoxonidineData$cE, cN = MoxonidineData$cN, 
            tE = MoxonidineData$tE, tN = MoxonidineData$tN, mu = 0.15, sigma = 0.1)

fit = stan(model_code=bern.stan, data=data)


interimassfunc <- function(n1, n2, fit){
  zvec <- vector(length=1000)
  for (i in 1:1000){
    control <- rbinom(1, n1, fit@sim[["samples"]][[1]][[1]][i])
    treatment <- rbinom(1, n2, fit@sim[["samples"]][[1]][[3]][i])
    MoxonidineData <- data.frame(RaisedcTni = c(control, treatment),
                                 NotRaisedcTni = c(n1 - control, n2 - treatment), 
                                 row.names = c("Control", "Treatment"))
    test <- chisq.test(MoxonidineData, correct = F)
    zvec[i] <- test$p.value < 0.05
  }
  mean(zvec)
}

N1vec <- seq(6,496, by=5)
N2vec <- seq(6,501, by=5)
interimass <- vector(length = 100)
for (i in 1:100){
  interimass[i] <- interimassfunc(N1vec[i], N2vec[i], fit = fit)
}
lo <- loess(interimass~N1vec)
lines(N1vec, predict(lo), col="green")
return(lo)
interimass

length(predict(lo))






InterimAssuranceFunc <- function(mu, nu){
  mu = 0.15
  nu= 0.01
  #Find the posterior dist. for theta1 and theta2
  library(rstan)
  
  bern.stan =
    "
data {
  int<lower=0> cE;               // control events
  int<lower=0> cN;               // number of patients in control group
  int<lower=0> tE;               // treatment events
  int<lower=0> tN;               // number of patients in treatment group
  real mu;                       // mean of rho
  real sigma;                    // sd of rho
}
parameters {
  real<lower=0, upper=1> theta1; // chance of success in control group
  real rho;                      // difference between treatment and control (theta1-theta2)
}
transformed parameters{
  real<lower=0, upper=1> theta2; // chance of success in treatment group
  theta2 = theta1 - rho;
}
model {
  theta1 ~ beta(10.7, 13.1);       
  rho ~ normal(mu, sigma);
  cE ~ binomial(cN, theta1);        
  tE ~ binomial(tN, theta2);
}
"
MoxonidineData <- data.frame(cE = 23, tE = 31, cN = 63, tN = 78)
data <- list(cE = MoxonidineData$cE, cN = MoxonidineData$cN, 
             tE = MoxonidineData$tE, tN = MoxonidineData$tN, mu = mu, sigma = sqrt(nu))

fit <- stan(model_code=bern.stan, data=data)


#For a given n1 and n2, compute the assurance
interimassfunc <- function(n1, n2, fit){
  zvec <- vector(length=1000)
  for (i in 1:1000){
    control <- rbinom(1, n1, fit@sim[["samples"]][[1]][[1]][i])
    treatment <- rbinom(1, n2, fit@sim[["samples"]][[1]][[3]][i])
    MoxonidineData <- data.frame(RaisedcTni = c(control, treatment),
                                 NotRaisedcTni = c(n1 - control, n2 - treatment), 
                                 row.names = c("Control", "Treatment"))
    test <- chisq.test(MoxonidineData, correct = F)
    zvec[i] <- test$p.value < 0.05
  }
  mean(zvec)
}

#Compute the assurance for a range of n1 and n2 values
N1vec <- seq(10,180, by=2)
N2vec <- seq(10,180, by=2)
interimass <- vector(length = 86)
for (i in 1:86){
  interimass[i] <- interimassfunc(N1vec[i], N2vec[i], fit = fit)
}
lo <- loess(interimass~N1vec)
return(lo)
}




#######################R markdown
powerfunc <- function(){
  N1vec <- seq(1,180)
  N2vec <- seq(1,180)
  powervec <- vector(length = 180)
  for (i in 1:180){
    power <- power.prop.test(n=N1vec[i],p1=0.45,p2=0.3)
    powervec[i] <- power$power
  }
  return(powervec)
}

#Assurance function for the moxonidine trial - at the design stage

#We input mu and nu (the parameters for the rho distribution)
DesignAssuranceFunc <- function(mu, nu){
  
  #We calculate the assurance for user-specified n1 and n2
  assurancefunc <- function(n1, n2, mu, nu){
    n <- 10e1
    zvec <- vector(length=n)
    for (i in 1:n){
      theta1 <- rbeta(1, 10.7, 13.1) #Control
      rho <- rnorm(1, mu, sd = sqrt(nu))
      theta2 <- theta1 - rho #Treatment
      if (theta2<0){
        zvec[i] <- NA
      } else{
        control <- rbinom(1, n1, theta1)
        treatment <- rbinom(1, n2, theta2)
        MoxonidineData <- data.frame(RaisedcTni = c(control, treatment),
                                     NotRaisedcTni = c(n1 - control, n2 - treatment), 
                                     row.names = c("Control", "Treatment"))
        test <- chisq.test(MoxonidineData, correct = F)
        zvec[i] <- test$p.value<0.05
      }
    }
    zvec <- na.omit(zvec)
    mean(zvec)
  }
  
  #We compute the assurance for a range of different n1, n2 values
  assurance <- vector(length = 180)
  N1vec <- 1:180
  N2vec <- 1:180
  for (i in 1:180){
    assurance[i] <- assurancefunc(N1vec[i], N2vec[i], mu = mu, nu = nu)
  }
  
  #Smooth the output
  lo <- loess(assurance~N1vec)
  return(lo)
}


#Assurance at the interim analysis


InterimAssuranceFunc <- function(mu, nu){
  
  #Find the posterior dist. for theta1 and theta2
  library(rstan)
  
  bern.stan =
    "
data {
  int<lower=0> cE;               // control events
  int<lower=0> cN;               // number of patients in control group
  int<lower=0> tE;               // treatment events
  int<lower=0> tN;               // number of patients in treatment group
  real mu;                       // mean of rho
  real sigma;                    // sd of rho
}
parameters {
  real<lower=0, upper=1> theta1; // chance of success in control group
  real rho;                      // difference between treatment and control (theta1-theta2)
}
transformed parameters{
  real<lower=0, upper=1> theta2; // chance of success in treatment group
  theta2 = theta1 - rho;
}
model {
  theta1 ~ beta(10.7, 13.1);       
  rho ~ normal(mu, sigma);
  cE ~ binomial(cN, theta1);        
  tE ~ binomial(tN, theta2);
}
"
MoxonidineData <- data.frame(cE = 23, tE = 31, cN = 63, tN = 78)
data <- list(cE = MoxonidineData$cE, cN = MoxonidineData$cN, 
             tE = MoxonidineData$tE, tN = MoxonidineData$tN, mu = mu, sigma = sqrt(nu))

fit <- stan(model_code=bern.stan, data=data)


#For a given n1 and n2, compute the assurance
interimassfunc <- function(n1, n2, fit){
  zvec <- vector(length=1000)
  for (i in 1:1000){
    control <- rbinom(1, n1, fit@sim[["samples"]][[1]][[1]][i])
    treatment <- rbinom(1, n2, fit@sim[["samples"]][[1]][[3]][i])
    MoxonidineData <- data.frame(RaisedcTni = c(control, treatment),
                                 NotRaisedcTni = c(n1 - control, n2 - treatment), 
                                 row.names = c("Control", "Treatment"))
    test <- chisq.test(MoxonidineData, correct = F)
    zvec[i] <- test$p.value < 0.05
  }
  mean(zvec)
}

#Compute the assurance for a range of n1 and n2 values
N1vec <- seq(10,180, by=5)
N2vec <- seq(10,180, by=5)
interimass <- vector(length = 35)
for (i in 1:35){
  interimass[i] <- interimassfunc(N1vec[i], N2vec[i], fit = fit)
}
lo <- loess(interimass~N1vec)
return(lo)
}


#Plotting power, assurance before trial and interim assurance for different scenarios
MoxonidineTrialPowerAndAssurance <- function(mu, nu, scenario){
  plot(powerfunc(), type = "l", ylim = c(0,1), xaxt = "n", xlab = "Total sample size",
       ylab = "Power/Assurance", lty=2,
       main = bquote(phantom(0)["Scenario "*.(scenario)* ":" ~ rho  ~ "~ N(" *.(mu)*","*.(nu)*")"~.]))
  axis(1, at = seq(0, 200,by = 50), labels = seq(0, 400, by = 100))
  lines(predict(DesignAssuranceFunc(mu, nu)), col="blue", lty=3)
  N1vec <- seq(10,180, by=5)
  lines(N1vec, predict(InterimAssuranceFunc(mu, nu)), col='red', lty=4)
  legend("bottomright", legend = c("Power", "Design assurance", "Interim assurance"), 
         col=c("black", "blue", "red"), lty=2:4)
}


#Scenario 1
MoxonidineTrialPowerAndAssurance(mu = 0.15, nu = 0.0001, scenario = 1)

#Scenario 2
MoxonidineTrialPowerAndAssurance(mu = 0.15, nu = 0.01, scenario = 2)

#Scenario 3
MoxonidineTrialPowerAndAssurance(mu = 0.10, nu = 0.01, scenario = 3)



############Non-informative prior

#We input mu and nu (the parameters for the rho distribution)
DesignAssuranceFunc <- function(){
  
  #We calculate the assurance for user-specified n1 and n2
  assurancefunc <- function(n1, n2){
    n <- 10e1
    zvec <- vector(length=n)
    for (i in 1:n){
      theta1 <- rbeta(1, 10.7, 13.1) #Control
      rho <- runif(1)
      theta2 <- theta1 - rho #Treatment
      if (theta2<0){
        zvec[i] <- NA
      } else{
        control <- rbinom(1, n1, theta1)
        treatment <- rbinom(1, n2, theta2)
        MoxonidineData <- data.frame(RaisedcTni = c(control, treatment),
                                     NotRaisedcTni = c(n1 - control, n2 - treatment), 
                                     row.names = c("Control", "Treatment"))
        test <- chisq.test(MoxonidineData, correct = F)
        zvec[i] <- test$p.value<0.05
      }
    }
    zvec <- na.omit(zvec)
    mean(zvec)
  }
  
  #We compute the assurance for a range of different n1, n2 values
  assurance <- vector(length = 180)
  N1vec <- 1:180
  N2vec <- 1:180
  for (i in 1:180){
    assurance[i] <- assurancefunc(N1vec[i], N2vec[i])
  }
  
  #Smooth the output
  lo <- loess(assurance~N1vec)
  return(lo)
}
x = DesignAssuranceFunc()
plot(predict(x))
