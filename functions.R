

simulateDTEWeibullData <- function(bigT, lambda2, gamma2, lambda1, gamma1, n1, n2){
  controldata <-  rweibull(n1, gamma2, 1/lambda2)
  
  CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
  u <- runif(n2)
  
  suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
  
  return(list(controldata = controldata, treatmentdata = z))
}



