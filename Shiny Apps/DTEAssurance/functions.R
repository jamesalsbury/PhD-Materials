CensFunc <- function(dataCombined, censEvents = NULL, censTime = NULL){
  
  if (is.null(censEvents) && is.null(censTime)) {
    stop("Either censEvents or censTime must be specified")
  }
  
  if (!is.null(censEvents)){
    dataCombined <- dataCombined[order(dataCombined$pseudoTime),]
    censTime <- dataCombined$pseudoTime[censEvents]
  }
  
  
  
  dataCombined$status <- dataCombined$pseudoTime <= censTime
  dataCombined$status <- dataCombined$status * 1
  dataCombined$enrolled <- dataCombined$recTime < censTime
  dataCombined <- dataCombined[dataCombined$enrolled, ]
  dataCombined$survival_time <- ifelse(dataCombined$pseudoTime > censTime,
                                       censTime - dataCombined$recTime,
                                       dataCombined$time)
  
  return(list(dataCombined = dataCombined, censEvents = censEvents, censTime = censTime))
}


calculateAssurance <- function(n_C, n_E, lambda_C, HRStar, HRStarDist = "hist", gamma_C, gamma_E, delayT, delayTDist = "hist",
                               P_S = 1, P_DTE = 0, censEvents = NULL, censTime = NULL, rec_method, rec_period=NULL, rec_power=NULL, rec_rate=NULL, rec_duration=NULL,
                               analysis_method, rho = 0, gamma = 0, nSims=1e4){
  
  control_n <- length(lambda_C)
  
  
  assVec <- rep(NA, nSims)
  censVec <- rep(NA, nSims)
  
  for (i in 1:nSims){
    
    #Sample the control parameters
    u <- sample(1:control_n, size = 1)
    sampled_lambdac <- lambda_C[u]
    sampled_gammac <- gamma_C[u]
    
    #Sample the treatment effect parameters
    sampled_HRStar <- SHELF::sampleFit(HRStar, n = 1)[,HRStarDist]
    sampled_delayT <- SHELF::sampleFit(delayT, n = 1)[,delayTDist]
    
    #Make the simplifications
    sampled_gammae <- sampled_gammac
    
    
    dataCombined <- SimDTEDataSet(n_C, n_E, sampled_lambdac, sampled_HRStar, sampled_gammac, sampled_gammae, sampled_delayT, P_S, P_DTE,
                                  rec_method, rec_period, rec_power, rec_rate, rec_duration)
    
    
    censoredDF <- CensFunc(dataCombined, censEvents, censTime)
    
    dataCombined <- censoredDF$dataCombined
    
    censVec[i] <- censoredDF$censTime
    
    assVec[i] <- survivalAnalysis(dataCombined, analysis_method, alpha = 0.05, rho, gamma)
    
    
  }
  
  pHat <- mean(assVec)
  
  return(list(assurance = pHat, duration = mean(censVec), LBAssurance = pHat - 1.96*sqrt(pHat*(1-pHat)/nSims),
              UBAssurance = pHat + 1.96*sqrt(pHat*(1-pHat)/nSims)))
  
}


SimDTEDataSet <- function(n_C, n_E, lambda_C, HRStar, gamma_C, gamma_E, delayT, P_S = 1, P_DTE = 0,
                          rec_method, rec_period=NULL, rec_power=NULL, rec_rate=NULL, rec_duration=NULL) {
  
  #Simulates the control data
  u <- runif(n_C)
  controlTimes <- (1/lambda_C)*(-log(u))^(1/gamma_C)
  
  #Make the simplification
  gamma_E <- gamma_C
  
  #Takes into account the probabilities of the curves separating (and then being subject to a DTE)
  if (runif(1)>P_S){
    delayT <- 0
    HRStar <- 1
  } else if (runif(1)>P_DTE){
    delayT <- 0
  }
  
  #Transforming the post-delay HR to lambda_E
  lambda_E <- lambda_C*HRStar^(1/gamma_C)
  
  #Simulates the treatment data
  CP <- exp(-(lambda_C*delayT)^gamma_C)
  u <- runif(n_E)
  treatmentTimes <- ifelse(u>CP, (1/lambda_C)*(-log(u))^(1/gamma_C),
                           (1/(lambda_E^gamma_E)*((lambda_E*delayT)^gamma_E-log(u)-(lambda_C*delayT)^gamma_C))^(1/gamma_E))
  
  
  #Combines the control and treatment data
  dataCombined <- data.frame(time = c(controlTimes, treatmentTimes),
                             group = c(rep("Control", n_C), rep("Treatment", n_E)))
  
  n_total <- n_C + n_E
  
  if (rec_method=="power"){
    
    dataCombined$recTime <- rec_period * stats::runif(n_total)^(1/rec_power)
    
    dataCombined$pseudoTime <- dataCombined$time + dataCombined$recTime
  }
  
  if (rec_method=="PWC"){
    if(any(rec_rate<0)){stop("rec_rate should be non-negative")}
    if(length(rec_rate)==1){#simple case with only one rate
      rec<-cumsum(stats::rexp(n=n_total,rate=rec_rate))
    }else{#piecewise
      if(length(rec_duration)!=length(rec_rate)){stop("Lengths of rec_duration and rec_rate should match")}
      n_periods<-length(rec_duration)
      df<-data.frame(rate=rec_rate,
                     duration=rec_duration,
                     period=1:n_periods,
                     finish=cumsum(rec_duration),
                     lambda=rec_duration*rec_rate,
                     origin=c(0,cumsum(rec_duration)[-n_periods]))
      df$N<-sapply(df$lambda,function(x){stats::rpois(n=1,lambda = x)})
      if (sum(df$N)==0){
        if (df$rate[n_periods]==0) stop("Please specify positive rec_rate for the last period; otherwise enrollment cannot finish.")
        rec<-c(cumsum(stats::rexp(n_total,rate=df$rate[n_periods]))+df$finish[n_periods])
      }else{
        rec<-unlist(apply(df,1,function(x){sort(stats::runif(n=x[["N"]],min=x[["origin"]],max=x[["finish"]]))}))
        if (length(rec) >= n_total){rec<-rec[1:n_total]} # if n already achieved, return first n observations
        # stop with error message if enrollment has not finished but enrollment rate for last period is less or equal with 0
        else{if (df$rate[n_periods]==0){stop("Please specify positive rec_rate for the last period; otherwise enrollment cannot finish.")}
          # Otherwise, return inter-arrival exponential times
          rec<-c(rec, cumsum(stats::rexp(n_total-nrow(rec),rate=df$rate[n_periods]))+df$finish[n_periods])
        }
      }
      
      dataCombined$recTime <- rec
      
    }
    
    dataCombined$pseudoTime <- dataCombined$time + dataCombined$recTime
  }
  
  
  
  return(dataCombined)
}

normal.error <-
  function(parameters, values, probabilities, weights){
    sum(weights * (pnorm(values, parameters[1], exp(parameters[2])) - probabilities)^2)
  }

t.error <-
  function(parameters, values, probabilities, weights, degreesfreedom){
    sum(weights * (pt((values-parameters[1]) / exp(parameters[2]),
                      degreesfreedom) - probabilities)^2)
    
  }

gamma.error <-
  function(parameters, values, probabilities, weights){
    sum(weights * (pgamma(values, exp(parameters[1]), exp(parameters[2])) -probabilities)^2)
  }

lognormal.error <-
  function(parameters, values, probabilities, weights){
    sum(weights * (plnorm(values, parameters[1], exp(parameters[2])) - probabilities)^2)
  }

logt.error <-
  function(parameters, values, probabilities, weights, degreesfreedom){
    sum(weights * (pt((log(values) - parameters[1]) / exp(parameters[2]), degreesfreedom) - probabilities)^2)
  }
beta.error <-
  function(parameters, values, probabilities, weights){
    sum(weights * (pbeta(values, exp(parameters[1]), exp(parameters[2])) - probabilities)^2)
  }


fitdistdelayT <-
  function(vals, probs, lower = -Inf,
           upper = Inf, weights = 1, tdf = 3,
           expertnames = NULL,
           excludelogt = FALSE){
    
    if(is.matrix(vals)==F){vals<-matrix(vals, nrow = length(vals), ncol = 1)}
    if(is.matrix(probs)==F){probs <- matrix(probs, nrow = nrow(vals), ncol = ncol(vals))}
    if(is.matrix(weights)==F){weights <- matrix(weights, nrow = nrow(vals), ncol = ncol(vals))}
    if(length(lower)==1){lower <- rep(lower, ncol(vals))}
    if(length(upper)==1){upper <- rep(upper, ncol(vals))}
    if(length(tdf)==1){tdf <- rep(tdf, ncol(vals))}
    
    
    
    n.experts <- ncol(vals)
    normal.parameters <- matrix(NA, n.experts, 2)
    t.parameters <- matrix(NA, n.experts, 3)
    mirrorgamma.parameters <- gamma.parameters <-
      matrix(NA, n.experts, 2)
    mirrorlognormal.parameters <-
      lognormal.parameters <- matrix(NA, n.experts, 2)
    mirrorlogt.parameters <- logt.parameters <-
      matrix(NA, n.experts, 3)
    beta.parameters <- matrix(NA, n.experts, 2)
    ssq<-matrix(NA, n.experts, 9)
    
    colnames(ssq) <- c("normal", "t",
                       "gamma", "lognormal", "logt", "beta",
                       "mirrorgamma",
                       "mirrorlognormal",
                       "mirrorlogt")
    
    
    if(n.experts > 1 & n.experts < 27 & is.null(expertnames)){
      expertnames <- paste("expert.", LETTERS[1:n.experts], sep="")
    }
    
    if(n.experts > 27 & is.null(expertnames)){
      expertnames <- paste("expert.", 1:n.experts, sep="")
    }
    
    limits <- data.frame(lower = lower, upper = upper)
    row.names(limits) <- expertnames
    
    
    for(i in 1:n.experts){
      if (min(probs[,i]) > 0.4 ){stop("smallest elicited probability must be less than 0.4")}
      if (min(probs[,i]) < 0 | max(probs[,i]) > 1 ){stop("probabilities must be between 0 and 1")}
      if (max(probs[,i]) < 0.6 ){stop("largest elicited probability must be greater than 0.6")}
      if (min(vals[,i]) < lower[i]){stop("elicited parameter values cannot be smaller than lower parameter limit")}
      if (max(vals[,i]) > upper[i]){stop("elicited parameter values cannot be greater than upper parameter limit")}
      if (tdf[i] <= 0 ){stop("Student-t degrees of freedom must be greater than 0")}
      if (min(probs[-1,i] - probs[-nrow(probs),i]) < 0 ){stop("probabilities must be specified in ascending order")}
      if (min(vals[-1,i] - vals[-nrow(vals),i]) <= 0 ){stop("parameter values must be specified in ascending order")}
      
      
      # Need to exclude any probability judgements
      # P(X<=x) = 0 or P(X<=x) = 1
      # Should enforce these probabilities via the parameter limits
      
      inc <- (probs[, i] > 0) & (probs[, i] < 1)
      
      minprob <- min(probs[inc, i])
      maxprob <- max(probs[inc, i])
      minvals <- min(vals[inc, i])
      maxvals <- max(vals[inc, i])
      
      q.fit <- approx(x = probs[inc,i], y = vals[inc,i],
                      xout = c(0.4, 0.5, 0.6))$y
      l <- q.fit[1] # estimated 40th percentile on original scale
      u <- q.fit[3] # estimated 60th percentile on original scale
      
      # if(minprob > 0 & maxprob < 1){
      
      minq <- qnorm(minprob)
      maxq <- qnorm(maxprob)
      # Estimate m and v assuming X~N(m,v)
      
      # Obtain m by solving simultaneously:
      # m + Z_l \sqrt{v} = X_l
      # m + Z_u \sqrt{v} = X_u
      # where Z_a is a-th quantile from N(0, 1), X_a is a-th quantile of X
      m <- (minvals * maxq - maxvals * minq) / (maxq - minq)
      v <- ((maxvals - minvals) / (maxq - minq))^2
      # }else{
      #  minq <- qnorm(min(probs[probs[, i] > 0, i]))
      #  maxq <- qnorm(max(probs[probs[, i] < 1, i]))
      #  m <- q.fit[2] # Estimated median on original scale
      #  v<- (u - l)^2 / 0.25 # Estimated variance on original scale
      # }
      
      # Symmetric distribution fits ----
      
      normal.fit <- optim(c(m, 0.5*log(v)),
                          normal.error, values = vals[inc,i],
                          probabilities = probs[inc,i],
                          weights = weights[inc,i])
      normal.parameters[i,] <- c(normal.fit$par[1], exp(normal.fit$par[2]))
      ssq[i, "normal"] <- normal.fit$value
      
      # starting values: c(m, log((u - m)/ qt(0.6, tdf[i])))
      
      t.fit <- optim(c(m, 0.5*log(v)), t.error,
                     values = vals[inc,i],
                     probabilities = probs[inc,i],
                     weights = weights[inc,i],
                     degreesfreedom = tdf[i])
      t.parameters[i, 1:2] <- c(t.fit$par[1], exp(t.fit$par[2]))
      t.parameters[i, 3] <- tdf[i]
      ssq[i, "t"] <- t.fit$value
      
      # Positive skew distribution fits ----
      
      
      if(lower[i] > -Inf){
        vals.scaled1 <- vals[inc,i] - lower[i]
        m.scaled1 <- m - lower[i]
        
        gamma.fit<-optim(c(log(m.scaled1^2/v), log(m.scaled1/v)),
                         gamma.error, values = vals.scaled1,
                         probabilities = probs[inc,i],
                         weights = weights[inc,i])
        gamma.parameters[i,] <- exp(gamma.fit$par)
        ssq[i, "gamma"] <- gamma.fit$value
        
        std<-((log(u - lower[i])-log(l - lower[i]))/1.35)
        
        mlog <- (log(minvals - lower[i]) *
                   maxq - log(maxvals - lower[i]) * minq) /
          (maxq - minq)
        
        lognormal.fit <- optim(c(mlog,
                                 log(std)),
                               lognormal.error,
                               values = vals.scaled1,
                               probabilities = probs[inc,i],
                               weights = weights[inc,i])
        lognormal.parameters[i, 1:2] <- c(lognormal.fit$par[1],
                                          exp(lognormal.fit$par[2]))
        ssq[i, "lognormal"] <- lognormal.fit$value
        
        logt.fit <- optim(c(log(m.scaled1), log(std)),
                          logt.error,
                          values = vals.scaled1,
                          probabilities = probs[inc,i],
                          weights = weights[inc,i],
                          degreesfreedom = tdf[i])
        logt.parameters[i,1:2] <- c(logt.fit$par[1], exp(logt.fit$par[2]))
        logt.parameters[i,3] <- tdf[i]
        ssq[i, "logt"] <- logt.fit$value
      }
      
      # Beta distribution fits ----
      
      if((lower[i] > -Inf) & (upper[i] < Inf)){
        vals.scaled2 <- (vals[inc,i] - lower[i]) / (upper[i] - lower[i])
        m.scaled2 <- (m - lower[i]) / (upper[i] - lower[i])
        v.scaled2 <- v / (upper[i] - lower[i])^2
        
        alp <- abs(m.scaled2 ^3 / v.scaled2 * (1/m.scaled2-1) - m.scaled2)
        bet <- abs(alp/m.scaled2 - alp)
        if(identical(probs[inc, i],
                     (vals[inc, i] - lower[i]) / (upper[i] - lower[i]))){
          alp <- bet <- 1
        }
        beta.fit <- optim(c(log(alp), log(bet)),
                          beta.error,
                          values = vals.scaled2,
                          probabilities = probs[inc,i],
                          weights = weights[inc,i])
        beta.parameters[i,] <- exp(beta.fit$par)
        ssq[i, "beta"] <- beta.fit$value
        
      }
      
      # Negative skew distribution fits ----
      
      if(upper[i] < Inf){
        
        # Distributions are fitted to Y:= Upper limit - X
        
        valsMirrored <- upper[i] - vals[inc, i]
        probsMirrored <- 1 - probs[inc, i]
        mMirrored <- upper[i] - m
        
        # Mirror gamma
        
        
        
        mirrorgamma.fit<-optim(c(log(mMirrored^2/v), log(mMirrored/v)),
                               gamma.error, values = valsMirrored,
                               probabilities = probsMirrored,
                               weights = weights[inc,i])
        mirrorgamma.parameters[i,] <- exp(mirrorgamma.fit$par)
        ssq[i, "mirrorgamma"] <- mirrorgamma.fit$value
        
        # Mirror log normal
        
        
        # Obtain mlogMirror by solving simultaneously:
        # m + Z_l \sqrt{v} = Y_l
        # m + Z_u \sqrt{v} = Y_u
        # where Z_a is a-th quantile from N(0, 1),
        # Y_a is a-th quantile of Y
        # and we model Y = log(upper - X) ~ N(mlogMirror, stdMirror^2)
        
        
        mlogMirror <- (log(upper[i] - maxvals) *
                         (1 - minq) -
                         log(upper[i] - minvals) * (1-maxq)) /
          (maxq - minq)
        
        stdMirror <-((log(upper[i] - l)-log(upper[i] - u))/1.35)
        
        
        mirrorlognormal.fit <- optim(c(mlogMirror,
                                       log(stdMirror)),
                                     lognormal.error,
                                     values = valsMirrored,
                                     probabilities = probsMirrored,
                                     weights = weights[inc,i])
        mirrorlognormal.parameters[i, 1:2] <-
          c(mirrorlognormal.fit$par[1],
            exp(mirrorlognormal.fit$par[2]))
        ssq[i, "mirrorlognormal"] <- mirrorlognormal.fit$value
        
        # Mirror log t
        
        mirrorlogt.fit <- optim(c(log(mMirrored), log(stdMirror)),
                                logt.error,
                                values = valsMirrored,
                                probabilities = probsMirrored,
                                weights = weights[inc,i],
                                degreesfreedom = tdf[i])
        mirrorlogt.parameters[i,1:2] <- c(mirrorlogt.fit$par[1],
                                          exp(mirrorlogt.fit$par[2]))
        mirrorlogt.parameters[i,3] <- tdf[i]
        ssq[i, "mirrorlogt"] <- mirrorlogt.fit$value
        
      }
    }
    dfn <- data.frame(normal.parameters)
    names(dfn) <-c ("mean", "sd")
    row.names(dfn) <- expertnames
    
    dft <- data.frame(t.parameters)
    names(dft) <-c ("location", "scale", "df")
    row.names(dft) <- expertnames
    
    dfg <- data.frame(gamma.parameters)
    names(dfg) <-c ("shape", "rate")
    row.names(dfg) <- expertnames
    
    dfmirrorg <- data.frame(mirrorgamma.parameters)
    names(dfmirrorg) <-c ("shape", "rate")
    row.names(dfmirrorg) <- expertnames
    
    dfln <- data.frame(lognormal.parameters)
    names(dfln) <-c ("mean.log.X", "sd.log.X")
    row.names(dfln) <- expertnames
    
    dfmirrorln <- data.frame(mirrorlognormal.parameters)
    names(dfmirrorln) <-c ("mean.log.X", "sd.log.X")
    row.names(dfmirrorln) <- expertnames
    
    dflt <- data.frame(logt.parameters)
    names(dflt) <-c ("location.log.X", "scale.log.X", "df.log.X")
    row.names(dflt) <- expertnames
    
    dfmirrorlt <- data.frame(mirrorlogt.parameters)
    names(dfmirrorlt) <-c ("location.log.X", "scale.log.X", "df.log.X")
    row.names(dfmirrorlt) <- expertnames
    
    dfb <- data.frame(beta.parameters)
    names(dfb) <-c ("shape1", "shape2")
    row.names(dfb) <- expertnames
    
    ssq <- data.frame(ssq)
    row.names(ssq) <- expertnames
    
    ssqT <- ssq[,c("gamma", "lognormal", "beta")]
    
    if(excludelogt){
      reducedssq <- ssq[, c("normal", "t", "gamma",
                            "lognormal", "beta",
                            "mirrorgamma",
                            "mirrorlognormal")]
      index <- apply(reducedssq, 1, which.min)
      best.fitting <- data.frame(best.fit=
                                   names(reducedssq)[index])}
    else{
      index <- apply(ssqT, 1, which.min)
      
      best.fitting <- data.frame(best.fit=names(ssqT)[index])
    }
    
    
    
    row.names(best.fitting) <- expertnames
    
    vals <- data.frame(vals)
    names(vals) <- expertnames
    
    probs <- data.frame(probs)
    names(probs) <- expertnames
    
    fit <- list(Normal = dfn, Student.t = dft,
                Gamma = dfg, Log.normal = dfln,
                Log.Student.t = dflt, Beta = dfb,
                mirrorgamma = dfmirrorg,
                mirrorlognormal = dfmirrorln,
                mirrorlogt = dfmirrorlt,
                ssq = ssq,
                best.fitting = best.fitting, vals = t(vals),
                probs = t(probs), limits = limits)
    class(fit) <- "elicitation"
    fit
  }










