
library(nleqslv)
library(shiny)
library(survival)
library(shinydashboard)
library(rsconnect)


ui <- fluidPage(
  
  titlePanel("Cancer survival times - Weibull parameterisation"),
  
  sidebarLayout(
    
    sidebarPanel(
      fluidRow(
        box(width = 14, title = "When does the treatment begin to take effect (months)?", 
            splitLayout(
              numericInput("T1Mean", "Mean", value=5),
              numericInput("T1Var", "Variance", value=1, min=0)
            )
        )
      ),
      fluidRow(
        box(width = 14, title = "What is the hazard ratio after the change-point?", 
            splitLayout(
              numericInput("T2HRMean", "Mean", value=0.65),
              numericInput("T2HRVar", "Variance", value=0.00001, min=0)
            )
        )
      ),
      checkboxGroupInput("showfeedback", "Add to plot", choices = c("Median survival line",
                                                                    "Hazard Ratio & 95% CI's", "95% CI for T", "Simulation curves")),
      fluidRow(
        box(width = 14, title = "Ratio of patients in each group?",
            splitLayout(
              numericInput("n1", "Control", value=1, min=1),
              numericInput("n2", "Treatment", value=2, min=1)
            )
        )
      ),
    ),
    
    mainPanel(
      
      tabsetPanel(type="tabs",
                  tabPanel("Elicitation",
                           plotOutput("plotBestFit"),
                           htmlOutput("feedback"), 
                  ), 
                  tabPanel("Assurance", 
                           plotOutput("plotAssurance")))
    )
  )
)


server <- function(input, output, session) {
  
  
  
  #lambda2 <- 0.06
  #gamma2 <- 0.8
  
  #Simulate data for the control, in pratice the data would just be given to us
  simdata <<- data.frame(time = rweibull(10000, 0.8, 1/0.06), cens = rep(1, 10000))
  
  #Determine lambda2 and gamma2 from the control data
  fitcontrol <<- survreg(Surv(time, cens)~1, dist="weibull", data = simdata)
  lambda2 <- exp(-fitcontrol$coefficients)
  gamma2 <- 1/fitcontrol$scale
  
  #Let gamma1 = gamma2 here
  gamma1 <<- gamma2
  
  #Calculate lambda1 from the inputs
  findlambda1 <<- reactive({
    lambda1 <- exp((log(input$T2HRMean)+gamma1*log(lambda2))/gamma1)
    return(lambda1)
  })
  
  estBetaParams <<- function(mu, var) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
  }
  
  proposed <- reactive({
    lambda1 <- findlambda1()
    
    #Plot the proposed treatment curve (not interested in variation here)
    #Plot the treatment curve until the change-point (just same as control) 
    effectt <- seq(0, input$T1Mean, by=0.01)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    
    
    #Plot the treatment curve after the change-point 
    aftereffectt <- seq(input$T1Mean, max(simdata$time)*1.1, by=0.01)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*input$T1Mean)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-input$T1Mean^gamma1))
    
    combinedt <- c(effectt, aftereffectt)
    combinedy <- c(effecty, aftereffecty)
    
    list(effectt = effectt, effecty = effecty, aftereffectt = aftereffectt, 
         aftereffecty = aftereffecty, combinedt = combinedt, combinedy = combinedy)
  })
  
  fivecurve <- reactive({
    
    #Calculating lambda2 and gamma2
    lambda2 <- exp(-fitcontrol$coefficients)
    gamma2 <- gamma1
    
    #Finding the 5% and 95% estimates for T and HR (given the inputs from the user)
    bigT <- qnorm(0.025, mean = input$T1Mean, sd = sqrt(input$T1Var))
    HRdist <- estBetaParams(mu = input$T2HRMean , var = sqrt(input$T2HRVar))
    HR <- qbeta(0.025, HRdist$alpha, HRdist$beta)
    
    
    #Calculating lambda1 in the case of HR being at 5%
    lambda1 <- exp((log(HR)+gamma1*log(lambda2))/gamma1)
    
    #Plotting the 5% lines (when T and HR is at 5%)
    #Before CP
    effectt <- seq(0, bigT, by=0.01)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    
    
    #After CP
    aftereffectt <- seq(bigT, max(simdata$time)*1.1, by=0.01)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*bigT)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-bigT^gamma1))
    
    combinedt <- c(effectt, aftereffectt)
    combinedy <- c(effecty, aftereffecty)
    
    list(effectt = effectt, effecty = effecty, aftereffectt = aftereffectt, 
         aftereffecty = aftereffecty, combinedt = combinedt, combinedy = combinedy,
         bigT = bigT, HR = HR)
  })
  
  ninetyfivecurve <- reactive({
    
    #Calculating lambda2 and gamma2
    lambda2 <- exp(-fitcontrol$coefficients)
    gamma2 <- gamma1
    
    #Finding the 5% and 95% estimates for T and HR (given the inputs from the user)
    bigT <- qnorm(0.975, mean = input$T1Mean, sd = sqrt(input$T1Var))
    HRdist <- estBetaParams(mu = input$T2HRMean , var = sqrt(input$T2HRVar))
    HR <- qbeta(0.975, HRdist$alpha, HRdist$beta)
    
    
    #Calculating lambda1 in the case of HR being at 5%
    lambda1 <- exp((log(HR)+gamma1*log(lambda2))/gamma1)
    
    #Plotting the 5% lines (when T and HR is at 5%)
    #Before CP
    effectt <- seq(0, bigT, by=0.01)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    
    
    #After CP
    aftereffectt <- seq(bigT, max(simdata$time)*1.1, by=0.01)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*bigT)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-bigT^gamma1))
    
    
    combinedt <- c(effectt, aftereffectt)
    combinedy <- c(effecty, aftereffecty)
    
    list(effectt = effectt, effecty = effecty, aftereffectt = aftereffectt, 
         aftereffecty = aftereffecty, combinedt = combinedt, combinedy = combinedy,
         bigT = bigT, HR = HR)
  })
  
  
  drawsimlines <- reactive({
    lambda2 <- exp(-fitcontrol$coefficients)
    gamma2 <- gamma1
    linelist <- list()
    for (i in 1:10){
      bigT <- rnorm(1, mean = input$T1Mean, sd = sqrt(input$T1Var))
      HRdist <- estBetaParams(mu = input$T2HRMean , var = sqrt(input$T2HRVar))
      HR <- rbeta(1, HRdist$alpha, HRdist$beta)
      lambda1 <- exp((log(HR)+gamma1*log(lambda2))/gamma1)
      aftereffectt <- seq(bigT, max(simdata$time)*1.1, by=0.01)
      aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*bigT)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-bigT^gamma1))
      linelist[[(i*2)-1]] <- aftereffectt
      linelist[[i*2]] <- aftereffecty
    }
    list(linelist = linelist)
  })
  
  
  output$plotBestFit <- renderPlot({
    
    #Plotting the control data using the Weibull parameters found
    plot(x = predict(fitcontrol, type = "quantile", p = seq(0.001, 0.999, by=.001))[1,],
         y = rev(seq(0.001, 0.999, by = 0.001)), type="l", xlab="Time (months)", ylab="Survival",
         col = "blue", xlim=c(0,sort(simdata$time)[length(simdata$time)*0.99]), main="Elicitation of treatment curve")
    
    
    lines(proposed()$effectt, proposed()$effecty, col="green", lty=2)
    lines(proposed()$aftereffectt, proposed()$aftereffecty, col="red", lwd=1.5)
    #abline(h = 0.5)
    
    legend("topright", legend = c("Weibull fit to control data", "Control + Treatment both Weibull",
                                  "Proposed treatment survival curve"),
           col=c("blue", "green", "red"), lty=c(1, 3, 1), cex=0.75)
    
    
    addfeedback <- input$showfeedback 
    
    if (is.null(addfeedback)){
      CalcSim <<- TRUE
    }
    
    if (!is.null(addfeedback)){
      CheckSim <<- FALSE
      for (i in 1:length(addfeedback)){
        if (addfeedback[i]=="Simulation curves"){
          CheckSim <- TRUE
        }
      }
      if (CheckSim == FALSE){
        CalcSim <<- TRUE
      }
    }
    
    if (!is.null(addfeedback)){
      for (i in 1:length(addfeedback)){
        if (addfeedback[i]=="Median survival line"){
          lines(seq(0, 120, length=2), rep(0.5, 2), lty=3)
        } else if (addfeedback[i]=="Hazard Ratio & 95% CI's"){
          #Top line
          lines(seq(0, input$T1Mean, length=2), rep(1, 2))
          #Vertical line
          lines(rep(input$T1Mean, 2), seq(input$T2HRMean, 1, length=2))
          #Bottom line
          lines(seq(input$T1Mean, 120,length=2), rep(input$T2HRMean, 2))
          #Conf intervals
          lines(seq(input$T1Mean, 120,length=2), rep(fivecurve()$HR, 2), lty=2)
          lines(seq(input$T1Mean, 120,length=2), rep(ninetyfivecurve()$HR, 2), lty=2)
        } else if (addfeedback[i]=="95% CI for T"){
          x <- predict(fitcontrol, type = "quantile", p = seq(0.001, 0.999, by=.001))[1,]
          p <- rev(seq(0.001, 0.999, by=.001))
          points(fivecurve()$bigT, p[sum(fivecurve()$bigT > x)], cex=1.5, col="orange", pch=19)
          points(ninetyfivecurve()$bigT, p[sum(ninetyfivecurve()$bigT > x)], cex=1.5, col="orange", pch=19)
        } else if (addfeedback[i]=="Simulation curves"){
          if (CalcSim){
            Curves <- drawsimlines()$linelist
            for (i in 1:10){
              lines(Curves[[(i*2)-1]], Curves[[i*2]], col="purple", lwd=0.25, lty=2)
            }
            CalcSim <- FALSE
          }
        }
      }
    }
    
    
    
  })
  
  
  output$feedback <- renderUI({
    str1 <- paste0("The median survival time for the control curve is ", round(predict(fitcontrol, type="quantile", p = 0.5)[[1]], 1), " months")
    str2 <- paste0("The median survival time for the proposed treatment curve is ", round(proposed()$combinedt[sum(!proposed()$combinedy<0.5)], 1), " months")
    str3 <- paste0("The 95% CI for T is (", round(fivecurve()$bigT, 4), ", ", round(ninetyfivecurve()$bigT, 4), ")")
    str4 <- paste0("The 95% CI for HR is (", round(fivecurve()$HR, 4), ", ", round(ninetyfivecurve()$HR, 4), ")")
    HTML(paste(str1, str2, str3, str4, sep = '<br/>'))
  })
  
  
  
  output$plotAssurance <- renderPlot({
    
    AssFunc <- function(n1, n2){
      assnum <- 50
      assvec <- rep(NA, assnum)
      
      for (i in 1:length(assvec)){
        
        
        lambda2 <- exp(-fitcontrol$coefficients)
        gamma2 <- gamma1
        bigT <- rnorm(1, mean = input$T1Mean, sd = sqrt(input$T1Var))
        HRdist <- estBetaParams(mu = input$T2HRMean , var = sqrt(input$T2HRVar))
        HR <- rbeta(1, HRdist$alpha, HRdist$beta)
        lambda1 <- exp((log(HR)+gamma1*log(lambda2))/gamma1)
        
        #Simulate data for the control group
        simcontrol <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))
        
        
        CP <- exp(-(lambda2*bigT)^gamma2)[[1]]
        u <- runif(n2)
        suppressWarnings(z <- ifelse(u>CP, (1/lambda2)*exp(1/gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(lambda2*bigT)^gamma2+lambda1^gamma1*bigT*gamma1)))))
        
        
        SimTreatment <- data.frame(time = z)
        DataCombined <- data.frame(time = c(simcontrol$time, SimTreatment$time), 
                                   group = c(rep("Control", n1), rep("Treatment", n2)), cens = rep(1, n1+n2))
        test <<- survdiff(Surv(time, cens)~group, data = DataCombined)
        assvec[i] <- test$chisq > qchisq(0.95, 1)
        
      }
      
      return(sum(assvec)/assnum)  
    }
    
    n1n2sum <- input$n1+input$n2
    n1vec <- floor(seq(input$n1*10, floor((1000/n1n2sum)*input$n1), length=20))
    n2vec <- floor(seq(input$n2*10, floor((1000/n1n2sum)*input$n2), length=20))
    assvec <- rep(NA, length(n1vec))
    
    for (i in 1:length(n1vec)){
      assvec[i] <- AssFunc(n1vec[i], n2vec[i])
    }
    
    sumvec <- n1vec+n2vec
    asssmooth <- loess(assvec~sumvec)
    plot(sumvec, predict(asssmooth), type="l", lty=2, ylim=c(0,1),
         xlab="Total sample size", ylab="Assurance")
    
    y <- predict(asssmooth,newdata=seq(10, 1000, by=10))
    y <- y>0.8
    x <- seq(10, 1000, by=0.1)
    print(x[min(which(y == TRUE))])
  })
  
}

shinyApp(ui, server)













