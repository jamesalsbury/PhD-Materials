
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
              numericInput("T1Mean", "Mean", value=10),
              numericInput("T1Var", "Variance", value=0.1, min=0)
            )
        )
      ),
      fluidRow(
        box(width = 14, title = "What is the hazard ratio after the change-point?", 
            splitLayout(
              numericInput("T2HRMean", "Mean", value=0.5),
              numericInput("T2HRVar", "Variance", value=0.001, min=0)
            )
        )
      ),
      # fluidRow(
      #   box(width = 14, title = "How many patients in each group?", 
      #       splitLayout(
      #         numericInput("n1", "Control", value=100, min=1),
      #         numericInput("n2", "Treatment", value=100, min=1)
      #       )
      #   )
      # ),
    ),
    
    mainPanel(
      
      tabsetPanel(type="tabs",
                  tabPanel("Elicitation",
                           plotOutput("plotBestFit"),
                           htmlOutput("control"),
                           htmlOutput("treatment"),
                           htmlOutput("controlparams"),
                           htmlOutput("treatmentparams"),
                           htmlOutput("tparams"),
                           htmlOutput("hazard") ), 
                  tabPanel("Assurance", 
                           #plotOutput("plotsimData"),
                           htmlOutput("assurance"),
                           plotOutput("plotAssurance")))
      
      
    )
  )
)


server <- function(input, output, session) {
  
  #lambda2 <- 0.06
  #gamma2 <- 0.8
  simdata <<- data.frame(time = rweibull(1000, 0.8, 1/0.06), cens = rep(1, 1000))
  fitcontrol <<- survreg(Surv(time, cens)~1, dist="weibull", data = simdata)
  lambda2 <- exp(-fitcontrol$coefficients)
  gamma2 <- 1/fitcontrol$scale
  
  gamma1 <<- gamma2
  findlambda1 <<- reactive({
    lambda1 <- exp((log(input$T2HRMean)+gamma1*log(lambda2))/gamma1)
    return(lambda1)
  })
  
  
  
  output$plotBestFit <- renderPlot({
    
    plot(x = predict(fitcontrol, type = "quantile", p = seq(0.001, 0.999, by=.001))[1,],
         y = rev(seq(0.001, 0.999, by = 0.001)), type="l", xlab="Time (months)", ylab="Survival",
         col = "blue", xlim=c(0,sort(simdata$time)[length(simdata$time)*0.99]))
    
    lambda1 <- findlambda1()
    
    effectt <- seq(0, input$T1Mean, by=0.01)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    lines(effectt, effecty, col="green", lty=3)
    
    aftereffectt <- seq(input$T1Mean, max(simdata$time)*1.1, by=0.01)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*input$T1Mean)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-input$T1Mean^gamma1))
    lines(aftereffectt, aftereffecty, col="red", lwd=1.5, lty=2) 
    
  
    
    estBetaParams <- function(mu, var) {
      alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
      beta <- alpha * (1 / mu - 1)
      return(params = list(alpha = alpha, beta = beta))
    }
    
    lambda2 <- exp(-fitcontrol$coefficients)
    gamma2 <- gamma1
    bigT <- qnorm(c(0.05, 0.95), mean = input$T1Mean, sd = sqrt(input$T1Var))
    HRdist <- estBetaParams(mu = input$T2HRMean , var = sqrt(input$T2HRVar))
    HR <- qbeta(c(0.05, 0.95), HRdist$alpha, HRdist$beta)
    
    
    lambda1 <- exp((log(HR[1])+gamma1*log(lambda2))/gamma1)
    
    
    #Fitting hypothetical lines for the treatment; before changepoint
    effectt <- seq(0, bigT[1], by=0.01)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    lines(effectt, effecty, col="green", lty=3)
    
    #Fitting hypothetical lines for the treatment; after changepoint
    aftereffectt <- seq(bigT[1], max(simdata$time)*1.1, by=0.01)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*bigT[1])^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-bigT[1]^gamma1))
    lines(aftereffectt, aftereffecty, col="red", lty=2)
    
    lambda1 <- exp((log(HR[2])+gamma1*log(lambda2))/gamma1)
    
    
    #Fitting hypothetical lines for the treatment; before changepoint
    effectt <- seq(0, bigT[2], by=0.01)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    lines(effectt, effecty, col="green", lty=3)
    
    #Fitting hypothetical lines for the treatment; after changepoint
    aftereffectt <- seq(bigT[2], max(simdata$time)*1.1, by=0.01)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*bigT[2])^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-bigT[2]^gamma1))
    lines(aftereffectt, aftereffecty, col="red", lty=2)

    
    
    
    
    
    
    
    legend("topright", legend = c("Weibull fit to control data", "Proposed treatment survival curve", "Control + Treatment both Weibull", "Hazard Ratio"),
           col=c("blue", "red", "green", "black"), lty=1:3, cex=0.75)
    
  })
  
  
  # output$control <- renderUI({
  #   withMathJax(paste0("We have parameterised the survival for the control as: $$S_c(t) = \\textrm{exp}\\{-(\\lambda_2t)^{\\gamma_2}\\}$$"))
  # })
  # 
  # output$treatment <- renderUI({
  #   withMathJax(paste0("We have parameterised the survival for the treatment as: $$S_t(t)=\\begin{cases}
  #              \\textrm{exp}\\{-(\\lambda_2t)^{\\gamma_2}\\},  & t\\leq T \\\\
  #              \\textrm{exp}\\{-(\\lambda_2T)^{\\gamma_2} - \\lambda_1^{\\gamma_1}(t^{\\gamma_1}- T^{\\gamma_1})\\}, & t > T
  #              \\end{cases}\\!$$"))
  # })
  # 
  # output$controlparams <- renderUI({
  #   withMathJax(paste0("The parameters seen in the plot above are:$$\\lambda_1 = ",signif(findlambda1(), 2),  "  ,\\gamma_1 = ",signif(gamma1, 2) , "$$"))
  # })
  # 
  # output$treatmentparams <- renderUI({
  #   withMathJax(paste0("$$\\lambda_2 = ", signif(as.numeric(exp(-fitcontrol$coefficients)), 2),  "  ,\\gamma_2 = ",signif(1/fitcontrol$scale, 2) , "$$"))
  # })
  # 
  # output$tparams <- renderUI({
  #   withMathJax(paste0("$$T = ", input$T1Mean, "$$"))
  # })
  
  
  # output$plotsimData <- renderPlot({
  #   
  #   estBetaParams <- function(mu, var) {
  #     alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  #     beta <- alpha * (1 / mu - 1)
  #     return(params = list(alpha = alpha, beta = beta))
  #   }
  #   
  #   lambda2 <- exp(-fitcontrol$coefficients)
  #   gamma2 <- gamma1
  #   bigT <- rnorm(1, mean = input$T1Mean, sd = sqrt(input$T1Var))
  #   HRdist <- estBetaParams(mu = input$T2HRMean , var = sqrt(input$T2HRVar))
  #   HR <- rbeta(1, HRdist$alpha, HRdist$beta)
  #   lambda1 <- exp((log(HR)+gamma1*log(lambda2))/gamma1)
  #   
  #   #Simulate data for the control group
  #   simcontrol <- data.frame(time = rweibull(input$n1, gamma2, 1/lambda2), cens = rep(1, input$n1))
  #   #Plot this on a KM curve
  #   fitcontrolKM <- survfit(Surv(time, cens)~1, data = simcontrol)
  #   plot(fitcontrolKM, conf.int = F, xlim=c(0,sort(simdata$time)[length(simdata$time)*0.99]), ylab="Survival", xlab="Time (months)", col="blue",
  #        main = "Simulated data")
  #   legend("topright", legend = c("KM curve to the control data","KM curve to the treatment data") , col=c("blue", "red"), lty=1)
  #   
  #   
  #   
  #   #Fitting hypothetical lines for the treatment; before changepoint
  #   effectt <- seq(0, bigT, by=0.01)
  #   effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
  #   #lines(effectt, effecty, col="green", lty=3)
  #   
  #   #Fitting hypothetical lines for the treatment; after changepoint
  #   aftereffectt <- seq(bigT, max(simcontrol$time)*1.1, by=0.01)
  #   aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*bigT)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-bigT^gamma1))
  #   #lines(aftereffectt, aftereffecty, col="red", lty=2)
  #   
  #   #Combining both before and after the changepoint
  #   treatmentcurvet <- c(effectt, aftereffectt)
  #   treatmentcurvey <- c(effecty, aftereffecty)
  #   
  #   
  #   #Simulating data for the treatment group
  #   SimTreatmentFunc <<- function(n2){
  #     
  #     split <- seq(1, 0, by=-(1/(n2/10)))
  #     
  #     events <- rep(NA, n2)
  #     for (i in 1:(length(split)-1)){
  #       if (i==1){
  #         events[((i-1)*10+1):(i*10)] <- runif(10, 0, treatmentcurvet[sum(split[i+1]<treatmentcurvey)])
  #       }
  #       else{
  #         events[((i-1)*10+1):(i*10)] <- runif(10, treatmentcurvet[sum(split[i]<treatmentcurvey)], treatmentcurvet[sum(split[i+1]<treatmentcurvey)])
  #       }
  #     }
  #     
  #     data <- data.frame(time = events, cens = rep(1, n2))
  #   }
  #   
  #   SimTreatment <- SimTreatmentFunc(input$n2)
  #   fittreatmentKM <- survfit(Surv(time, cens)~1, data = SimTreatment)
  #   lines(fittreatmentKM, conf.int = F, col="red")
  # })
  # 
  
  output$plotAssurance <- renderPlot({
    
    AssFunc <- function(n1, n2){
      assnum <- 50
      assvec <- rep(NA, assnum)
      
      for (i in 1:length(assvec)){
        
        estBetaParams <- function(mu, var) {
          alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
          beta <- alpha * (1 / mu - 1)
          return(params = list(alpha = alpha, beta = beta))
        }
        
        lambda2 <- exp(-fitcontrol$coefficients)
        gamma2 <- gamma1
        bigT <- rnorm(1, mean = input$T1Mean, sd = sqrt(input$T1Var))
        HRdist <- estBetaParams(mu = input$T2HRMean , var = sqrt(input$T2HRVar))
        HR <- rbeta(1, HRdist$alpha, HRdist$beta)
        lambda1 <- exp((log(HR)+gamma1*log(lambda2))/gamma1)
        
        #Simulate data for the control group
        simcontrol <- data.frame(time = rweibull(n1, gamma2, 1/lambda2))
        
        #Fitting hypothetical lines for the treatment; before changepoint
        effectt <- seq(0, bigT, by=0.01)
        effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
        #lines(effectt, effecty, col="green", lty=3)
        
        #Fitting hypothetical lines for the treatment; after changepoint
        aftereffectt <- seq(bigT, max(simcontrol$time)*1.1, by=0.01)
        aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*bigT)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-bigT^gamma1))
        #lines(aftereffectt, aftereffecty, col="red", lty=2)
        
        #Combining both before and after the changepoint
        treatmentcurvet <- c(effectt, aftereffectt)
        treatmentcurvey <- c(effecty, aftereffecty)
        
        
        #Simulating data for the treatment group
        SimTreatmentFunc <<- function(n2){
          
          split <- seq(1, 0, by=-(1/(n2/10)))
          
          events <- rep(NA, n2)
          for (i in 1:(length(split)-1)){
            if (i==1){
              events[((i-1)*10+1):(i*10)] <- runif(10, 0, treatmentcurvet[sum(split[i+1]<treatmentcurvey)])
            }
            else{
              events[((i-1)*10+1):(i*10)] <- runif(10, treatmentcurvet[sum(split[i]<treatmentcurvey)], treatmentcurvet[sum(split[i+1]<treatmentcurvey)])
            }
          }
          
          data <- data.frame(time = events)
        }
        
        SimTreatment <- SimTreatmentFunc(n2)
        DataCombined <- data.frame(time = c(simcontrol$time, SimTreatment$time), 
                                   group = c(rep("Control", n1), rep("Treatment", n2)), cens = rep(1, n1+n2))
        test <<- survdiff(Surv(time, cens)~group, data = DataCombined)
        assvec[i] <- test$chisq > qchisq(0.95, 1)
        
      }
      
      return(sum(assvec)/assnum)  
    }
    
    n1vec <- seq(10, 500, by=20)
    n2vec <- seq(10, 500, by=20)
    assvec <- rep(NA, length(n1vec))
    
    for (i in 1:length(n1vec)){
      assvec[i] <- AssFunc(n1vec[i], n2vec[i])
    }
    
    asssmooth <- loess(assvec~n1vec)
    plot(n1vec, predict(asssmooth), type="l", lty=2, ylim=c(0,1), xaxt = "n", 
         xlab="Total sample size", ylab="Assurance")
    axis(1, at = seq(0, 500,by = 125), labels = seq(0, 1000, by = 250))
    
  })
  
}

shinyApp(ui, server)


