


library(shiny)
library(survival)

ui <- fluidPage(
  
  titlePanel("Cancer survival times - Weibull parameterisation"),
  
  sidebarLayout(
    
    sidebarPanel(
      numericInput("cp", label="When does the treatment begin to work (months)?", value = 3),
      numericInput("lambda1", label= withMathJax(paste0("$$\\lambda_1$$ (Decreasing value normally shifts line upwards)")), value=0.06, min=0),
      numericInput("gamma1", label= withMathJax(paste0("$$\\gamma_1$$ (Decreasing value normally shifts line upwards)")), value = 0.8, min=0),
      actionButton("reset", label="Fit to control/reset")
    ),
    
    mainPanel(
      plotOutput("plotData"),
      plotOutput("plotBestFit"),
      htmlOutput("control"),
      htmlOutput("treatment"),
      htmlOutput("controlparams"),
      htmlOutput("treatmentparams"),
      htmlOutput("tparams"),
      plotOutput("simulatedData")
    )
  )
)


server <- function(input, output, session) {
  

  output$plotData <- renderPlot({

    lambda2 <- 0.06
    gamma2 <- 0.8

    simdata <<- data.frame(time = rweibull(10000, gamma2, 1/lambda2), cens = rep(1, 10000))
    fitcontrolKM <- survfit(Surv(time, cens)~1, data = simdata)
    plot(fitcontrolKM, conf.int = F, xlim=c(0,sort(simdata$time)[length(simdata$time)*0.99]), ylab="Survival", xlab="Time (months)", col="blue",
          main = "The historical data for the control")
    legend("topright", legend = "Kaplan-Meier curve to the control data", col="blue", lty=1)

  
  
  })
  
  output$plotBestFit <- renderPlot({
     
    
    fitcontrol <<- survreg(Surv(time, cens)~1, dist="weibull", data = simdata)
    lambda2 <<- exp(-fitcontrol$coefficients)
    gamma2 <<- 1/fitcontrol$scale
    plot(x = predict(fitcontrol, type = "quantile", p = seq(0.001, 0.999, by=.001))[1,],
          y = rev(seq(0.001, 0.999, by = 0.001)), type="l", xlab="Time (months)", ylab="Survival",
          col = "blue", xlim=c(0,sort(simdata$time)[length(simdata$time)*0.99]), main="Proposed survival and treatment curves")
                                                   

    effectt <- seq(0, input$cp, by=0.01)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    lines(effectt, effecty, col="green", lty=3)

    aftereffectt <- seq(input$cp, sort(simdata$time)[length(simdata$time)*0.99], by=0.01)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*input$cp)^(1/fitcontrol$scale)-(input$lambda1^input$gamma1)*(aftereffectt^input$gamma1-input$cp^input$gamma1))
    lines(aftereffectt, aftereffecty, col="red", lty=2)
    
  
    t <- seq(input$cp, sort(simdata$time)[length(simdata$time)*0.99], by=0.001)
    HR <- (input$lambda1*input$gamma1*(input$lambda1*t)^(input$gamma1-1))/(lambda2*gamma2*(lambda2*t)^(gamma2-1))
    lines(t, HR)
    t1 <- seq(0, input$cp, by=0.01)
    HR1 <- rep(1, length(t1))
    lines(t1, HR1)
    if (HR[1]<1){
      t3 <- seq(HR[1], 1, by=0.01)
      HR3 <- rep(input$cp, length(t3))
      lines(HR3, t3)
    }
     
    
    legend("topright", legend = c("Weibull fit to control data", "Proposed treatment survival curve", "Control + Treatment both Weibull", "Hazard ratio"),
           col=c("blue", "red", "green", "black"), lty=1:3)
    
    
  })
  
  observeEvent(input$reset, {
    updateNumericInput(session, inputId = "lambda1", value = signif(as.numeric(exp(-fitcontrol$coefficients)), 2))
    updateNumericInput(session, inputId  = "gamma1", value = signif(1/fitcontrol$scale, 2))
  })
  
  
  output$control <- renderUI({
   withMathJax(paste0("We have parameterised the survival for the control as: $$S_c(t) = \\textrm{exp}\\{-(\\lambda_2t)^{\\gamma_2}\\}$$"))
   })
  
  output$treatment <- renderUI({
     withMathJax(paste0("We have parameterised the survival for the treatment as: $$S_t(t)=\\begin{cases}
               \\textrm{exp}\\{-(\\lambda_2t)^{\\gamma_2}\\},  & t\\leq T \\\\
               \\textrm{exp}\\{-(\\lambda_2T)^{\\gamma_2} - \\lambda_1^{\\gamma_1}(t^{\\gamma_1}- T^{\\gamma_1})\\}, & t > T
               \\end{cases}\\!$$"))
  })
  
  output$controlparams <- renderUI({
    withMathJax(paste0("The parameters seen in the plot above are:$$\\lambda_1 = ",input$lambda1,  "  ,\\gamma_1 = ",input$gamma1 , "$$"))
  })
  
  output$treatmentparams <- renderUI({
    withMathJax(paste0("$$\\lambda_2 = ", signif(as.numeric(exp(-fitcontrol$coefficients)), 2),  "  ,\\gamma_2 = ",signif(1/fitcontrol$scale, 2) , "$$"))
  })
  
  output$tparams <- renderUI({
    withMathJax(paste0("$$T = ", input$cp, "$$"))
  })
  
  output$simulatedData <- renderPlot({
    n1 <- 500
    n2 <- 500
    bigT <- input$cp
    lambda1 <- input$lambda1
    gamma1 <- input$gamma1
    
    
    
    simcontrol <- data.frame(time = rweibull(n1, gamma2, 1/lambda2), cens = rep(1, n1))
    #Plot this on a KM curve
    fitcontrolKM <- survfit(Surv(time, cens)~1, data = simcontrol)
    plot(fitcontrolKM, conf.int = F, xlim=c(0,sort(simdata$time)[length(simdata$time)*0.99]), ylab="Survival", xlab="Time (months)", col="blue",
         main = "Simulated data")
    legend("topright", legend = c("KM curve to the control data","KM curve to the treatment data") , col=c("blue", "red"), lty=1)
    
    #Fit a Weibull to the control data
    fitcontrol <- survreg(Surv(time, cens)~1, dist="weibull", data = simcontrol)
    
    
    #Fitting hypothetical lines for the treatment; before changepoint
    effectt <- seq(0, bigT, by=0.001)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    #lines(effectt, effecty, col="green", lty=3)
    
    #Fitting hypothetical lines for the treatment; after changepoint
    aftereffectt <- seq(bigT, max(simcontrol$time)*0.9, by=0.001)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*bigT)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-bigT^gamma1))
    #lines(aftereffectt, aftereffecty, col="red", lty=2)
    
    #Combining both before and after the changepoint
    treatmentcurvet <- c(effectt, aftereffectt)
    treatmentcurvey <- c(effecty, aftereffecty)
    
    
    #Simulating data for the treatment group
    SimTreatmentFunc <- function(n2){
      
      split <- seq(1, 0, by=-(1/(n2/15)))
      
      events <- rep(NA, n2)
      for (i in 1:(length(split)-1)){
        if (i==1){
          events[((i-1)*10+1):(i*10)] <- runif(10, 0, treatmentcurvet[sum(split[i+1]<treatmentcurvey)])
        } 
        else{
          events[((i-1)*10+1):(i*10)] <- runif(10, treatmentcurvet[sum(split[i]<treatmentcurvey)], treatmentcurvet[sum(split[i+1]<treatmentcurvey)])
        }
      }
      
      data <- data.frame(time = events, cens = rep(1, n2))
    }
    
    SimTreatment <- SimTreatmentFunc(n2)
    fittreatmentKM <- survfit(Surv(time, cens)~1, data = SimTreatment)
    lines(fittreatmentKM, conf.int = F, col="red")
  })
  
}

shinyApp(ui, server)





