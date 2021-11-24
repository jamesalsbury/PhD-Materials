


library(shiny)
library(survival)

ui <- fluidPage(
  
  titlePanel("Cancer survival times - Weibull parameterisation"),
  
  sidebarLayout(
    
    sidebarPanel(
      numericInput("cp", label="When does the treatment begin to work (months)?", value = 5),
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
      htmlOutput("hazard")
    )
  )
)


server <- function(input, output, session) {
  

  output$plotData <- renderPlot({

    lambda2 <- 0.06
    gamma2 <- 0.8

    simdata <<- data.frame(time = rweibull(10000, gamma2, 1/lambda2), cens = rep(1, 10000))
    fitcontrolKM <- survfit(Surv(time, cens)~1, data = simdata)
    plot(fitcontrolKM, conf.int = F, xlim=c(0,100), ylab="Survival", xlab="Time (months)", col="blue",
          main = "The historical data for the control")
    legend("topright", legend = "Kaplan-Meier curve to the control data", col="blue", lty=1)

    
  })
  
  output$plotBestFit <- renderPlot({
     
    
    fitcontrol <<- survreg(Surv(time, cens)~1, dist="weibull", data = simdata)
    plot(x = predict(fitcontrol, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,],
          y = rev(seq(0.01, 0.99, by = 0.01)), type="l", xlab="Time (months)", ylab="Survival",
          col = "blue", xlim=c(0,100))
                                                   

    effectt <- seq(0, input$cp, by=0.01)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    lines(effectt, effecty, col="green", lty=3)

    aftereffectt <- seq(input$cp, 100, by=0.01)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*input$cp)^(1/fitcontrol$scale)-(input$lambda1^input$gamma1)*(aftereffectt^input$gamma1-input$cp^input$gamma1))
    lines(aftereffectt, aftereffecty, col="red", lty=2)
    
    lambda2 <- exp(-fitcontrol$coefficients)
    gamma2 <- 1/fitcontrol$scale
    t <- seq(input$cp, 100, by=0.001)
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
    updateNumericInput(session, "cp", value = 5)
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
  
  output$hazard <- renderUI({
    # lambda2 <- exp(-fitcontrol$coefficients)
    # gamma2 <- 1/fitcontrol$scale
    # t <- seq(0, 50, by=0.01)
    # HR <- (input$lambda1*input$gamma1*(input$lambda1*t)^(input$gamma1-1))/(lambda2*gamma2*(lambda2*t)^(gamma2-1))
  })
  
}

shinyApp(ui, server)





