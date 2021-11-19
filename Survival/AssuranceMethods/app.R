
library(nleqslv)
library(shiny)
library(survival)

ui <- fluidPage(
  
  titlePanel("Cancer survival times - Weibull parameterisation"),
  
  sidebarLayout(
    
    sidebarPanel(
      numericInput("T1", label="When does the treatment begin to work (months)?", value = 5),
      numericInput("T2", label= "When is the maximum hazard ratio (months)?", value = 10),
      numericInput("T2hazard", label= "What is the hazard ratio here?", value = 0.5),
      numericInput("HazardOne", label= "When does the hazard equal 1 again (months)?", value = 40)
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
    
    simdata <<- data.frame(time = rweibull(1000, gamma2, 1/lambda2), cens = rep(1, 1000))
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
    
    
    effectt <- seq(0, input$T1, by=0.01)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    lines(effectt, effecty, col="green", lty=3)
    
    lambda2 <- exp(-fitcontrol$coefficients)
    gamma2 <- 1/fitcontrol$scale
    
    fn <- function(x) {

      first <- x[1]*x[2]*(x[1]*input$T2)^(x[2]-1) - input$T2hazard*(lambda2*gamma2*(lambda2*input$T2)^(gamma2-1))
      second <- x[1]*x[2]*(x[1]*input$HazardOne)^(x[2]-1) - (lambda2*gamma2*(lambda2*input$HazardOne)^(gamma2-1))

      return(c(first, second))

    }

    y <- nleqslv(c(1,5), fn)
    lambda1 <- y[[1]][1]
    gamma1 <- y[[1]][2]

    aftereffectt <- seq(input$T1, 100, by=0.01)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*input$T1)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-input$T1^gamma1))
    lines(aftereffectt, aftereffecty, col="red", lty=2)

    legend("topright", legend = c("Weibull fit to control data", "Proposed treatment survival curve", "Control + Treatment both Weibull"),
           col=c("blue", "red", "green"), lty=1:3)
  
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
    withMathJax(paste0("$$T = ", input$T1, "$$"))
  })
  
  output$hazard <- renderUI({
    # lambda2 <- exp(-fitcontrol$coefficients)
    # gamma2 <- 1/fitcontrol$scale
    # t <- seq(0, 50, by=0.01)
    # HR <- (input$lambda1*input$gamma1*(input$lambda1*t)^(input$gamma1-1))/(lambda2*gamma2*(lambda2*t)^(gamma2-1))
  })
  
}

shinyApp(ui, server)



