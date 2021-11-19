
library(nleqslv)
library(shiny)
library(survival)

ui <- fluidPage(
  
  titlePanel("Cancer survival times - Weibull parameterisation"),
  
  sidebarLayout(
    
    sidebarPanel(
      numericInput("T1", label="When does the treatment begin to take effect (months)?", value=5),
      numericInput("T2", label= "How long after the treatment takes effect is the maximum hazard ratio (months)?", value=5),
      numericInput("T2HR", label= "What is the hazard ratio here?", value=0.5),
      numericInput("T3", label= "How long after the maxmimum hazard ratio does the hazard ratio equal one again (months)?", value=20)
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
    gamma2 <- 2
    
    simdata <<- data.frame(time = rweibull(1000, gamma2, 1/lambda2), cens = rep(1, 1000))
    fitcontrolKM <- survfit(Surv(time, cens)~1, data = simdata)
    plot(fitcontrolKM, conf.int = F, xlim=c(0,max(simdata$time)*1.1), ylab="Survival", xlab="Time (months)", col="blue",
         main = "The historical data for the control")
    legend("topright", legend = "Kaplan-Meier curve to the control data", col="blue", lty=1)
    
    
  })
  
  output$plotBestFit <- renderPlot({
    
    
    fitcontrol <<- survreg(Surv(time, cens)~1, dist="weibull", data = simdata)
    plot(x = predict(fitcontrol, type = "quantile", p = seq(0.001, 0.999, by=.001))[1,],
         y = rev(seq(0.001, 0.999, by = 0.001)), type="l", xlab="Time (months)", ylab="Survival",
         col = "blue", xlim=c(0,max(simdata$time)*1.1))
    
    
    effectt <- seq(0, input$T1, by=0.01)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    lines(effectt, effecty, col="green", lty=3)
    
    lambda2 <- exp(-fitcontrol$coefficients)
    gamma2 <- 1/fitcontrol$scale
    
    fn <- function(x) {

      first <- x[1]*x[2]*(x[1]*(input$T1+input$T2))^(x[2]-1) - input$T2HR*(lambda2*gamma2*(lambda2*(input$T1+input$T2))^(gamma2-1))
      second <- x[1]*x[2]*(x[1]*(input$T1+input$T2+input$T3))^(x[2]-1) - (lambda2*gamma2*(lambda2*(input$T1+input$T2+input$T3))^(gamma2-1))

      return(c(first, second))

    }

    y <- nleqslv(c(1,5), fn)
    lambda1 <<- y[[1]][1]
    gamma1 <<- y[[1]][2]

    aftereffectt <- seq(input$T1, max(simdata$time)*1.1, by=0.01)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*input$T1)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-input$T1^gamma1))
    lines(aftereffectt, aftereffecty, col="red", lty=2)

    legend("topright", legend = c("Weibull fit to control data", "Proposed treatment survival curve", "Control + Treatment both Weibull"),
           col=c("blue", "red", "green"), lty=1:3)
  
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
    withMathJax(paste0("The parameters seen in the plot above are:$$\\lambda_1 = ", signif(lambda1, 2),  "  ,\\gamma_1 = ", signif(gamma1, 2) , "$$"))
  })
  
  output$treatmentparams <- renderUI({
    withMathJax(paste0("$$\\lambda_2 = ", signif(as.numeric(exp(-fitcontrol$coefficients)), 2),  "  ,\\gamma_2 = ",signif(1/fitcontrol$scale, 2) , "$$"))
  })
  
  output$tparams <- renderUI({
    withMathJax(paste0("$$T = ", input$T1, "$$"))
  })
  
  
}

shinyApp(ui, server)



