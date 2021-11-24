
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
        box(width = 14, title = "How long after the treatment takes effect is the maximum hazard ratio (months)?", 
            splitLayout(
              numericInput("T2Mean", "Mean", value=5),
              numericInput("T2Var", "Variance", value=0.1, min=0)
            )
        )
      ),
      fluidRow(
        box(width = 14, title = "What is the hazard ratio here?", 
            splitLayout(
              numericInput("T2HRMean", "Mean", value=0.5),
              numericInput("T2HRVar", "Variance", value=0.01, min=0)
            )
        )
      ),
      fluidRow(
        box(width = 14, title = "How long after the maxmimum hazard ratio does the hazard ratio equal one again (months)?", 
            splitLayout(
              numericInput("T3Mean", "Mean", value=20),
              numericInput("T3Var", "Variance", value=0.1, min=0)
            )
        )
      )
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
    
    lambda2 <- exp(-fitcontrol$coefficients)
    gamma2 <- 1/fitcontrol$scale
    
    fn <- function(x) {
      
      first <- x[1]*x[2]*(x[1]*(input$T1Mean+input$T2Mean))^(x[2]-1) - input$T2HRMean*(lambda2*gamma2*(lambda2*(input$T1Mean+input$T2Mean))^(gamma2-1))
      second <- x[1]*x[2]*(x[1]*(input$T1Mean+input$T2Mean+input$T3Mean))^(x[2]-1) - (lambda2*gamma2*(lambda2*(input$T1Mean+input$T2Mean+input$T3Mean))^(gamma2-1))
      
      return(c(first, second))
      
    }
    
    # t1 <- seq(0, input$T1Mean, by=0.01)
    # HR1 <- rep(1, length(t1))
    y <- nleqslv(c(1,5), fn)
    lambda1 <- y[[1]][1]
    gamma1 <- y[[1]][2]
    # t2 <- seq(input$T1Mean, max(simdata$time)*1.1, by=0.01)
    # HR2 <- (lambda1*gamma1*(lambda1*t2)^(gamma1-1))/(lambda2*gamma2*(lambda2*t2)^(gamma2-1))
    # lines(t1, HR1)
    # lines(t2, HR2)
    # t3 <- seq(HR2[1], 1, by=0.01)
    # HR3 <- rep(input$T1Mean, length(t3))
    # lines(HR3, t3)
    
    effectt <- seq(0, input$T1Mean, by=0.01)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    lines(effectt, effecty, col="green", lty=3)
    
    aftereffectt <- seq(input$T1Mean, max(simdata$time)*1.1, by=0.01)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*input$T1Mean)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-input$T1Mean^gamma1))
    lines(aftereffectt, aftereffecty, col="red", lwd=1.5, lty=2) 
    
    #Need to draw samples from the distributions for T1, T2, T2HR and T3
    
    # for (i in 1:10){
    #   
    #   
    #   T1 <- rnorm(1, mean = input$T1Mean, sd = sqrt(input$T1Var))
    #   
    #   T2 <- rnorm(1, mean = input$T2Mean, sd = sqrt(input$T2Var))
    #   
    #   estBetaParams <- function(mu, var) {
    #     alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    #     beta <- alpha * (1 / mu - 1)
    #     return(params = list(alpha = alpha, beta = beta))
    #   }
    #   
    #   EstT2HR <- estBetaParams(input$T2HRMean, input$T2HRVar)
    #   
    #   T2HR <- rbeta(1, EstT2HR$alpha, EstT2HR$beta)
    #   
    #   T3 <- rnorm(1, mean = input$T3Mean, sd = sqrt(input$T3Var))
    #   
    #   
    #   y <- nleqslv(c(1,5), fn)
    #   lambda1 <- y[[1]][1]
    #   gamma1 <- y[[1]][2]
    #   
    #   effectt <- seq(0, T1, by=0.01)
    #   effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    #   lines(effectt, effecty, col="green", lty=3)
    #   
    #   aftereffectt <- seq(T1, max(simdata$time)*1.1, by=0.01)
    #   aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*T1)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-T1^gamma1))
    #   lines(aftereffectt, aftereffecty, col="red", lty=2) 
    # }
    
    
    legend("bottomleft", legend = c("Weibull fit to control data", "Proposed treatment survival curve", "Control + Treatment both Weibull", "Hazard Ratio"),
           col=c("blue", "red", "green", "black"), lty=1:3, cex=0.75)
    

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
  
  
}

shinyApp(ui, server)



