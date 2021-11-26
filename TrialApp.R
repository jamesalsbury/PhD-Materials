
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
              numericInput("T2HRVar", "Variance", value=0.01, min=0)
            )
        )
      ),
      fluidRow(
        box(width = 14, title = "How many patients in each group?", 
            splitLayout(
              numericInput("n1", "Control", value=100, min=1),
              numericInput("n2", "Treatment", value=100, min=1)
            )
        )
      ),
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
                  tabPanel("Assurance"))
      
      
    )
  )
)


server <- function(input, output, session) {
  

  
  output$plotBestFit <- renderPlot({
    
    
    fitcontrol <<- survreg(Surv(time, cens)~1, dist="weibull", data = simdata)
    plot(x = predict(fitcontrol, type = "quantile", p = seq(0.001, 0.999, by=.001))[1,],
         y = rev(seq(0.001, 0.999, by = 0.001)), type="l", xlab="Time (months)", ylab="Survival",
         col = "blue", xlim=c(0,sort(simdata$time)[length(simdata$time)*0.99]))
    
    lambda2 <- exp(-fitcontrol$coefficients)
    gamma2 <- 1/fitcontrol$scale
    
    gamma1 <<- gamma2
    lambda1 <<-  exp((log(input$T2HRMean)+gamma1*log(lambda2))/gamma1)
    
    
    
    effectt <- seq(0, input$T1Mean, by=0.01)
    effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    lines(effectt, effecty, col="green", lty=3)
    
    aftereffectt <- seq(input$T1Mean, max(simdata$time)*1.1, by=0.01)
    aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*input$T1Mean)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-input$T1Mean^gamma1))
    lines(aftereffectt, aftereffecty, col="red", lwd=1.5, lty=2) 
    
    #Need to draw samples from the distributions for T1, T2, T2HR and T3

    # for (i in 1:50){
    # 
    # 
    #   T1 <- rnorm(1, mean = input$T1Mean, sd = sqrt(input$T1Var))
    # 
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
    # 
    #   
    #   gamma1 <- gamma2
    #   lambda1 <-  exp((log(T2HR)+gamma1*log(lambda2))/gamma1)
    # 
    #   effectt <- seq(0, T1, by=0.01)
    #   effecty <- exp(-((exp(-fitcontrol$coefficients))*effectt)^(1/fitcontrol$scale))
    #   lines(effectt, effecty, col="green", lty=3)
    # 
    #   aftereffectt <- seq(T1, max(simdata$time)*1.1, by=0.01)
    #   aftereffecty <- exp(-(exp(-fitcontrol$coefficients)*T1)^(1/fitcontrol$scale)-(lambda1^gamma1)*(aftereffectt^gamma1-T1^gamma1))
    #   lines(aftereffectt, aftereffecty, col="red", lty=2)
    # }
    
    
    legend("topright", legend = c("Weibull fit to control data", "Proposed treatment survival curve", "Control + Treatment both Weibull", "Hazard Ratio"),
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
  
  output$controlparams <- renderUI({
    withMathJax(paste0("The parameters seen in the plot above are:$$\\lambda_1 = ",signif(lambda1, 2),  "  ,\\gamma_1 = ",signif(gamma1, 2) , "$$"))
  })
  
  output$treatmentparams <- renderUI({
    withMathJax(paste0("$$\\lambda_2 = ", signif(as.numeric(exp(-fitcontrol$coefficients)), 2),  "  ,\\gamma_2 = ",signif(1/fitcontrol$scale, 2) , "$$"))
  })
  
  output$tparams <- renderUI({
    withMathJax(paste0("$$T = ", input$T1Mean, "$$"))
  })
  
}

shinyApp(ui, server)



