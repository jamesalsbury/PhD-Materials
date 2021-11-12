library(shiny)
library(survival)

ui <- fluidPage(
  
  titlePanel("Cancer survival times - Weibull parameterisation"),
  
  sidebarLayout(
    
    sidebarPanel(
      numericInput("n1", label="How many patients in the control group?", value = 100),
      numericInput("n2", label="How many patients in the treatment group?", value = 200),
      numericInput("cp", label="When does the treatment begin to work (days)?", value = 50),
      numericInput("percentagealive", label="When the treatment begins to work, what percentage of patients are still alive?", value = 0.5),
      numericInput("rate", label="When the treatment begins to work, what rate do patients die (compared to control?)", value = 2),
      actionButton("go", label="Generate the data")
    ),
    
    mainPanel(
        plotOutput("plotData"),
        plotOutput("plotWeib"),
        htmlOutput("params")
      )
    )
  )




server <- function(input, output) {
  
  moxonidineData <- eventReactive(input$go, {
    control <- data.frame(time = c(runif(input$n1*(1-input$percentagealive), 0, input$cp), 
                runif(input$n1*input$percentagealive, input$cp, input$cp +
                        (input$n1*input$percentagealive*input$cp)/(input$n1*(1-input$percentagealive)))), cens = rep(1, input$n1))
    
    treatment <- data.frame(time = c(runif(input$n2*(1-input$percentagealive), 0, input$cp), 
                  runif(input$n2*input$percentagealive, input$cp, input$cp +
                          (input$rate*input$n2*input$percentagealive*input$cp)/(input$n1*(1-input$percentagealive)))), cens = rep(1, input$n2))
    
    moxonidineData <- list(control, treatment)
    return(moxonidineData)
  },
)
  
  output$plotData <- renderPlot({
    simMoxonidineData <- moxonidineData()
    fitc <- survfit(Surv(time, cens)~1, data = simMoxonidineData[[1]])
    fitt <- survfit(Surv(time, cens)~1, data = simMoxonidineData[[2]])
    plot(fitc, conf.int = F, col=c("blue"), xlim=c(0, max(simMoxonidineData[[2]])), 
         main="Kaplan-Meier curve for the data")
    lines(fitt, conf.int = F, col="red", lty=2)
    legend("topright", legend = c("Control", "Treatment"), lty=1:2, col=c("blue", "red"))
  })
  
  output$plotWeib <- renderPlot({
    simMoxonidineData <- moxonidineData()
    fitcontrol <<- survreg(Surv(time, cens)~1, dist="weibull", data = simMoxonidineData[[1]])
    plot(x = predict(fitcontrol, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,],
         y = rev(seq(0.01, 0.99, by = 0.01)), type="l", xlab="Time", ylab="Survival",
         col = "blue", xlim=c(0, max(simMoxonidineData[[2]])),
         main="Best estimate of the survival curves using Weibull distributions")
    legend("topright", legend = c("Control+Treatment", "Control", "Treatment"), lty=c(1, 1, 2),
           col=c("purple", "blue", "red"))
    
    PredictControl <- predict(fitcontrol, type = "quantile", p = seq(0.01, 0.99, by=.01))[1,]
    for (i in 1:length(PredictControl)){
      if ((PredictControl[i]<input$cp)==F){
        break
      }
    }
    
    probs <- seq(0.01, 0.01*i, by=.01)
    
    
    #The treatment line is the same as the control line up until the changepoint
    lines(x = predict(fitcontrol, type = "quantile", p =probs)[1,],
          y = rev(seq(1-0.01*length(probs), 0.99, by = 0.01)),
          col = "purple")
    
    TreatmentKMCP <- survfit(Surv(time, cens)~1, data = simMoxonidineData[[2]])
    TreatmentKMCPSurv <- summary(TreatmentKMCP)$surv[(input$cp+1):input$n2]
    TreatmentKMCPTime <- summary(TreatmentKMCP)$time[(input$cp+1):input$n2]

    ScaleGuess <<- 160*(input$n2/input$n1)*(input$percentagealive/0.5)*(input$rate/2)*(input$cp/50)
    ShapeGuess <<- 1.9+(0.2*(50/input$cp))*(1.3*(0.5/input$percentagealive))*(0.07*(input$rate/2))
   
    x <- input$cp:floor(max(TreatmentKMCPTime))
    y <- dweibull(x, shape=ShapeGuess, scale=ScaleGuess)
    z <- cumsum(y)
    lines(x, 1-z-max(probs), col="red", lty=2)
  })
  
  
  output$params <-  renderUI({
    simMoxonidineData <- moxonidineData()
    str1 <- paste0("The parameters for the control group: Weibull(scale = ", round(exp(fitcontrol$coefficients), 2),
          ", shape = ", round(1/fitcontrol$scale, 2), ")")
    str2 <- paste0("The parameters for the treatment group when treatment starts to take effect: Weibull(scale = ", round(ScaleGuess, 2),
                   ", shape = ", round(ShapeGuess, 2), ")")
    HTML(paste(str1, str2, sep = '<br/>'))
  })
  

}
shinyApp(ui, server)






