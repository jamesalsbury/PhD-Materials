library(SHELF)
library(shiny)
library(survival)
library(shinydashboard)
library(readxl)
library(rsconnect)
library(ggplot2)
library(ggfortify)
library(dplyr)
library(shinyjs)

ui <- fluidPage(
  
  # Application title
  titlePanel("Delayed Treatment Effects - Changing Parameters"),
  
  # sidebarLayout(
  mainPanel(tags$style(type="text/css",
                       ".shiny-output-error { visibility: hidden; }",
                       ".shiny-output-error:before { visibility: hidden; }"
  ),
  
    # Control UI ---------------------------------
    
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 useShinyjs(),
                 numericInput("lambda1", "lambda1", value = 0.05),
                 numericInput("gamma1", "gamma1", value = 1),
                 numericInput("lambda2", "lambda2", value=0.05),
                 numericInput("gamma2", "gamma2", value=1),
                 numericInput("numofpatients", "How many patients could you enrol into the trial?", value=1000),
                 numericInput("rectime", "How long would it take to enrol all of these patients?", value=6),
                 
                 box(width = 10, title = "Ratio of patients in each group?",
                     splitLayout(
                       numericInput("n1", "Control", value=1, min=1),
                       numericInput("n2", "Treatment", value=1, min=1)
                     )
                 ),
                 numericInput("chosenLength", "How long do you want to run the trial for? (Including recruitment time)", value=60),
                 actionButton("drawAssurance", "Produce plot")
               ), 
               mainPanel = mainPanel(
                 plotOutput("plotSurvival"),
                 plotOutput("assurancePlot")
               )
             ),
  )
)
  




server = function(input, output, session) {
  
  
  output$plotSurvival <- renderPlot({

    controltime <- seq(0, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
    controlsurv <- exp(-(input$lambda2*controltime)^input$gamma2)
    controldf <- data.frame(controltime = controltime,
                            controlsurv = controlsurv)
    theme_set(theme_grey(base_size = 12))
    p1 <- ggplot(data=controldf, aes(x=controltime, y=controlsurv)) +
      geom_line(colour="blue") + xlab("Time") + ylab("Survival") + ylim(0,1)
    
    print(p1)


  })
  
  
  
  # calculateAssurance <- eventReactive(input$drawAssurance, {
  #   
  #   gamma1 <- input$gamma2
  #   conc.probs <- matrix(0, 2, 2)
  #   conc.probs[1, 2] <- 0.5
  #   assnum <- 100
  #   assvec <- rep(NA, assnum)
  #   eventsvec <- rep(NA, assnum)
  #   controlevents <- rep(NA, assnum)
  #   treatmentevents <- rep(NA, assnum)
  #   
  #   
  #   assFunc <- function(n1, n2){
  #     
  #     mySample <- data.frame(copulaSample(myfit1(), myfit2(), cp = conc.probs, n = assnum, d = c(input$dist1, input$dist2)))
  #     
  #     
  #     for (i in 1:assnum){
  #       mySample <- data.frame(copulaSample(myfit1(), myfit2(), cp = conc.probs, n = assnum, d = c(input$dist1, input$dist2)))
  #       
  #       bigT <- mySample[i,1]
  #       HR <- mySample[i,2]
  #       
  #       lambda1 <- exp((log(HR)/input$gamma2)+log(input$lambda2))
  #       
  #       controldata <- data.frame(time = rweibull(n1, input$gamma2, 1/input$lambda2))
  #       
  #       CP <- exp(-(input$lambda2*bigT)^input$gamma2)[[1]]
  #       u <- runif(n2)
  #       
  #       suppressWarnings(z <- ifelse(u>CP, (1/input$lambda2)*exp(1/input$gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(input$lambda2*bigT)^input$gamma2+lambda1^gamma1*bigT*gamma1)))))
  #       
  #       DataCombined <- data.frame(time = c(controldata$time, z),
  #                                  group = c(rep("Control", n1), rep("Treatment", n2)))
  #       
  #       
  #       DataCombined$time <- DataCombined$time + runif(n1+n2, min = 0, max = input$rectime)
  #       
  #       DataCombined$cens <- DataCombined$time < input$chosenLength
  #       
  #       DataCombined$cens <- DataCombined$cens*1
  #       
  #       test <- survdiff(Surv(time, cens)~group, data = DataCombined)
  #       assvec[i] <- test$chisq > qchisq(0.95, 1)
  #       
  #       eventsseen <- DataCombined %>%
  #         filter(time < input$chosenLength)
  #       
  #       controlevents[i] <- sum(eventsseen$group=="Control")
  #       
  #       treatmentevents[i] <- sum(eventsseen$group=="Treatment")
  #       
  #       eventsvec[i] <- sum(DataCombined$cens==1)
  #     }
  #     
  #     return(list(assvec = mean(assvec), eventvec = mean(eventsvec), controlevents = mean(controlevents), treatmentevents = mean(treatmentevents)))
  #   }
  #   
  #   
  #   samplesizevec <- seq(30, input$numofpatients, length=10)
  #   
  #   n1vec <- floor(input$n1*(samplesizevec/(input$n1+input$n2)))
  #   n2vec <- ceiling(input$n2*(samplesizevec/(input$n1+input$n2)))
  #   calcassvec <- rep(NA, length = length(samplesizevec))
  #   
  #   
  #   withProgress(message = "Calculating assurance", value = 0, {
  #     for (i in 1:length(n1vec)){
  #       calcassvec[i] <- assFunc(n1vec[i], n2vec[i])$assvec
  #       incProgress(1/length(n1vec))
  #     }
  #   })
  #   
  #   
  #   eventsseen <- assFunc(n1vec[length(samplesizevec)], n2vec[length(samplesizevec)])$eventvec
  #   
  #   controlevents <- assFunc(n1vec[length(samplesizevec)], n2vec[length(samplesizevec)])$controlevents
  #   
  #   treatmentevents <- assFunc(n1vec[length(samplesizevec)], n2vec[length(samplesizevec)])$treatmentevents
  #   
  #   asssmooth <- loess(calcassvec~samplesizevec)
  #   
  #   return(list(calcassvec = calcassvec, asssmooth = asssmooth, samplesizevec = samplesizevec, 
  #               eventsseen = eventsseen, controlevents = controlevents, treatmentevents = treatmentevents))
  #   
  #   
  # })    
  
  # output$assurancePlot <- renderPlot({
  #   
  #   theme_set(theme_grey(base_size = input$fs))
  #   assurancedf <- data.frame(x = calculateAssurance()$samplesizevec, y = predict(calculateAssurance()$asssmooth))
  #   p1 <- ggplot(data = assurancedf) + geom_line(aes(x = x, y = y), linetype="dashed") + xlab("Number of patients") +
  #     ylab("Assurance") + ylim(0, 1.05)
  #   print(p1) 
  #   
  #   
  #   
  # })
  
  

}

shinyApp(ui, server)



