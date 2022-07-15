library(shiny)
library(shinyjs)
library(shinydashboard)
library(ggplot2)
library(survival)
library(tidyverse)
ui <- fluidPage(
  
  # Application title
  titlePanel("Delayed Treatment Effects - Interim Analyses"),
  
  # sidebarLayout(
  mainPanel(tags$style(type="text/css",
                       ".shiny-output-error { visibility: hidden; }",
                       ".shiny-output-error:before { visibility: hidden; }"
  ),
  
  tabsetPanel(
    # Set up UI ---------------------------------
    
    tabPanel("Set up", 
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 numericInput("lambda2", "lambda2", value=0.08),
                 numericInput("gamma2", "gamma2", value=0.8),
                 box(width = 10, title = "Delay",
                     splitLayout(
                       numericInput("DelayMean", "Mean", value=6, min=0),
                       numericInput("DelaySD", "SD", value=1, min=0)
                     )
                 ),
                 box(width = 10, title = "Post-delay HR",
                     splitLayout(
                       numericInput("HRa", "a", value=10, min=1),
                       numericInput("HRb", "b", value=6, min=1)
                     )
                 ),
                 checkboxGroupInput("showfeedback", "Add to plot", choices = c("Median survival line", "95% CI for T", "CI for Treatment Curve (0.1 and 0.9)")),
                 
               ), 
               mainPanel = mainPanel(
                 plotOutput("plotFeedback")
               )
             ),
    ),
    
    # Assurance UI ---------------------------------
    tabPanel("Assurance",
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 useShinyjs(),
                 numericInput("numofpatients", "How many patients could you enrol into the trial?", value=1000),
                 numericInput("chosenLength", "How long do you want to run the trial for?", value=60),
                 actionButton("drawAssurance", "Produce plot"),
                 numericInput("chosenAss", "What assurance do you want to aim for?", value=0.8),
                 numericInput("chosenIA", "When do you want to perform the interim analysis?", value=0.75)
               ), 
               mainPanel = mainPanel(
                 plotOutput("assurancePlot"),
                 htmlOutput("assuranceText")
               )
             ),
             
    ),
  ),
  )
)


server = function(input, output, session) {
  

  # Functions for the Set up tab ---------------------------------
  
  drawsimlines <- reactive({
    
    gamma1 <- input$gamma2
    mySample <- data.frame(t = rnorm(500, input$DelayMean, sd = input$DelaySD), HR = rbeta(500, input$HRa, input$HRb))
    
    time <- seq(0, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
    SimMatrix <- matrix(NA, nrow = 500, ncol=length(time))
    
    for (i in 1:500){
      bigT <- mySample[i,1]
      HR <- mySample[i,2]
      lambda1 <- exp((log(HR)/input$gamma2)+log(input$lambda2))
      
      controltime <- seq(0, bigT, by=0.01)
      controlsurv <- exp(-(input$lambda2*controltime)^input$gamma2)
      
      treatmenttime <- seq(bigT, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
      treatmentsurv <- exp(-(input$lambda2*bigT)^input$gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))
      
      timecombined <- c(controltime, treatmenttime)[1:length(time)]
      survcombined <- c(controlsurv, treatmentsurv)[1:length(time)]
      
      
      SimMatrix[i,] <- survcombined
      
    }
    
    lowerbound <- rep(NA, length(time))
    upperbound <- rep(NA, length(time))
    for (j in 1:length(time)){
      lowerbound[j] <- quantile(SimMatrix[,j], 0.1)
      upperbound[j] <- quantile(SimMatrix[,j], 0.9)
    }
    
    return(list(lowerbound=lowerbound, upperbound=upperbound, time=time))
    
  })
  
  output$plotFeedback <- renderPlot({
    
    gamma1 <- input$gamma2
    controltime <- seq(0, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
    controlcurve <- exp(-(input$lambda2*controltime)^input$gamma2)
    controldf <- data.frame(controltime = controltime, controlcurve = controlcurve)
    p1 <- ggplot(data=controldf, aes(x=controltime, y=controlcurve)) +
      geom_line(colour="blue") + xlab("Time") + ylab("Survival") + ylim(0,1)

      bigTMedian <- input$DelayMean
      HRMedian <- (input$HRa)/(input$HRa+input$HRb)
      lambda1 <- exp((log(HRMedian)/input$gamma2)+log(input$lambda2))

      treatmenttime1 <- seq(0, bigTMedian, by=0.01)
      treatmentsurv1 <- exp(-(input$lambda2*treatmenttime1)^input$gamma2)
      treatmenttime1df <- data.frame(treatmenttime1 = treatmenttime1, treatmentsurv1 = treatmentsurv1)
      p1 <-  p1 + geom_line(data = treatmenttime1df, aes(x = treatmenttime1, y = treatmentsurv1), colour = "green")


      treatmenttime2 <- seq(bigTMedian, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
      treatmentsurv2 <- exp(-(input$lambda2*bigTMedian)^input$gamma2 - lambda1^gamma1*(treatmenttime2^gamma1-bigTMedian^gamma1))
      treatmenttime2df <- data.frame(treatmenttime2 = treatmenttime2, treatmentsurv2 = treatmentsurv2)
      p1 <-  p1 + geom_line(data = treatmenttime2df, aes(x = treatmenttime2, y = treatmentsurv2), colour = "red")
      
      print(p1)
      
      
      addfeedback <- input$showfeedback 
      
      if (!is.null(addfeedback)){
        for (i in 1:length(addfeedback)){
          if (addfeedback[i]=="Median survival line"){
            if (exp(-(input$lambda2*bigTMedian)^input$gamma2)<0.5){
              mediandf <- data.frame(x = seq(0, controltime[sum(controlcurve>0.5)], length=2), y = rep(0.5, 2))
              mediandf1 <- data.frame(x = rep(controltime[sum(controlcurve>0.5)], 2), y = seq(0, 0.5, length=2))
              p1 <- p1 + geom_line(data = mediandf, aes(x = x, y=y), linetype = "dashed") + geom_line(data = mediandf1, aes(x = x, y=y), linetype="dashed") 
            } else {
              mediandf <- data.frame(x = seq(0, treatmenttime2[sum(treatmentsurv2>0.5)], length=2), y = rep(0.5, 2))
              mediandf1 <- data.frame(x = rep(treatmenttime2[sum(treatmentsurv2>0.5)], 2), y = seq(0, 0.5, length=2))
              mediandf2 <- data.frame(x = rep(controltime[sum(controlcurve>0.5)], 2), y = seq(0, 0.5, length=2))
              p1 <- p1 + geom_line(data = mediandf, aes(x = x, y=y), linetype = "dashed") + geom_line(data = mediandf1, aes(x = x, y=y), linetype="dashed") +
                geom_line(data = mediandf2, aes(x = x, y=y), linetype="dashed")
            }
          } else if (addfeedback[i]=="95% CI for T"){
            p1 <- p1 + geom_point(aes(x =  qnorm(0.025, input$DelayMean, sd = input$DelaySD), y = controlcurve[sum(controltime< qnorm(0.025, input$DelayMean, sd = input$DelaySD))]), colour="orange", size = 4) +
              geom_point(aes(x =  qnorm(0.975, input$DelayMean, sd = input$DelaySD), y = controlcurve[sum(controltime< qnorm(0.975, input$DelayMean, sd = input$DelaySD))]), colour="orange", size = 4)
          } else if (addfeedback[i]=="CI for Treatment Curve (0.1 and 0.9)"){
            #Function to draw simulated lines
            simlineslower <- data.frame(x = drawsimlines()$time, y = drawsimlines()$lowerbound)
            simlinesupper <- data.frame(x = drawsimlines()$time, y = drawsimlines()$upperbound)
            p1 <- p1 + geom_line(data = simlineslower, aes(x=x, y=y), linetype="dashed")+
              geom_line(data = simlinesupper, aes(x=x, y=y), linetype="dashed")
          }
      }
      print(p1)
    }
  })
  

  
  # Functions for the Assurance tab ---------------------------------
  
  
  calculateAssurance <- eventReactive(input$drawAssurance, {
    
    gamma1 <- input$gamma2
    assnum <- 200
    
    
    assFunc <- function(n1, n2){
      
      assvec <- rep(NA, assnum)
      eventsvec <- rep(NA, assnum)
      mySample <- data.frame(t = rnorm(assnum, input$DelayMean, sd = input$DelaySD), HR = rbeta(assnum, input$HRa, input$HRb))
      
      
      for (i in 1:assnum){
        
        bigT <- mySample[i,1]
        HR <- mySample[i,2]
        
        lambda1 <- exp((log(HR)/input$gamma2)+log(input$lambda2))
        
        controldata <- data.frame(time = rweibull(n1, input$gamma2, 1/input$lambda2))
        
        CP <- exp(-(input$lambda2*bigT)^input$gamma2)[[1]]
        u <- runif(n2)
        
        suppressWarnings(z <- ifelse(u>CP, (1/input$lambda2)*exp(1/input$gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(input$lambda2*bigT)^input$gamma2+lambda1^gamma1*bigT*gamma1)))))
        
        DataCombined <- data.frame(time = c(controldata$time, z),
                                   group = c(rep("Control", n1), rep("Treatment", n2)))
        
        
        DataCombined$cens <- DataCombined$time < input$chosenLength
        
        DataCombined$cens <- DataCombined$cens*1
        
        
        
        if (sum(DataCombined$cens)==(n1+n2)){
          
        } else {
          DataCombined[DataCombined$cens==0,]$time <- input$chosenLength
        }
        

        test <- survdiff(Surv(time, cens)~group, data = DataCombined)
        assvec[i] <- test$chisq > qchisq(0.95, 1)
        
        eventsvec[i] <- sum(DataCombined$cens==1)
      }
      
      return(list(assvec = mean(assvec), eventvec = mean(eventsvec)))
    }
    
    
    samplesizevec <- seq(30, input$numofpatients, length=10)
    
    n1vec <- round(samplesizevec/2)
    n2vec <- round(samplesizevec/2)
    calcassvec <- rep(NA, length = length(samplesizevec))
    eventvec <- rep(NA, length(samplesizevec))
    
    withProgress(message = "Calculating assurance", value = 0, {
      for (i in 1:length(samplesizevec)){
        output <- assFunc(n1vec[i], n2vec[i])
        calcassvec[i] <- output$assvec
        eventvec[i] <- output$eventvec
        incProgress(1/length(n1vec))
      }
    })
    
    asssmooth <- loess(calcassvec~samplesizevec)
    eventsmooth <- loess(eventvec~samplesizevec)
    
    return(list(calcassvec = calcassvec, asssmooth = asssmooth, samplesizevec = samplesizevec, 
                eventsmooth = eventsmooth))
    
    
  })    
  
  output$assurancePlot <- renderPlot({
    
    assurancedf <- data.frame(x = calculateAssurance()$samplesizevec, y = predict(calculateAssurance()$asssmooth))
    p1 <- ggplot(data = assurancedf) + geom_line(aes(x = x, y = y), linetype="dashed") + xlab("Number of patients") +
      ylab("Assurance") + ylim(0, 1.05)
    print(p1) 
    
    
    
  })
  
  
  output$assuranceText  <- renderUI({
    
    for (i in 50:1000){
      if (predict(calculateAssurance()$asssmooth, newdata = i)>input$chosenAss){
        break
      }
    }

    npatients <- round_any(i, 2)

    events <- round(predict(calculateAssurance()$eventsmooth, newdata = npatients))


    str1 <- paste0("For ", input$chosenAss*100, "% assurance and for a trial which lasts ",input$chosenLength ," months, we require ", npatients, " patients.")
    str2 <- paste0("On average, ", events, " events will be seen.")
    str3 <- paste0("An IA will be performed at ", input$chosenIA*100, "% of the information fraction, this is when ", events*input$chosenIA, " events have occured.")
    str4 <- paste0("We will invesigate the optimal time to perform an IA when the delay ranges from ", round(qnorm(0.025, input$DelayMean, sd = input$DelaySD),2), " to ", round(qnorm(0.975, input$DelayMean, sd = input$DelaySD),2),
                   " months and the post-delay HR ranges from ", round(qbeta(0.025, input$HRa, input$HRb),2), " to ", round(qbeta(0.975, input$HRa, input$HRb), 2))
    HTML(paste(str1, str2, str3, str4, sep = '<br/>'))
    
    
  })

  
}

shinyApp(ui, server)
