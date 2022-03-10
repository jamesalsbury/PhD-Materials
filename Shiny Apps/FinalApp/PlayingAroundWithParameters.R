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
  
    
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 useShinyjs(),
                 fileInput("uploadSample", "Upload your data", accept = ".rds"),
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
  
  
  inputData <- reactive({
    chosenFile <<- input$uploadSample
    if (is.null(chosenFile)){
      return(NULL)
    } else {
      
      chosenFile <- readRDS(chosenFile$datapath)
      
      #Find the parameters from the control group; lambda2, gamma2
      controldata <- chosenFile[chosenFile$group=="control",]
      weibcontrolfit <- survreg(Surv(time, cens)~1, data = controldata, dist = "weibull")
      
      gamma2 <- as.numeric(exp(-weibcontrolfit$icoef[2]))
      lambda2 <- as.numeric(1/(exp(weibcontrolfit$icoef[1])))
      
      #Extract the treatment data
      treatmentdata <- chosenFile[chosenFile$group=="treatment",]
      
      kmtreatment <- survfit(Surv(time, cens)~1, data = treatmentdata)

      #Optimise lambda1 in the case where we are setting gamma1 = gamma2
      optimfunc1 <- function(par){
        
        gamma1 <- gamma2
        se <- 0
        
        for (i in 1:150){
          if (kmtreatment$time[i]<3){
            y <- exp(-(lambda2*kmtreatment$time[i])^gamma2)
            se <- se + (y-kmtreatment$surv[i])^2
          } else {
            y <- exp(-(lambda2*3)^gamma2-par^gamma1*(kmtreatment$time[i]^gamma1-3^gamma1))
            se <- se + (y-kmtreatment$surv[i])^2
          }
        }
        return(se)
        
      }
      
      optimoutput1 <-  optimize(f = optimfunc1, interval = c(0,1))
      
      
      lambda1only <- optimoutput1$minimum
      
      
      optimfunc2 <- function(par){
        
        se <- 0
        
        for (i in 1:150){
          if (kmtreatment$time[i]<3){
            y <- exp(-(lambda2*kmtreatment$time[i])^gamma2)
            se <- se + (y-kmtreatment$surv[i])^2
          } else {
            y <- exp(-(lambda2*3)^gamma2-par[1]^par[2]*(kmtreatment$time[i]^par[2]-3^par[2]))
            se <- se + (y-kmtreatment$surv[i])^2
          }
        }
        return(se)
        
      }
      
      
      optimoutput2 <-  optim(par = c(0, 1), fn = optimfunc2)
      lambda1 <- optimoutput2$par[1]
      gamma1 <- optimoutput2$par[2]
      
      
      return(list(lambda2 = lambda2, gamma2 = gamma2, lambda1only = lambda1only, lambda1 = lambda1, 
                  gamma1 = gamma1, controltime = controldata$time, controlcens = controldata$cens,
                  treatmenttime = treatmentdata$time, treatmentcens = treatmentdata$cens))
    }
    
  })
  
  
  output$plotSurvival <- renderPlot({
    
    
    if (is.null(inputData())){

    } else {
      controltime <- seq(0, exp((1.527/inputData()$gamma2)-log(inputData()$lambda2))*2, by=0.01)
      controlsurv <- exp(-(inputData()$lambda2*controltime)^inputData()$gamma2)

      plot(controltime, controlsurv, type="l", col="blue", xlab = "Time", ylab = "Survival")
      
      controldata <- data.frame(time = inputData()$controltime, cens =  inputData()$controlcens)
      
      kmcontrolfit <- survfit(Surv(time, cens)~1, data = controldata)
      
      lines(kmcontrolfit, conf.int = F, col="blue")
      
      
      
      #Firstly plot before the changepoint
      treatmenttime1 <- seq(0, 3, by=0.01)
      treatmentsurv1 <- exp(-(inputData()$lambda2*treatmenttime1)^inputData()$gamma2)
      lines(treatmenttime1,treatmentsurv1, col="green")
      
      #Do the Kaplan-Meier plot of the treatment data
      treatmentdata <- data.frame(time = inputData()$treatmenttime, cens =  inputData()$treatmentcens)
      
      kmtreatmentfit <- survfit(Surv(time, cens)~1, data = treatmentdata)
      
      lines(kmtreatmentfit, conf.int = F, col="red")
      
      
      #The case where gamma1 = gamma2
      treatmenttime2 <- seq(3, exp((1.527/inputData()$gamma2)-log(inputData()$lambda2))*2, by=0.01)
      treatmentsurv2 <- exp(-(inputData()$lambda2*3)^inputData()$gamma2 - inputData()$lambda1only^inputData()$gamma2*(treatmenttime2^inputData()$gamma2-3^inputData()$gamma2))
      lines(treatmenttime2, treatmentsurv2, col="red")
      
      #The case where we allow gamma1 to vary
      
      
      treatmenttime3 <- seq(3, exp((1.527/inputData()$gamma2)-log(inputData()$lambda2))*2, by=0.01)
      treatmentsurv3 <- exp(-(inputData()$lambda2*3)^inputData()$gamma2 - inputData()$lambda1^inputData()$gamma1*(treatmenttime3^inputData()$gamma1-3^inputData()$gamma1))
      lines(treatmenttime3, treatmentsurv3, col="red")
      
    }
    
      
  })
  
  
  calculateAssurance1param <- eventReactive(input$drawAssurance, {

    gamma1 <- inputData()$gamma2
    
    assnum <- 100
    assvec <- rep(NA, assnum)
    eventsvec <- rep(NA, assnum)
    controlevents <- rep(NA, assnum)
    treatmentevents <- rep(NA, assnum)


    assFunc <- function(n1, n2){


      for (i in 1:assnum){

        bigT <- 3
        
        lambda1 <- inputData()$lambda1only

        controldata <- data.frame(time = rweibull(n1, inputData()$gamma2, 1/inputData()$lambda2))

        CP <- exp(-(inputData()$lambda2*bigT)^inputData()$gamma2)[[1]]
        u <- runif(n2)

        suppressWarnings(z <- ifelse(u>CP, (1/inputData()$lambda2)*exp(1/inputData()$gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(inputData()$lambda2*bigT)^inputData()$gamma2+lambda1^gamma1*bigT*gamma1)))))

        DataCombined <- data.frame(time = c(controldata$time, z),
                                   group = c(rep("Control", n1), rep("Treatment", n2)))


        DataCombined$time <- DataCombined$time + runif(n1+n2, min = 0, max = input$rectime)

        DataCombined$cens <- DataCombined$time < input$chosenLength

        DataCombined$cens <- DataCombined$cens*1

        test <- survdiff(Surv(time, cens)~group, data = DataCombined)
        assvec[i] <- test$chisq > qchisq(0.95, 1)

        eventsseen <- DataCombined %>%
          filter(time < input$chosenLength)

        controlevents[i] <- sum(eventsseen$group=="Control")

        treatmentevents[i] <- sum(eventsseen$group=="Treatment")

        eventsvec[i] <- sum(DataCombined$cens==1)
      }

      return(list(assvec = mean(assvec), eventvec = mean(eventsvec), controlevents = mean(controlevents), treatmentevents = mean(treatmentevents)))
    }


    samplesizevec <- seq(30, input$numofpatients, length=10)

    n1vec <- floor(input$n1*(samplesizevec/(input$n1+input$n2)))
    n2vec <- ceiling(input$n2*(samplesizevec/(input$n1+input$n2)))
    calcassvec <- rep(NA, length = length(samplesizevec))


    withProgress(message = "Calculating assurance (1/2)", value = 0, {
      for (i in 1:length(n1vec)){
        calcassvec[i] <- assFunc(n1vec[i], n2vec[i])$assvec
        incProgress(1/length(n1vec))
      }
    })


    eventsseen <- assFunc(n1vec[length(samplesizevec)], n2vec[length(samplesizevec)])$eventvec

    controlevents <- assFunc(n1vec[length(samplesizevec)], n2vec[length(samplesizevec)])$controlevents

    treatmentevents <- assFunc(n1vec[length(samplesizevec)], n2vec[length(samplesizevec)])$treatmentevents

    asssmooth <- loess(calcassvec~samplesizevec)

    return(list(calcassvec = calcassvec, asssmooth = asssmooth, samplesizevec = samplesizevec,
                eventsseen = eventsseen, controlevents = controlevents, treatmentevents = treatmentevents))


  })
  
  calculateAssurance2params <- eventReactive(input$drawAssurance, {
    
    gamma1 <- inputData()$gamma1
    
    assnum <- 100
    assvec <- rep(NA, assnum)
    eventsvec <- rep(NA, assnum)
    controlevents <- rep(NA, assnum)
    treatmentevents <- rep(NA, assnum)
    
    
    assFunc <- function(n1, n2){
      
      
      for (i in 1:assnum){
        
        bigT <- 3
        
        lambda1 <- inputData()$lambda1
        
        controldata <- data.frame(time = rweibull(n1, inputData()$gamma2, 1/inputData()$lambda2))
        
        CP <- exp(-(inputData()$lambda2*bigT)^inputData()$gamma2)[[1]]
        u <- runif(n2)
        
        suppressWarnings(z <- ifelse(u>CP, (1/inputData()$lambda2)*exp(1/inputData()$gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(inputData()$lambda2*bigT)^inputData()$gamma2+lambda1^gamma1*bigT*gamma1)))))
        
        DataCombined <- data.frame(time = c(controldata$time, z),
                                   group = c(rep("Control", n1), rep("Treatment", n2)))
        
        
        DataCombined$time <- DataCombined$time + runif(n1+n2, min = 0, max = input$rectime)
        
        DataCombined$cens <- DataCombined$time < input$chosenLength
        
        DataCombined$cens <- DataCombined$cens*1
        
        test <- survdiff(Surv(time, cens)~group, data = DataCombined)
        assvec[i] <- test$chisq > qchisq(0.95, 1)
        
        eventsseen <- DataCombined %>%
          filter(time < input$chosenLength)
        
        controlevents[i] <- sum(eventsseen$group=="Control")
        
        treatmentevents[i] <- sum(eventsseen$group=="Treatment")
        
        eventsvec[i] <- sum(DataCombined$cens==1)
      }
      
      return(list(assvec = mean(assvec), eventvec = mean(eventsvec), controlevents = mean(controlevents), treatmentevents = mean(treatmentevents)))
    }
    
    
    samplesizevec <- seq(30, input$numofpatients, length=10)
    
    n1vec <- floor(input$n1*(samplesizevec/(input$n1+input$n2)))
    n2vec <- ceiling(input$n2*(samplesizevec/(input$n1+input$n2)))
    calcassvec <- rep(NA, length = length(samplesizevec))
    
    
    withProgress(message = "Calculating assurance (2/2)", value = 0, {
      for (i in 1:length(n1vec)){
        calcassvec[i] <- assFunc(n1vec[i], n2vec[i])$assvec
        incProgress(1/length(n1vec))
      }
    })
    
    
    eventsseen <- assFunc(n1vec[length(samplesizevec)], n2vec[length(samplesizevec)])$eventvec
    
    controlevents <- assFunc(n1vec[length(samplesizevec)], n2vec[length(samplesizevec)])$controlevents
    
    treatmentevents <- assFunc(n1vec[length(samplesizevec)], n2vec[length(samplesizevec)])$treatmentevents
    
    asssmooth <- loess(calcassvec~samplesizevec)
    
    return(list(calcassvec = calcassvec, asssmooth = asssmooth, samplesizevec = samplesizevec,
                eventsseen = eventsseen, controlevents = controlevents, treatmentevents = treatmentevents))
    
    
  })

  output$assurancePlot <- renderPlot({

    theme_set(theme_grey(base_size = 12))
    assurancedf1param <- data.frame(x = calculateAssurance1param()$samplesizevec, y = predict(calculateAssurance1param()$asssmooth))
    p1 <- ggplot(data = assurancedf1param) + geom_line(aes(x = x, y = y), linetype="dashed") + xlab("Number of patients") +
      ylab("Assurance") + ylim(0, 1.05)
    
    assurancedf2params <- data.frame(x = calculateAssurance2params()$samplesizevec, y = predict(calculateAssurance2params()$asssmooth))
    
    
    p1 <- p1 + geom_line(data = assurancedf2params, aes(x = x, y = y), colour = "green") 
    print(p1)



  })
  
  

}

shinyApp(ui, server)

##Need to write-up about how all this works
#Tidy up some of the plots, maybe put a legend and stuff?


