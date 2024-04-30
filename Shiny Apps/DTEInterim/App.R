library(shiny)
library(shinyjs)
 library(survival)
library(doParallel)
library(ggplot2)
library(DT)
# library(dplyr)
# library(rjags)
# library(magrittr)
# library(future.apply)
#library(foreach)
# library(future)
# library(promises)
# library(ipc)


source("functions.R")

# UI definition
ui <- fluidPage(
  
  tags$style(HTML("
    .error-message {
      color: red;
      display: none;
    }
  ")),
  
  
  withMathJax(),
  shinyjs::useShinyjs(),
  
  # Application title
  titlePanel("Bayesian Futility Analyses: Delayed Treatment Effects"),
  
  mainPanel(
    tabsetPanel(
      # Data UI ---------------------------------
      tabPanel("Data", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   numericInput("numPatients", 'Number of Patients (in each group, 1:1)', value=340, min=0),
                   numericInput("numEvents", 'Number of Events', value=512, min=0),
                   numericInput("lambdac", '\\( \\lambda_c \\)', value = round(log(2)/12, 4)),
                   numericInput("recTime", 'Recruitment time', value = 34),
                   fileInput("samplesFile", "Upload samples")
                 ), 
                 mainPanel = mainPanel(
                   textOutput("P_S"),
                   textOutput("P_DTE"),
                   plotOutput("TDist"),
                   plotOutput("HRDist")
                   
                 )
               )
      ),
      # One Look UI ---------------------------------
      
      tabPanel("One Look", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   wellPanel(
                   numericInput("OneLookLB", "IA1, from:", value = 0.2),
                   numericInput("OneLookUB", "to:", value = 0.8),
                   numericInput("OneLookBy", "by:", value = 0.1),
                   textOutput("OneLookText")),
                   numericInput("OneLookHR", "Stop if observed HR > ", value = 1),
                   actionButton("calcFutilityOneLook", label  = "Calculate", disabled = T)
                 ), 
                 mainPanel = mainPanel(
                  tableOutput("futilityTable"),
                  tableOutput("finalAssTable"),
                  plotOutput("oneLookROCPlot"),
                  plotOutput("oneLookPlotDuration"),
                  plotOutput("oneLookPlotSS")
                 )
               )
      ),
      
      # Two Looks UI ---------------------------------
      
      tabPanel("Two Looks", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   wellPanel(
                   numericInput("TwoLooksLB1", "IA1, from:", value = 0.2),
                   numericInput("TwoLooksUB1", "to:", value = 0.8),
                   numericInput("TwoLooksBy1", "by:", value = 0.2),
                   numericInput("TwoLooksHR1", "Stop if observed HR > ", value = 1)),
                   wellPanel(
                    numericInput("TwoLooksLB2", "IA2, from:", value = 0.3),
                   numericInput("TwoLooksUB2", "to:", value = 0.9),
                   numericInput("TwoLooksBy2", "by:", value = 0.2),
                   numericInput("TwoLooksHR2", "Stop if observed HR > ", value = 1)),
                   textOutput("TwoLooksText1"),
                   textOutput("TwoLooksText2"),
                   div(id = "twoLooksErrorMessage", class = "error-message", textOutput("twoLooksErrorMessage")),
                   actionButton("calcFutilityTwoLooks", label  = "Calculate", disabled = T)
                 ), 
                 mainPanel = mainPanel(
                   hidden(uiOutput("AssText")),
                   tableOutput("twoLooksAss"),
                   hidden(uiOutput("SSText")),
                   tableOutput("twoLooksSS"),
                   hidden(uiOutput("DurationText")),
                   tableOutput("twoLooksDuration"),
                   hidden(uiOutput("IATime1Text")),
                   tableOutput("IATime1Table"),
                   hidden(uiOutput("IATime2Text")),
                   tableOutput("IATime2Table"),
                   hidden(uiOutput("PercentStopText")),
                   tableOutput("PercentStopTable"),
                   hidden(uiOutput("PercentStopLook1Text")),
                   tableOutput("PercentStopLook1Table"),
                   hidden(uiOutput("PercentStopLook2Text")),
                   tableOutput("PercentStopLook2Table"),
                   hidden(uiOutput("falselyStopLook1Text")),
                   tableOutput("falselyStopLook1Table"),
                   hidden(uiOutput("correctlyStopLook1Text")),
                   tableOutput("correctlyStopLook1Table"),
                   hidden(uiOutput("falselyContinueLook1Text")),
                   tableOutput("falselyContinueLook1Table"),
                   hidden(uiOutput("correctlyContinueLook1Text")),
                   tableOutput("correctlyContinueLook1Table"),
                   hidden(uiOutput("falselyStopLook2Text")),
                   tableOutput("falselyStopLook2Table"),
                   hidden(uiOutput("correctlyStopLook2Text")),
                   tableOutput("correctlyStopLook2Table"),
                   hidden(uiOutput("falselyContinueLook2Text")),
                   tableOutput("falselyContinueLook2Table"),
                   hidden(uiOutput("correctlyContinueLook2Text")),
                   tableOutput("correctlyContinueLook2Table"),
                   hidden(uiOutput("falselyStopTotalText")),
                   tableOutput("falselyStopTotalTable"),
                   hidden(uiOutput("correctlyStopTotalText")),
                   tableOutput("correctlyStopTotalTable"),
                   hidden(uiOutput("falselyContinueTotalText")),
                   tableOutput("falselyContinueTotalTable"),
                   hidden(uiOutput("correctlyContinueTotalText")),
                   tableOutput("correctlyContinueTotalTable"),
                   tableOutput("finalAssTable2Looks"),
                   tableOutput("proposedTable2Looks"),
                   plotOutput("rocCurveTwoLooks1"),
                   plotOutput("rocCurveTwoLooks2"),
                   plotOutput("rocCurveTwoLooksTotal"),
                   plotOutput("twoLooksPlotDuration"),
                   plotOutput("twoLooksPlotSS")
                   
                 )
               )
      ),
      
      
      # Bayesian UI ---------------------------------
      
      tabPanel("Bayesian", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   numericInput("IFBayesian", "Information Fraction", value = 0.5),
                   actionButton("calcFutilityBayesian", label  = "Calculate", disabled = T)
                 ), 
                 mainPanel = mainPanel(
                   plotOutput("BayesianPlot"),
                   tableOutput("BayesianSS"),
                   tableOutput("BayesianDuration")
                   
                 )
               )
      )
      
      
    ), style='width: 1000px; height: 600px'
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Data Logic ---------------------------------
  
  data <- reactiveValues(treatmentSamplesDF = NULL)
  
  observeEvent(input$samplesFile, {
    inFile  <- input$samplesFile
    if (!is.null(inFile)) {
      # Read the uploaded file into treatmentSamplesDF
      data$treatmentSamplesDF <-  readRDS(inFile$datapath)

       
      output$TDist <- renderPlot({
        SHELF::plotfit(data$treatmentSamplesDF$fit1, d = data$treatmentSamplesDF$d[1])
      })

      output$HRDist <- renderPlot({
        SHELF::plotfit(data$treatmentSamplesDF$fit2, d = data$treatmentSamplesDF$d[2])
      })

      output$P_S <- renderText({
        paste0("P_S is: ", data$treatmentSamplesDF$P_S)
      })

      output$P_DTE <- renderText({
        paste0("P_DTE is: ", data$treatmentSamplesDF$P_DTE)
      })
    }
  })
  
  observe({
    if (is.null(data$treatmentSamplesDF)) {
      updateActionButton(session, "calcFutilityOneLook", disabled = TRUE)
      updateActionButton(session, "calcFutilityTwoLooks", disabled = TRUE)
      updateActionButton(session, "calcFutilityBayesian", disabled = TRUE)
    } else {
      updateActionButton(session, "calcFutilityOneLook", disabled = FALSE)
      updateActionButton(session, "calcFutilityTwoLooks", disabled = FALSE)
      updateActionButton(session, "calcFutilityBayesian", disabled = FALSE)
    }
  })
  
  
  
  # One Look Logic ---------------------------------
  
  output$OneLookText <- renderText({
    OneLookSeq <- seq(input$OneLookLB, input$OneLookUB, by = input$OneLookBy)
    value <- paste0("We look at: ", paste(OneLookSeq, collapse = " "))
    return(value)
  })
  
  observeEvent(input$calcFutilityOneLook, {
    NRep <- 500
    futilityVec <- seq(input$OneLookLB, input$OneLookUB, by = input$OneLookBy)
    
    iterationArray <- array(NA,  dim = c(3, length(futilityVec)+1, NRep))
    
    conc.probs <- matrix(0, 2, 2)
    conc.probs[1, 2] <- 0.5
    
    treatmentSamplesDF <- SHELF::copulaSample(data$treatmentSamplesDF$fit1, data$treatmentSamplesDF$fit2,
                                              cp = conc.probs, n = 1e4, d = data$treatmentSamplesDF$d)
    
    withProgress(message = 'Calculating', value = 0, {
      for (i in 1:NRep){
        
        #Compute treatment times
        HRStar <- sample(treatmentSamplesDF[,2], 1)
        bigT <- sample(treatmentSamplesDF[,1], 1)
        
        #Simulate control and treatment data
        dataCombined <- SimDTEDataSet(input$numPatients, input$lambdac, bigT, HRStar, input$recTime)  
        
        #Perform futility look at different Information Fractions
        finalDF <- CensFunc(dataCombined, input$numEvents)
        test <- survdiff(Surv(survival_time, status)~group, data = finalDF$dataCombined)
        coxmodel <- coxph(Surv(survival_time, status)~group, data = finalDF$dataCombined)
        deltad <- as.numeric(exp(coef(coxmodel)))
        
        iterationArray[1, length(futilityVec)+1, i] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
        iterationArray[2, length(futilityVec)+1, i] <- finalDF$censTime
        iterationArray[3, length(futilityVec)+1, i] <- finalDF$SS
        
        
        for (k in 1:length(futilityVec)){
          futilityCens <- CensFunc(dataCombined, input$numEvents*futilityVec[k])
          futilityLook <- interimLookFunc(futilityCens$dataCombined, input$OneLookHR)
          
          iterationArray[1, k, i] <- futilityLook
          iterationArray[2, k, i] <- futilityCens$censTime
          iterationArray[3, k, i] <- futilityCens$SS
          
        }
        
        incProgress(1/NRep)
      }
    })
    
    
    powerDF <- data.frame(matrix(NA, nrow = NRep, ncol = futilityVec))
    durationDF <- data.frame(matrix(NA, nrow = NRep, ncol = futilityVec))
    ssDF <- data.frame(matrix(NA, nrow = NRep, ncol = futilityVec))
    iaTimeDF <- data.frame(matrix(NA, nrow = NRep, ncol = futilityVec))
    stopDF <- data.frame(matrix(NA, nrow = NRep, ncol = futilityVec))
    falselyStop <- rep(NA, length(futilityVec))
    correctlyStop <- rep(NA, length(futilityVec))
    falselyContinue <- rep(NA, length(futilityVec))
    correctlyContinue <- rep(NA, length(futilityVec))
    
    
    for (i in 1:NRep){
      for (k in 1:length(futilityVec)){
        powerDF[i, k] <- iterationArray[1,k,i]*iterationArray[1,length(futilityVec)+1,i]
        iaTimeDF[i,k] <- iterationArray[2,k,i]
        stopDF[i,k] <- iterationArray[1,k,i]

        if (iterationArray[1,k,i]==0){
          durationDF[i, k] <- iterationArray[2,k,i]
          ssDF[i, k] <- iterationArray[3,k,i]
        } else {
          durationDF[i, k] <- iterationArray[2,length(futilityVec)+1,i]
          ssDF[i, k] <- iterationArray[3,length(futilityVec)+1,i]
        }
      }
    }
    
    #Calculating falsely stopping
    for (k in 1:length(futilityVec)){
      indices <- which(iterationArray[1, k, ] == 0)
      selected_layers <- iterationArray[,,indices]
      if (length(dim(selected_layers))==2){
        falselyStop[k] <- mean(selected_layers[1, length(futilityVec)+1])
      } else {
        if (dim(selected_layers)[3]==0){
          falselyStop[k] <- NA
        } else{
          falselyStop[k] <- mean(selected_layers[1, length(futilityVec)+1, ])
        }
      }
    }
    
    #Calculating correctly stopping
    for (k in 1:length(futilityVec)){
      indices <- which(iterationArray[1, k, ] == 0)
      selected_layers <- iterationArray[,,indices]
      if (length(dim(selected_layers))==2){
        correctlyStop[k] <- 1-mean(selected_layers[1, length(futilityVec)+1])
      } else {
        if (dim(selected_layers)[3]==0){
          correctlyStop[k] <- NA
        } else{
          correctlyStop[k] <- 1-mean(selected_layers[1, length(futilityVec)+1, ])
        }
      }
    }
    
    #Calculating falsely continuing
    for (k in 1:length(futilityVec)){
      indices <- which(iterationArray[1, k, ] == 1)
      selected_layers <- iterationArray[,,indices]
      if (length(dim(selected_layers))==2){
        falselyContinue[k] <- 1-mean(selected_layers[1, length(futilityVec)+1])
      } else {
        if (dim(selected_layers)[3]==0){
          falselyContinue[k] <- NA
        } else{
          falselyContinue[k] <- 1-mean(selected_layers[1, length(futilityVec)+1, ])
        }
      }
    }
    
    #Calculating correctly continuing
    for (k in 1:length(futilityVec)){
      indices <- which(iterationArray[1, k, ] == 1)
      selected_layers <- iterationArray[,,indices]
      if (length(dim(selected_layers))==2){
        correctlyContinue[k] <- mean(selected_layers[1, length(futilityVec)+1])
      } else {
        if (dim(selected_layers)[3]==0){
          correctlyContinue[k] <- NA
        } else{
          correctlyContinue[k] <- mean(selected_layers[1, length(futilityVec)+1, ])
        }
      }
    }
    
    
    
    
    
    
    futilityDF <- data.frame(IF = futilityVec,
                             IATime = colMeans(iaTimeDF),
                             Assurance = colMeans(powerDF),
                             Stop = 1-colMeans(stopDF),
                             falselyStop = falselyStop,
                             correctlyStop = correctlyStop,
                             falselyContinue = falselyContinue,
                             correctlyContinue = correctlyContinue,
                             FPR = 1 - correctlyStop/(correctlyStop+falselyStop),
                             TPR = correctlyContinue/(correctlyContinue+falselyContinue),
                             Duration = colMeans(durationDF),
                             SS = colMeans(ssDF))
    
    #Continuing is positive!
    
    colnames(futilityDF) <- c("Information Fraction", "IA Time", "Assurance", "% stop",
                              "Falsely Stop", "Correctly Stop", "Falsely Continue", "Correctly Continue",
                              "False Positive Rate", "True Positive Rate",
                               "Duration", "Sample Size")
    
    
    output$futilityTable <- renderTable({
      futilityDF
    }, digits = 3)
    
    
    output$finalAssTable <- renderTable({
      
      FinalAss <- mean(iterationArray[1, length(futilityVec)+1, ])
      FinalDuration <- mean(iterationArray[2, length(futilityVec)+1, ])
      FinalSS <- mean(iterationArray[3, length(futilityVec)+1, ])
      
      
      FinalAss <- data.frame(Assurance = FinalAss, Duration = FinalDuration, SS = FinalSS)
      
      colnames(FinalAss) <- c("Assurance", "Duration", "Sample Size")
      
      FinalAss
    }, digits = 3)
    
    output$oneLookROCPlot <- renderPlot({
      
      sens <- correctlyContinue/(correctlyContinue+falselyContinue)
      spec <- correctlyStop/(correctlyStop+falselyStop)
      
      plot(1-spec, sens, ylim = c(0, 1), xlim = c(0,1), ylab = "True Positive Rate", xlab = "False Positive Rate", 
           main  = "ROC curve")
      abline(a = 0, b = 1, lty = 2)
      
    })
    
    output$oneLookPlotDuration <- renderPlot({
      
      plot(futilityDF$Assurance, futilityDF$Duration, xlab = "Assurance", xlim = c(0,1), ylab = "Duration",
           main = "Assurance vs Duration for the different stopping rules")
      
    })
    
    output$oneLookPlotSS <- renderPlot({
      
      plot(futilityDF$Assurance, futilityDF$`Sample Size`, xlab = "Assurance", xlim = c(0,1), ylab = "Sample size",
           main = "Assurance vs Sample Size for the different stopping rules")
      
    })
    
    
    
  })
  
  
  # Two Looks Logic ---------------------------------
  
  output$TwoLooksText1 <- renderText({
    TwoLooksSeq1 <- seq(input$TwoLooksLB1, input$TwoLooksUB1, by = input$TwoLooksBy1)
    value <- paste0("We perform IA1 at: ", paste(TwoLooksSeq1, collapse = " "))
    return(value)
  })
  
  output$TwoLooksText2 <- renderText({
    TwoLooksSeq2 <- seq(input$TwoLooksLB2, input$TwoLooksUB2, by = input$TwoLooksBy2)
    value <- paste0("We perform IA2 at: ", paste(TwoLooksSeq2, collapse = " "))
    return(value)
  })
  

  
  observe({
    
   # Get the value from the textbox
      TwoLooksSeq1 <- seq(input$TwoLooksLB1, input$TwoLooksUB1, by = input$TwoLooksBy1)
      TwoLooksSeq2 <- seq(input$TwoLooksLB2, input$TwoLooksUB2, by = input$TwoLooksBy2)

      errorMessage <- NULL

      if (TwoLooksSeq1[1] >= TwoLooksSeq2[1]){
        errorMessage <- "Your first interim analysis needs to be IA1"
      }

      if (TwoLooksSeq1[length(TwoLooksSeq1)] >= TwoLooksSeq2[length(TwoLooksSeq2)]){
        if(is.null(errorMessage)){
          errorMessage <- "Your last interim analysis needs to be IA2"
        } else {
          errorMessage <- "Your first interim analysis needs to be IA1 AND your last interim analysis needs to be IA2"
        }

      }
    
    
      if (!is.null(errorMessage)){
        shinyjs::show("twoLooksErrorMessage")
        output$twoLooksErrorMessage <- renderText({
            errorMessage
        })
        updateActionButton(session, "calcFutilityTwoLooks", disabled = TRUE)
      } else {
        shinyjs::hide("twoLooksErrorMessage")
        if (!is.null(data$treatmentSamplesDF)){
          updateActionButton(session, "calcFutilityTwoLooks", disabled = F)
          
        }
        
      }
    
  })
  

  
  
  observeEvent(input$calcFutilityTwoLooks, {
    NRep <- 5
    TwoLooksSeq1 <- seq(input$TwoLooksLB1, input$TwoLooksUB1, by = input$TwoLooksBy1)
    TwoLooksSeq2 <- seq(input$TwoLooksLB2, input$TwoLooksUB2, by = input$TwoLooksBy2)
    
    TwoLooksCombined <- unique(c(round(TwoLooksSeq1, 2), round(TwoLooksSeq2, 2)))
    

    conc.probs <- matrix(0, 2, 2)
    conc.probs[1, 2] <- 0.5
    
    treatmentSamplesDF <- SHELF::copulaSample(data$treatmentSamplesDF$fit1, data$treatmentSamplesDF$fit2,
                                              cp = conc.probs, n = 1e4, d = data$treatmentSamplesDF$d)
    
    iterationArray <- array(NA,  dim = c(3, length(TwoLooksCombined)+1, NRep))
    
    dimnames(iterationArray) <- list(c("Power", "Duration", "Sample Size"), c(TwoLooksCombined, "Final"), 1:NRep)
    
    proposedDF <- data.frame(matrix(NA, ncol = 9, nrow = NRep))
    
    colnames(proposedDF) <- c("Look1Power", "Look2Power", "FinalLookPower", 
                              "Look1SS", "Look2SS", "FinalLookSS",
                              "Look1Duration", "Look2Duration", "FinalLookDuration")
    

    withProgress(message = 'Calculating', value = 0, {
      for (i in 1:NRep){
        
        #Compute treatment times
        HRStar <- sample(treatmentSamplesDF[,2], 1)
        bigT <- sample(treatmentSamplesDF[,1], 1)
        
        #Simulate control and treatment data
        dataCombined <- SimDTEDataSet(input$numPatients, input$lambdac, bigT, HRStar, input$recTime)  
        
        #Do the proposed rule on this data set
        proposedDF[i,] <- proposedRuleFunc(dataCombined, input$numEvents, 3, 2/3) 
        
        
        #Perform futility look at different Information Fractions
        finalDF <- CensFunc(dataCombined, input$numEvents)
        test <- survdiff(Surv(survival_time, status)~group, data = finalDF$dataCombined)
        coxmodel <- coxph(Surv(survival_time, status)~group, data = finalDF$dataCombined)
        deltad <- as.numeric(exp(coef(coxmodel)))
        
        iterationArray[1, length(TwoLooksCombined)+1, i] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
        iterationArray[2, length(TwoLooksCombined)+1, i] <- finalDF$censTime
        iterationArray[3, length(TwoLooksCombined)+1, i] <- finalDF$SS
        
        
        for (k in 1:length(TwoLooksCombined)){
          futilityCens <- CensFunc(dataCombined, input$numEvents*TwoLooksCombined[k])
          futilityLook <- interimLookFunc(futilityCens$dataCombined, input$OneLookHR)
          
          iterationArray[1, k, i] <- futilityLook
          iterationArray[2, k, i] <- futilityCens$censTime
          iterationArray[3, k, i] <- futilityCens$SS
          
        }
        
        incProgress(1/NRep)
      }
      
    })
      

      PowerArray <- array(NA,  dim = c(length(TwoLooksSeq2), length(TwoLooksSeq1), NRep))
      DurationArray <- array(NA,  dim = c(length(TwoLooksSeq2), length(TwoLooksSeq1), NRep))
      SSArray <- array(NA,  dim = c(length(TwoLooksSeq2), length(TwoLooksSeq1), NRep))
      IATime1 <- array(NA,  dim = c(length(TwoLooksSeq2), length(TwoLooksSeq1), NRep))
      IATime2 <- array(NA,  dim = c(length(TwoLooksSeq2), length(TwoLooksSeq1), NRep))
      PercentStop <- array(NA,  dim = c(length(TwoLooksSeq2), length(TwoLooksSeq1), NRep))
      PercentStopLook1 <- array(NA,  dim = c(length(TwoLooksSeq2), length(TwoLooksSeq1), NRep))
      PercentStopLook2 <- array(NA,  dim = c(length(TwoLooksSeq2), length(TwoLooksSeq1), NRep))
      falselyStopLook1 <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
      correctlyStopLook1 <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
      falselyContinueLook1 <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
      correctlyContinueLook1 <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
      falselyStopLook2 <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
      correctlyStopLook2 <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
      falselyContinueLook2 <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
      correctlyContinueLook2 <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
      falselyStopTotal <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
      correctlyStopTotal <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
      falselyContinueTotal <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
      correctlyContinueTotal <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
      FinalPowerDF <- data.frame(matrix(NA, nrow = NRep, ncol = 3))
      
      #Do proposed rule logic
      proposedDF$power <- proposedDF$Look1Power*proposedDF$Look2Power*proposedDF$FinalLookPower
      proposedDF$SS <- ifelse(proposedDF$Look1Power==0, proposedDF$Look1SS, ifelse(proposedDF$Look2Power==0, 
                                                                                   proposedDF$Look2SS, proposedDF$FinalLookSS))
      proposedDF$Duration <- ifelse(proposedDF$Look1Power==0, proposedDF$Look1Duration, ifelse(proposedDF$Look2Power==0, 
                                                                                               proposedDF$Look2Duration, proposedDF$FinalLookDuration))
      
      for (l in 1:NRep){
        for(i in 1:length(TwoLooksSeq2)) {
          for(j in 1:length(TwoLooksSeq1)) {
            if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
              Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), l]
              Look2Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq2[i])), l]
              FinalPower <- iterationArray[1, length(TwoLooksCombined)+1, l]
              PowerArray[i,j, l] <- Look1Power*Look2Power*FinalPower
              
            
              
              #Sample size and duration logic
              if (Look1Power==0){
                DurationArray[i,j,l] <- iterationArray[2,which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), l]
                SSArray[i,j,l] <- iterationArray[3,which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), l]
              } else if (Look2Power==0){
                DurationArray[i,j,l] <- iterationArray[2,which(colnames(iterationArray) == as.character(TwoLooksSeq2[i])), l]
                SSArray[i,j,l] <- iterationArray[3,which(colnames(iterationArray) == as.character(TwoLooksSeq2[i])), l]
              } else {
                
                DurationArray[i,j,l] <- iterationArray[2,length(TwoLooksCombined)+1, l]
                SSArray[i,j,l] <- iterationArray[3,length(TwoLooksCombined)+1, l]
              }
              
             
              #Interim analysis time logic
              IATime1[i,j,l] <- iterationArray[2, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), l]
              IATime2[i,j,l] <- iterationArray[2, which(colnames(iterationArray) == as.character(TwoLooksSeq2[i])), l]
              
              #Percentage stop logic
              
              PercentStop[i,j,l] <- 0
              PercentStopLook1[i,j,l] <- 0
              PercentStopLook2[i,j,l] <- 0
              
              if (Look1Power==0 || Look2Power==0){
                PercentStop[i,j,l] <- 1
              } 
              
              if (Look1Power==0){
                PercentStopLook1[i,j,l] <- 1
              } else if (Look2Power==0){
                PercentStopLook2[i,j,l] <- 1
              }
              # 
              # if (Look1Power==1){
              #   PercentStopLook1[i,j,l] <- 0
              # } else if (Look2Power==1){
              #   PercentStopLook2[i,j,l] <- 0
              # }
            } 
            
          }
        } 
      }
      
    
      
      #Calculating falsely stopping at Look 1
      for(i in 1:length(TwoLooksSeq2)) {
        for(j in 1:length(TwoLooksSeq1)) {
          if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
            Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), ]
            selected_layers <- iterationArray[,,Look1Power==0]
            if (length(dim(selected_layers))==2){
              falselyStopLook1[i,j] <- mean(selected_layers[1, length(TwoLooksCombined)+1])
            } else {
              if (dim(selected_layers)[3]==0){
                falselyStopLook1[i,j] <- NA
              } else{
                falselyStopLook1[i,j] <- mean(selected_layers[1, length(TwoLooksCombined)+1,])
              }
            }
          } 
        }
      }
      
      #Calculating correctly stopping at Look 1
      for(i in 1:length(TwoLooksSeq2)) {
        for(j in 1:length(TwoLooksSeq1)) {
          if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
            Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), ]
            selected_layers <- iterationArray[,,Look1Power==0]
            if (length(dim(selected_layers))==2){
              correctlyStopLook1[i,j] <- 1- mean(selected_layers[1, length(TwoLooksCombined)+1])
            } else {
              if (dim(selected_layers)[3]==0){
                correctlyStopLook1[i,j] <- NA
              } else{
                correctlyStopLook1[i,j] <- 1- mean(selected_layers[1, length(TwoLooksCombined)+1,])
              }
            }
          } 
        }
      }
      
      #Calculating falsely continuing at Look 1
      for(i in 1:length(TwoLooksSeq2)) {
        for(j in 1:length(TwoLooksSeq1)) {
          if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
            Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), ]
            selected_layers <- iterationArray[,,Look1Power==1]
            if (length(dim(selected_layers))==2){
              falselyContinueLook1[i,j] <- 1- mean(selected_layers[1, length(TwoLooksCombined)+1])
            } else {
              if (dim(selected_layers)[3]==0){
                falselyContinueLook1[i,j] <- NA
              } else{
                falselyContinueLook1[i,j] <- 1- mean(selected_layers[1, length(TwoLooksCombined)+1,])
              }
            }
          } 
        }
      }
      
      #Calculating correctly continuing at Look 1
      for(i in 1:length(TwoLooksSeq2)) {
        for(j in 1:length(TwoLooksSeq1)) {
          if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
            Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), ]
            selected_layers <- iterationArray[,,Look1Power==1]
            if (length(dim(selected_layers))==2){
              correctlyContinueLook1[i,j] <- mean(selected_layers[1, length(TwoLooksCombined)+1])
            } else {
              if (dim(selected_layers)[3]==0){
                correctlyContinueLook1[i,j] <- NA
              } else{
                correctlyContinueLook1[i,j] <- mean(selected_layers[1, length(TwoLooksCombined)+1,])
              }
            }
          } 
        }
      }
      
    
      
      #Calculating falsely stopping at Look 2
      for(i in 1:length(TwoLooksSeq2)) {
        for(j in 1:length(TwoLooksSeq1)) {
          if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
            Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), ]
            Look2Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq2[i])), ]
            selected_layers <- iterationArray[,,Look1Power==1&Look2Power==0]
            if (length(dim(selected_layers))==2){
              falselyStopLook2[i,j] <- mean(selected_layers[1, length(TwoLooksCombined)+1])
            } else {
              if (dim(selected_layers)[3]==0){
                falselyStopLook2[i,j] <- NA
              } else{
                falselyStopLook2[i,j] <- mean(selected_layers[1, length(TwoLooksCombined)+1,])
              }
            }
          } 
        }
      }
      
      #Calculating correctly stopping at Look 2
      for(i in 1:length(TwoLooksSeq2)) {
        for(j in 1:length(TwoLooksSeq1)) {
          if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
            Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), ]
            Look2Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq2[i])), ]
            selected_layers <- iterationArray[,,Look1Power==1&Look2Power==0]
            if (length(dim(selected_layers))==2){
              correctlyStopLook2[i,j] <- 1- mean(selected_layers[1, length(TwoLooksCombined)+1])
            } else {
              if (dim(selected_layers)[3]==0){
                correctlyStopLook2[i,j] <- NA
              } else{
                correctlyStopLook2[i,j] <- 1- mean(selected_layers[1, length(TwoLooksCombined)+1,])
              }
            }
          } 
        }
      }
      
      #Calculating falsely continuing at Look 2
      for(i in 1:length(TwoLooksSeq2)) {
        for(j in 1:length(TwoLooksSeq1)) {
          if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
            Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), ]
            Look2Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq2[i])), ]
            selected_layers <- iterationArray[,,Look1Power==1&Look2Power==1]
            if (length(dim(selected_layers))==2){
              falselyContinueLook2[i,j] <- 1- mean(selected_layers[1, length(TwoLooksCombined)+1])
            } else {
              if (dim(selected_layers)[3]==0){
                falselyContinueLook2[i,j] <- NA
              } else{
                falselyContinueLook2[i,j] <- 1- mean(selected_layers[1, length(TwoLooksCombined)+1,])
              }
            }
          } 
        }
      }
      
      #Calculating correctly continuing at Look 2
      for(i in 1:length(TwoLooksSeq2)) {
        for(j in 1:length(TwoLooksSeq1)) {
          if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
            Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), ]
            Look2Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq2[i])), ]
            selected_layers <- iterationArray[,,Look1Power==1&Look2Power==1]
            if (length(dim(selected_layers))==2){
              correctlyContinueLook2[i,j] <- mean(selected_layers[1, length(TwoLooksCombined)+1])
            } else {
              if (dim(selected_layers)[3]==0){
                correctlyContinueLook2[i,j] <- NA
              } else{
                correctlyContinueLook2[i,j] <- mean(selected_layers[1, length(TwoLooksCombined)+1,])
              }
            }
          } 
        }
      }
      
      #Calculating falsely stopping at all
      for(i in 1:length(TwoLooksSeq2)) {
        for(j in 1:length(TwoLooksSeq1)) {
          if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
            Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), ]
            Look2Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq2[i])), ]
            selected_layers <- iterationArray[,,Look1Power==0|Look2Power==0]
            if (length(dim(selected_layers))==2){
              falselyStopTotal[i,j] <- mean(selected_layers[1, length(TwoLooksCombined)+1])
            } else {
              if (dim(selected_layers)[3]==0){
                falselyStopTotal[i,j] <- NA
              } else{
                falselyStopTotal[i,j] <- mean(selected_layers[1, length(TwoLooksCombined)+1,])
              }
            }
          } 
        }
      }
      
      #Calculating correctly stopping at all
      for(i in 1:length(TwoLooksSeq2)) {
        for(j in 1:length(TwoLooksSeq1)) {
          if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
            Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), ]
            Look2Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq2[i])), ]
            selected_layers <- iterationArray[,,Look1Power==0|Look2Power==0]
            if (length(dim(selected_layers))==2){
              correctlyStopTotal[i,j] <- 1- mean(selected_layers[1, length(TwoLooksCombined)+1])
            } else {
              if (dim(selected_layers)[3]==0){
                correctlyStopTotal[i,j] <- NA
              } else{
                correctlyStopTotal[i,j] <- 1- mean(selected_layers[1, length(TwoLooksCombined)+1,])
              }
            }
          } 
        }
      }
      
      #Calculating falsely continuing at all
      for(i in 1:length(TwoLooksSeq2)) {
        for(j in 1:length(TwoLooksSeq1)) {
          if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
            Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), ]
            Look2Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq2[i])), ]
            selected_layers <- iterationArray[,,Look1Power==1&Look2Power==1]
            if (length(dim(selected_layers))==2){
              falselyContinueTotal[i,j] <- 1- mean(selected_layers[1, length(TwoLooksCombined)+1])
            } else {
              if (dim(selected_layers)[3]==0){
                falselyContinueTotal[i,j] <- NA
              } else{
                falselyContinueTotal[i,j] <- 1- mean(selected_layers[1, length(TwoLooksCombined)+1,])
              }
            }
          } 
        }
      }
      
      #Calculating correctly continuing at Look 2
      for(i in 1:length(TwoLooksSeq2)) {
        for(j in 1:length(TwoLooksSeq1)) {
          if(TwoLooksSeq1[j] < TwoLooksSeq2[i]) {
            Look1Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq1[j])), ]
            Look2Power <- iterationArray[1, which(colnames(iterationArray) == as.character(TwoLooksSeq2[i])), ]
            selected_layers <- iterationArray[,,Look1Power==1&Look2Power==1]
            if (length(dim(selected_layers))==2){
              correctlyContinueTotal[i,j] <- mean(selected_layers[1, length(TwoLooksCombined)+1])
            } else {
              if (dim(selected_layers)[3]==0){
                correctlyContinueTotal[i,j] <- NA
              } else{
                correctlyContinueTotal[i,j] <- mean(selected_layers[1, length(TwoLooksCombined)+1,])
              }
            }
          } 
        }
      }
      
      
      
      
      PowerArray <- apply(PowerArray, c(1, 2), mean, na.rm = TRUE)
      SSArray <- apply(SSArray, c(1, 2), mean, na.rm = TRUE)
      DurationArray <- apply(DurationArray, c(1, 2), mean, na.rm = TRUE)
      
      output$twoLooksAss <- renderTable({
        
        PowerArray <- round(PowerArray, 3)
        
        colnames(PowerArray) <- TwoLooksSeq1
        rownames(PowerArray) <- TwoLooksSeq2
        
        PowerArray
      }, rownames = TRUE)
      
      
      output$twoLooksSS <- renderTable({
        
        
        SSArray <- round(SSArray, 3)
        
        colnames(SSArray) <- TwoLooksSeq1
        rownames(SSArray) <- TwoLooksSeq2
        
        SSArray
      }, rownames = T)
      
      output$twoLooksDuration <- renderTable({
        
        
        
        DurationArray <- round(DurationArray, 3)
        
        colnames(DurationArray) <- TwoLooksSeq1
        rownames(DurationArray) <- TwoLooksSeq2
        
        DurationArray
      }, rownames = T)
      
      output$IATime1Table <- renderTable({
        
        IATime1 <- apply(IATime1, c(1, 2), mean, na.rm = TRUE)
        
        IATime1 <- round(IATime1, 3)
        
        colnames(IATime1) <- TwoLooksSeq1
        rownames(IATime1) <- TwoLooksSeq2
        
        IATime1
      }, rownames = T)
      
      output$IATime2Table <- renderTable({
        
        IATime2 <- apply(IATime2, c(1, 2), mean, na.rm = TRUE)
        
        IATime2 <- round(IATime2, 3)
        
        colnames(IATime2) <- TwoLooksSeq1
        rownames(IATime2) <- TwoLooksSeq2
        
        IATime2
      }, rownames = T)
      
      
      output$PercentStopTable <- renderTable({
        
        PercentStop <- apply(PercentStop, c(1, 2), mean, na.rm = TRUE)
        
        PercentStop <- round(PercentStop, 3)
        
        colnames(PercentStop) <- TwoLooksSeq1
        rownames(PercentStop) <- TwoLooksSeq2
        
        PercentStop
        
        
      }, rownames = T)
      
      output$PercentStopLook1Table <- renderTable({
        
        PercentStopLook1 <- apply(PercentStopLook1, c(1, 2), mean, na.rm = TRUE)
        
        PercentStopLook1 <- round(PercentStopLook1, 3)
        
        colnames(PercentStopLook1) <- TwoLooksSeq1
        rownames(PercentStopLook1) <- TwoLooksSeq2
        
        PercentStopLook1
        
        
      }, rownames = T)
      
      output$PercentStopLook2Table <- renderTable({
        
        PercentStopLook2 <- apply(PercentStopLook2, c(1, 2), mean, na.rm = TRUE)
        
        PercentStopLook2 <- round(PercentStopLook2, 3)
        
        colnames(PercentStopLook2) <- TwoLooksSeq1
        rownames(PercentStopLook2) <- TwoLooksSeq2
        
        PercentStopLook2
        
        
      }, rownames = T)
      
      
      output$falselyStopLook1Table <- renderTable({
        
        colnames(falselyStopLook1) <- TwoLooksSeq1
        rownames(falselyStopLook1) <- TwoLooksSeq2
        
        falselyStopLook1
        
        
      }, rownames = T)
      
      output$correctlyStopLook1Table <- renderTable({
        
        colnames(correctlyStopLook1) <- TwoLooksSeq1
        rownames(correctlyStopLook1) <- TwoLooksSeq2
        
        correctlyStopLook1
        
        
      }, rownames = T)
      
      output$falselyContinueLook1Table <- renderTable({
        
        colnames(falselyContinueLook1) <- TwoLooksSeq1
        rownames(falselyContinueLook1) <- TwoLooksSeq2
        
        falselyContinueLook1
        
        
      }, rownames = T)
      
      output$correctlyContinueLook1Table <- renderTable({
        
        colnames(correctlyContinueLook1) <- TwoLooksSeq1
        rownames(correctlyContinueLook1) <- TwoLooksSeq2
        
        correctlyContinueLook1
        
        
      }, rownames = T)
      
      output$falselyStopLook2Table <- renderTable({
        
        colnames(falselyStopLook2) <- TwoLooksSeq1
        rownames(falselyStopLook2) <- TwoLooksSeq2
        
        falselyStopLook2
        
        
      }, rownames = T)
      
      output$correctlyStopLook2Table <- renderTable({
        
        colnames(correctlyStopLook2) <- TwoLooksSeq1
        rownames(correctlyStopLook2) <- TwoLooksSeq2
        
        correctlyStopLook2
        
        
      }, rownames = T)
      
      output$falselyContinueLook2Table <- renderTable({
        
        colnames(falselyContinueLook2) <- TwoLooksSeq1
        rownames(falselyContinueLook2) <- TwoLooksSeq2
        
        falselyContinueLook2
        
        
      }, rownames = T)
      
      output$correctlyContinueLook2Table <- renderTable({
        
        colnames(correctlyContinueLook2) <- TwoLooksSeq1
        rownames(correctlyContinueLook2) <- TwoLooksSeq2
        
        correctlyContinueLook2
        
        
      }, rownames = T)
      
      output$falselyStopTotalTable <- renderTable({
        
        colnames(falselyStopTotal) <- TwoLooksSeq1
        rownames(falselyStopTotal) <- TwoLooksSeq2
        
        falselyStopTotal
        
        
      }, rownames = T)
      
      output$correctlyStopTotalTable <- renderTable({
        
        colnames(correctlyStopTotal) <- TwoLooksSeq1
        rownames(correctlyStopTotal) <- TwoLooksSeq2
        
        correctlyStopTotal
        
        
      }, rownames = T)
      
      output$falselyContinueTotalTable <- renderTable({
        
        colnames(falselyContinueTotal) <- TwoLooksSeq1
        rownames(falselyContinueTotal) <- TwoLooksSeq2
        
        falselyContinueTotal
        
        
      }, rownames = T)
      
      output$correctlyContinueTotalTable <- renderTable({
        
        colnames(correctlyContinueTotal) <- TwoLooksSeq1
        rownames(correctlyContinueTotal) <- TwoLooksSeq2
        
        correctlyContinueTotal
        
        
      }, rownames = T)
      
      
      
      output$finalAssTable2Looks <- renderTable({
        
        FinalAss <- mean(iterationArray[1, length(TwoLooksCombined)+1, ])
        FinalDuration <- mean(iterationArray[2, length(TwoLooksCombined)+1, ])
        FinalSS <- mean(iterationArray[3, length(TwoLooksCombined)+1, ])
        
        
        FinalAss <- data.frame(Assurance = FinalAss, Duration = FinalDuration, SS = FinalSS)
        
        colnames(FinalAss) <- c("Assurance", "Duration", "Sample Size")
        
        FinalAss
      }, digits = 3)
      
      output$proposedTable2Looks <- renderTable({
        
        
        FinalProposedDF <- data.frame(Assurance = mean(proposedDF$power),
                                      Duration = mean(proposedDF$Duration),
                                      SS = mean(proposedDF$SS))
          
          colnames(FinalProposedDF) <- c("Assurance", "Duration", "Sample Size")
        
          FinalProposedDF
        
        
      }, digits = 3)
      
      output$rocCurveTwoLooks1 <- renderPlot({
        
        spec <- 1 - correctlyStopLook1/(correctlyStopLook1 + falselyStopLook1)
        sens <- correctlyContinueLook1/(correctlyContinueLook1 + falselyContinueLook1)
        
        
        spec[is.na(spec)] <- NaN
        sens[is.na(sens)] <- NaN
        
        spec <- unlist(spec)
        sens <- unlist(sens)
        
        plot(1-spec, sens, ylim = c(0, 1), xlim = c(0,1), ylab = "True Positive Rate", xlab = "False Positive Rate")
        abline(a = 0, b = 1, lty = 2)
        
        

        

      })
        
      
      output$rocCurveTwoLooks2 <- renderPlot({
        
        
        spec <- 1 - correctlyStopLook2/(correctlyStopLook2 + falselyStopLook2)
        sens <- correctlyContinueLook2/(correctlyContinueLook2 + falselyContinueLook2)
        
        spec[is.na(spec)] <- NaN
        sens[is.na(sens)] <- NaN
        
        spec <- unlist(spec)
        sens <- unlist(sens)
        
        plot(1-spec, sens, ylim = c(0, 1), xlim = c(0,1), ylab = "True Positive Rate", xlab = "False Positive Rate")
        abline(a = 0, b = 1, lty = 2)
        
        
      })
      
      output$rocCurveTwoLooksTotal <- renderPlot({
        
        spec <- 1 - correctlyStopTotal/(correctlyStopTotal + falselyStopTotal)
        sens <- correctlyContinueTotal/(correctlyContinueTotal + falselyContinueTotal)
        
        spec[is.na(spec)] <- NaN
        sens[is.na(sens)] <- NaN
        
        spec <- unlist(spec)
        sens <- unlist(sens)
        
        plot(1-spec, sens, ylim = c(0, 1), xlim = c(0,1), ylab = "True Positive Rate", xlab = "False Positive Rate")
        abline(a = 0, b = 1, lty = 2)
        
        
        
      })
        
      
      output$twoLooksPlotDuration <- renderPlot({
        
        plot(PowerArray, DurationArray, xlab = "Assurance", ylab = "Duration", xlim = c(0,1))
        
      })
      
      output$twoLooksPlotSS <- renderPlot({
        
        plot(PowerArray, SSArray,  xlab = "Assurance", ylab = "Sample size", xlim = c(0,1))
        
      })
        
    
    
    output$AssText <- renderUI({
      p(HTML("<b>Assurance table</b>"))
    })
    
    output$SSText <- renderUI({
      p(HTML("<b>Sample size table</b>"))
    })
    
    output$DurationText <- renderUI({
      p(HTML("<b>Duration table</b>"))
    })
    
    output$IATime1Text <- renderUI({
      p(HTML("<b>Interim Analysis 1 Time</b>"))
    })
    
    output$IATime2Text <- renderUI({
      p(HTML("<b>Interim Analysis 2 Time</b>"))
    })
    
    output$PercentStopText <- renderUI({
      p(HTML("<b>% stop (in total)</b>"))
    })
    
    output$PercentStopLook1Text <- renderUI({
      p(HTML("<b>% stop at Look 1</b>"))
    })
    
    output$PercentStopLook2Text <- renderUI({
      p(HTML("<b>% stop at Look 2</b>"))
    })
    
    output$falselyStopLook1Text <- renderUI({
      p(HTML("<b>% falsely stop at Look 1</b>"))
    })
    
    output$correctlyStopLook1Text <- renderUI({
      p(HTML("<b>% correctly stop at Look 1</b>"))
    })
    
    output$falselyContinueLook1Text <- renderUI({
      p(HTML("<b>% falsely continue at Look 1</b>"))
    })
    
    output$correctlyContinueLook1Text <- renderUI({
      p(HTML("<b>% correctly continue at Look 1</b>"))
    })
    
    output$falselyStopLook2Text <- renderUI({
      p(HTML("<b>% falsely stop at Look 2</b>"))
    })
    
    output$correctlyStopLook2Text <- renderUI({
      p(HTML("<b>% correctly stop at Look 2</b>"))
    })
    
    output$falselyContinueLook2Text <- renderUI({
      p(HTML("<b>% falsely continue at Look 2</b>"))
    })
    
    output$correctlyContinueLook2Text <- renderUI({
      p(HTML("<b>% correctly continue at Look 2</b>"))
    })
    
    output$falselyStopTotalText <- renderUI({
      p(HTML("<b>% falsely stop at all</b>"))
    })
    
    output$correctlyStopTotalText <- renderUI({
      p(HTML("<b>% correctly stop at all</b>"))
    })
    
    output$falselyContinueTotalText <- renderUI({
      p(HTML("<b>% falsely continue at all</b>"))
    })
    
    output$correctlyContinueTotalText <- renderUI({
      p(HTML("<b>% correctly continue at all</b>"))
    })
    
    shinyjs::show("AssText")
    shinyjs::show("SSText")
    shinyjs::show("DurationText")
    shinyjs::show("IATime1Text")
    shinyjs::show("IATime2Text")
    shinyjs::show("PercentStopText")
    shinyjs::show("PercentStopLook1Text")
    shinyjs::show("PercentStopLook2Text")
    shinyjs::show("falselyStopLook1Text")
    shinyjs::show("correctlyStopLook1Text")
    shinyjs::show("falselyContinueLook1Text")
    shinyjs::show("correctlyContinueLook1Text")
    shinyjs::show("falselyStopLook2Text")
    shinyjs::show("correctlyStopLook2Text")
    shinyjs::show("falselyContinueLook2Text")
    shinyjs::show("correctlyContinueLook2Text")
    shinyjs::show("falselyStopTotalText")
    shinyjs::show("correctlyStopTotalText")
    shinyjs::show("falselyContinueTotalText")
    shinyjs::show("correctlyContinueTotalText")
    
    
    
    
  })

  
  # Bayesian Logic ---------------------------------
  
  
  #Parallel: # Set up parallel processing
    cl <- makeCluster(detectCores())  # Use all available cores
  registerDoParallel(cl)
  

  observeEvent(input$calcFutilityBayesian, {
    NRep <- 50
    
    conc.probs <- matrix(0, 2, 2)
    conc.probs[1, 2] <- 0.5
    
    # Extract required input values
    numPatients <- input$numPatients
    lambdac <- input$lambdac
    recTime <- input$recTime
    numEvents <- input$numEvents
    IFBayesian <- input$IFBayesian
    
    treatmentSamplesDF <- SHELF::copulaSample(data$treatmentSamplesDF$fit1, data$treatmentSamplesDF$fit2,
                                              cp = conc.probs, n = 1e4, d = data$treatmentSamplesDF$d)
    
    BPPVec <- foreach(i = 1:NRep, .combine = c, .export = c("SimDTEDataSet", "CensFunc", "BPPFunc"),
                      .packages = c("survival", "rjags", "dplyr")) %dopar% {
      
      # Compute treatment times
      HRStar <- sample(treatmentSamplesDF[,2], 1)
      bigT <- sample(treatmentSamplesDF[,1], 1)
      
      # Simulate control and treatment data
      dataCombined <- SimDTEDataSet(numPatients, lambdac, bigT, HRStar, recTime)  
      
      # Perform futility look at different Information Fractions
      finalDF <- CensFunc(dataCombined, numEvents)
      
      test <- survdiff(Surv(survival_time, status) ~ group, data = finalDF$dataCombined)
      coxmodel <- coxph(Surv(survival_time, status) ~ group, data = finalDF$dataCombined)
      deltad <- as.numeric(exp(coef(coxmodel)))
      
      BPPOutcome <- BPPFunc(dataCombined, numPatients, numEvents * IFBayesian, numEvents, recTime)
      
      Success <- (test$chisq > qchisq(0.95, 1) & deltad<1)
      
      return(list(BPP = BPPOutcome$BPP, Success = Success))
    }
    
    stopCluster(cl)  # Stop the cluster
    
    output$BayesianPlot <- renderPlot({
      
      unregister_dopar <- function() {
        env <- foreach:::.foreachGlobals
        rm(list=ls(name=env), pos=env)
      }
      # print(BPPVec)
      # 
      BPPVec <- data.frame(BPP = unlist(BPPVec[seq(1, length(BPPVec), by = 2)]),
                           Success = unlist(BPPVec[seq(2, length(BPPVec), by = 2)]))
      # 
      # 
      # print(BPPVec)                     
      
      
      #plot(BPPVec$BPP, BPPVec$Success)
      
     # x <<- BPPVec
       
      #BPPVec <- data.frame(BPP = BPPVec$BPP, Success = BPPVec$Success)
      
     # print(BPPVec)
      
      # Plotting histogram colored by ColorVar
      ggplot(BPPVec, aes(x = BPP, fill = Success)) +
        geom_histogram(position = "identity", alpha = 0.5) 
        
      
      #hist(BPPVec$B, xlab = "Bayesian Predictive Probability", freq = F)
    })
  })
  
  
  
  
  
}

# Run the Shiny app
shinyApp(ui, server)
