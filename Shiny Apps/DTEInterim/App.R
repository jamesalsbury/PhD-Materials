library(shiny)
library(shinyjs)
library(ggplot2)
library(survival)
library(dplyr)
library(rjags)

source("functions.R")

# UI definition
ui <- fluidPage(
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
                   tableOutput("finalAssTable2Looks"),
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
    if (!is.null(file)) {
      # Read the uploaded file into treatmentSamplesDF
      data$treatmentSamplesDF <-  readRDS(inFile$datapath)

      # Enable the action button when a file is selected
      shinyjs::enable("calcFutilityOneLook")
      shinyjs::enable("calcFutilityTwoLooks")
      shinyjs::enable("calcFutilityBayesian")

       
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
    correctlystop <- rep(NA, length(futilityVec))
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
        correctlystop[k] <- 1-mean(selected_layers[1, length(futilityVec)+1])
      } else {
        if (dim(selected_layers)[3]==0){
          correctlystop[k] <- NA
        } else{
          correctlystop[k] <- 1-mean(selected_layers[1, length(futilityVec)+1, ])
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
                             correctlyStop = correctlystop,
                             falselyContinue = falselyContinue,
                             correctlyContinue = correctlyContinue,
                             Duration = colMeans(durationDF),
                             SS = colMeans(ssDF))
    
    colnames(futilityDF) <- c("Information Fraction", "IA Time", "Assurance", "% stop",
                              "Falsely Stop", "Correctly Stop", "Falsely Continue", "Correctly Continue",
                               "Duration", "Sample Size")
    
    
    output$futilityTable <- renderTable({
      futilityDF
    }, digits = 3)
    
   x <<- futilityDF

    output$finalAssTable <- renderTable({
      
      FinalAss <- mean(iterationArray[1, length(futilityVec)+1, ])
      FinalDuration <- mean(iterationArray[2, length(futilityVec)+1, ])
      FinalSS <- mean(iterationArray[3, length(futilityVec)+1, ])
      
      
      FinalAss <- data.frame(Assurance = FinalAss, Duration = FinalDuration, SS = FinalSS)
      
      colnames(FinalAss) <- c("Assurance", "Duration", "Sample Size")
      
      FinalAss
    }, digits = 3)
    
    output$oneLookPlotDuration <- renderPlot({
      
      plot(futilityDF$Assurance, futilityDF$Duration, xlab = "Assurance", xlim = c(0,1), ylab = "Duration")
      
    })
    
    output$oneLookPlotSS <- renderPlot({
      
      plot(futilityDF$Assurance, futilityDF$`Sample Size`, xlab = "Assurance", xlim = c(0,1), ylab = "Sample size")
      
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
  
  
  observeEvent(input$calcFutilityTwoLooks, {
    NRep <- 500
    TwoLooksSeq1 <- seq(input$TwoLooksLB1, input$TwoLooksUB1, by = input$TwoLooksBy1)
    TwoLooksSeq2 <- seq(input$TwoLooksLB2, input$TwoLooksUB2, by = input$TwoLooksBy2)
    TwoLooksCombined <- unique(c(round(TwoLooksSeq1, 2), round(TwoLooksSeq2, 2)))
    
    conc.probs <- matrix(0, 2, 2)
    conc.probs[1, 2] <- 0.5
    
    treatmentSamplesDF <- SHELF::copulaSample(data$treatmentSamplesDF$fit1, data$treatmentSamplesDF$fit2,
                                              cp = conc.probs, n = 1e4, d = data$treatmentSamplesDF$d)
    
    iterationArray <- array(NA,  dim = c(3, length(TwoLooksCombined)+1, NRep))
    
    dimnames(iterationArray) <- list(c("Power", "Duration", "Sample Size"), c(TwoLooksCombined, "Final"), 1:NRep)
    

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
      
      PowerArray <- array(NA,  dim = c(length(TwoLooksSeq1), length(TwoLooksSeq2), NRep))
      DurationArray <- array(NA,  dim = c(length(TwoLooksSeq1), length(TwoLooksSeq2), NRep))
      SSArray <- array(NA,  dim = c(length(TwoLooksSeq1), length(TwoLooksSeq2), NRep))
      IATime1 <- array(NA,  dim = c(length(TwoLooksSeq1), length(TwoLooksSeq2), NRep))
      IATime2 <- array(NA,  dim = c(length(TwoLooksSeq1), length(TwoLooksSeq2), NRep))
      PercentStop <- array(0,  dim = c(length(TwoLooksSeq1), length(TwoLooksSeq2), NRep))
      PercentStopLook1 <- array(0,  dim = c(length(TwoLooksSeq1), length(TwoLooksSeq2), NRep))
      PercentStopLook2 <- array(0,  dim = c(length(TwoLooksSeq1), length(TwoLooksSeq2), NRep))
      FinalPowerDF <- data.frame(matrix(NA, nrow = NRep, ncol = 3))
      
      

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
              if (Look1Power==0 || Look2Power==0){
                PercentStop[i,j,l] <- 1
              }
              
              if (Look1Power==0){
                PercentStopLook1[i,j,l] <- 1
              } else if (Look2Power==0){
                PercentStopLook2[i,j,l] <- 1
              }
            } 
            
          }
        } 
      }
      
      
      PowerArray <- apply(PowerArray, c(1, 2), mean, na.rm = TRUE)
      SSArray <- apply(SSArray, c(1, 2), mean, na.rm = TRUE)
      DurationArray <- apply(DurationArray, c(1, 2), mean, na.rm = TRUE)
      
      x <<- PowerArray
      y <<- DurationArray
      
      
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
      
      
      
      
      output$finalAssTable2Looks <- renderTable({
        
        FinalAss <- mean(iterationArray[1, length(TwoLooksCombined)+1, ])
        FinalDuration <- mean(iterationArray[2, length(TwoLooksCombined)+1, ])
        FinalSS <- mean(iterationArray[3, length(TwoLooksCombined)+1, ])
        
        
        FinalAss <- data.frame(Assurance = FinalAss, Duration = FinalDuration, SS = FinalSS)
        
        colnames(FinalAss) <- c("Assurance", "Duration", "Sample Size")
        
        FinalAss
      }, digits = 3)
      
      output$twoLooksPlotDuration <- renderPlot({
        
        plot(PowerArray, DurationArray, xlab = "Assurance", ylab = "Duration", xlim = c(0,1))
        
      })
      
      output$twoLooksPlotSS <- renderPlot({
        
        plot(PowerArray, SSArray,  xlab = "Assurance", ylab = "Sample size", xlim = c(0,1))
        
      })
        
      
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
  
    
    shinyjs::show("AssText")
    shinyjs::show("SSText")
    shinyjs::show("DurationText")
    shinyjs::show("IATime1Text")
    shinyjs::show("IATime2Text")
    shinyjs::show("PercentStopText")
    shinyjs::show("PercentStopLook1Text")
    shinyjs::show("PercentStopLook2Text")
    
    
    
    
  })

  
  # Bayesian Logic ---------------------------------
  
  
  observeEvent(input$calcFutilityBayesian, {
    NRep <- 10
    
    conc.probs <- matrix(0, 2, 2)
    conc.probs[1, 2] <- 0.5
    
    treatmentSamplesDF <- SHELF::copulaSample(data$treatmentSamplesDF$fit1, data$treatmentSamplesDF$fit2,
                                              cp = conc.probs, n = 1e4, d = data$treatmentSamplesDF$d)
    
    BPPVec <- rep(NA, NRep)
    
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
        
        #Censor at chosen IF
        #censoredDF <- CensFunc(dataCombined, input$numEvents*input$IFBayesian)
        
       # print(censoredDF)
        
        # iterationArray[1, length(TwoLooksCombined)+1, i] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
        # iterationArray[2, length(TwoLooksCombined)+1, i] <- finalDF$censTime
        # iterationArray[3, length(TwoLooksCombined)+1, i] <- finalDF$SS
        
        #y <<- dataCombined
        
        
        BPPOutcome <- BPPFunc(dataCombined, input$numPatients, input$numEvents*input$IFBayesian, input$numEvents, input$recTime)
        
        BPPVec[i] <- BPPOutcome$BPP
          
        #print(x)
        
        incProgress(1/NRep)
      }
      
    })
    
    
    
    
    output$BayesianPlot <- renderPlot({
      
      BPPDF <- data.frame(BPPVec = BPPVec)
      print("yes")
      print(BPPVec)
      
      p <- ggplot(BPPDF, aes(x = BPPVec)) +
        geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black", aes(y = after_stat(density))) +
        geom_density(alpha = 0.2, fill = "orange") +
        labs(title = "Histogram of Sample Data", x = "Values", y = "Density") +
        theme_minimal()
      
      p
      
    })
    
    
    
    
  })
  
  
}

# Run the Shiny app
shinyApp(ui, server)
