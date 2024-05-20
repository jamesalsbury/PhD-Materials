library(shiny)
library(shinyjs)
library(survival)
library(doParallel)
library(ggplot2)
library(DT)
library(rhandsontable)
library(rpact)
library(SHELF)
library(plotly)
# library(dplyr)
# library(rjags)
# library(magrittr)
# library(future.apply)
#library(foreach)
# library(future)
# library(promises)
# library(ipc)


source("functions.R")

rowCallback <- c(
  "function(row, data){",
  "  for(var i=0; i<data.length; i++){",
  "    if(data[i] === null){",
  "      $('td:eq('+i+')', row).html('NA')",
  "        .css({'color': 'rgb(151,151,151)', 'font-style': 'italic'});",
  "    }",
  "  }",
  "}"
)


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
  titlePanel("Bayesian Interim Analyses: Delayed Treatment Effects"),
  
  mainPanel(
    tabsetPanel(
      # Design UI ---------------------------------
      tabPanel("Design", 
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
      
      # No Look UI ---------------------------------
      
      tabPanel("No Look", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   actionButton("calcNoLook", label  = "Calculate", disabled = T)
                 ), 
                 mainPanel = mainPanel(
                   hidden(numericInput("noLookAssuranceValue", "Number of patients (altogether)", value = 10)),
                   textOutput("noLookAssuranceText"),
                   plotOutput("noLookAssurancePlot"),
                   tableOutput("finalAssTableNoLook"),
                   
                 )
               )
      ),
      
      # One Look UI ---------------------------------
      
      
      tabPanel("One Look", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   wellPanel(
                     fluidRow(
                       column(4, numericInput("OneLookLB", "IA1, from:", value = 0.25)),
                       column(4, numericInput("OneLookUB", "to:", value = 0.75)),
                       column(4, numericInput("OneLookBy", "by:", value = 0.25))
                     ),
                     textOutput("OneLookText"),
                     div(id = "oneLookErrorMessage", class = "error-message", textOutput("oneLookErrorMessage"))),
                   wellPanel(
                     rHandsontableOutput("spendingOneLook")),
                   actionButton("calcOneLook", label  = "Calculate", disabled = T)
                 ), 
                 mainPanel = mainPanel(
                   tabsetPanel(
                     tabPanel("Tables",
                              hidden(selectizeInput("selectedOptionsIATableOneLook", "Selected Metrics", 
                                                    choices = c("Interim Analysis Time", "Assurance", "Duration", "Sample Size",
                                                                "% Stop", "% Stop for Efficacy", "% Stop for Futility",
                                                                "% Correctly Stop", "% Correctly Stop for Efficacy", "% Correctly Stop for Futility",
                                                                "% Correctly Continue"), 
                                                    selected = c("Assurance", "Duration", "Sample Size"),
                                                    multiple = TRUE)), 
                              DTOutput("IATableOneLook"),
                              hidden(uiOutput("finalAssTable1LookText")),
                              tableOutput("noIATableOneLook")),
                     tabPanel("Plots",
                              selectInput("oneLookBoundaryIA", "Choose the IF (to view)", choices = NULL),
                              plotlyOutput("oneLookBoundaries"),
                              plotlyOutput("oneLookPlotDuration"),
                              plotlyOutput("oneLookPlotSS")),

                   ),
                   
                 )
               )
      ),
      
      # Two Looks UI ---------------------------------
      
      tabPanel("Two Looks", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   wellPanel(
                     fluidRow(
                       column(4, numericInput("TwoLooksLB1", "IA1, from:", value = 0.2)),
                       column(4, numericInput("TwoLooksUB1", "to:", value = 0.8)),
                       column(4, numericInput("TwoLooksBy1", "by:", value = 0.2))
                     ),
                     fluidRow(
                       column(4, numericInput("TwoLooksLB2", "IA2, from:", value = 0.3)),
                       column(4, numericInput("TwoLooksUB2", "to:", value = 0.9)),
                       column(4, numericInput("TwoLooksBy2", "by:", value = 0.2))
                     ),
                     textOutput("TwoLooksText1"),
                     textOutput("TwoLooksText2"),
                     div(id = "twoLooksErrorMessage", class = "error-message", textOutput("twoLooksErrorMessage"))),
                   #selectInput("sidedTwoLooks", "Test", choices = c("One-sided", "Two-sided"), selected = "One-sided"),
                   wellPanel(
                     rHandsontableOutput("spendingTwoLooks")),
                   actionButton("calcTwoLooks", label  = "Calculate", disabled = T)
                 ), 
                 mainPanel = mainPanel(
                   tabsetPanel(
                     tabPanel("Tables",
                              hidden(selectizeInput("selectedOptionstwoLooks", "Selected Metrics", 
                                                    choices = NULL, 
                                                    multiple = TRUE)),
                              uiOutput("metricsTableoutput")),
                     tabPanel("Plots",
                              selectInput("twoLooksBoundaryIA", "Choose the IFs (to view)", choices = NULL),
                              plotlyOutput("twoLooksBoundaries"),
                              plotlyOutput("twoLooksPlotDuration"),
                              plotlyOutput("twoLooksPlotSS"))
                   ),

                 )
               )
      ),
      
      
      # Bayesian UI ---------------------------------
      
      tabPanel("Bayesian", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   numericInput("IFBayesian", "Information Fraction", value = 0.5),
                   numericInput("tEffBayesian", "Target Effect", value = 0.8),
                   div(id = "bayesianErrorMessage", class = "error-message", textOutput("bayesianErrorMessage")),
                   actionButton("calcBayesian", label  = "Calculate", disabled = T)
                 ), 
                 mainPanel = mainPanel(
                   plotOutput("BayesianPlot"),
                   plotOutput("BayesianEffPlot"),
                   plotOutput("BayesianBPPvTE"),
                   tableOutput("BayesianSS"),
                   tableOutput("BayesianDuration")
                   
                 )
               )
      ),
      
      # Report UI ---------------------------------
      
      tabPanel("Report", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   fluidRow(
                     column(6, checkboxInput("checkDesign", "Design", value = FALSE)),
                     column(6, selectizeInput("checkDesignOptions", "Selected Metrics", 
                                                     choices = c("Interim Analysis Time", "Assurance", "Duration", "Sample Size",
                                                                 "% Stop", "% Stop for Efficacy", "% Stop for Futility",
                                                                 "% Correctly Stop", "% Correctly Stop for Efficacy", "% Correctly Stop for Futility",
                                                                 "% Correctly Continue"), 
                                                     selected = c("Assurance", "Duration", "Sample Size"),
                                                     multiple = TRUE))
                   ),
                   fluidRow(
                     column(6, checkboxInput("checkNoLook", "No Look", value = FALSE)),
                     column(6, selectizeInput("checkNoLookOptions", "Selected Metrics", 
                                              choices = c("Interim Analysis Time", "Assurance", "Duration", "Sample Size",
                                                          "% Stop", "% Stop for Efficacy", "% Stop for Futility",
                                                          "% Correctly Stop", "% Correctly Stop for Efficacy", "% Correctly Stop for Futility",
                                                          "% Correctly Continue"), 
                                              selected = c("Assurance", "Duration", "Sample Size"),
                                              multiple = TRUE))
                   ),
                   
                   fluidRow(
                     column(6, checkboxInput("checkOneLook", "One Look", value = FALSE)),
                     column(6, selectizeInput("checkOneLookOptions", "Selected Metrics", 
                                              choices = c("Interim Analysis Time", "Assurance", "Duration", "Sample Size",
                                                          "% Stop", "% Stop for Efficacy", "% Stop for Futility",
                                                          "% Correctly Stop", "% Correctly Stop for Efficacy", "% Correctly Stop for Futility",
                                                          "% Correctly Continue"), 
                                              selected = c("Assurance", "Duration", "Sample Size"),
                                              multiple = TRUE))
                   ),
                   fluidRow(
                     column(6, checkboxInput("checkTwoLooks", "Two Looks", value = FALSE)),
                     column(6, selectizeInput("checkTwoLooksOptions", "Selected Metrics", 
                                              choices = c("Interim Analysis Time", "Assurance", "Duration", "Sample Size",
                                                          "% Stop", "% Stop for Efficacy", "% Stop for Futility",
                                                          "% Correctly Stop", "% Correctly Stop for Efficacy", "% Correctly Stop for Futility",
                                                          "% Correctly Continue"), 
                                              selected = c("Assurance", "Duration", "Sample Size"),
                                              multiple = TRUE))
                   ),
                   fluidRow(
                     column(6, checkboxInput("checkBayesian", "Bayesian", value = FALSE)),
                     column(6, selectizeInput("checkBayesianOptions", "Selected Metrics", 
                                              choices = c("Interim Analysis Time", "Assurance", "Duration", "Sample Size",
                                                          "% Stop", "% Stop for Efficacy", "% Stop for Futility",
                                                          "% Correctly Stop", "% Correctly Stop for Efficacy", "% Correctly Stop for Futility",
                                                          "% Correctly Continue"), 
                                              selected = c("Assurance", "Duration", "Sample Size"),
                                              multiple = TRUE))
                   )
                   
                   
                 ), 
                 mainPanel = mainPanel(
                   downloadButton("downloadHTML", "Html"),
                   downloadButton("downloadPDF", "Pdf")
                   
                 )
               )
      )
      
      
      
    ), style='width: 1200px; height: 1000px'
  )
)


# Server logic
server <- function(input, output, session) {
  
  # Design Logic ---------------------------------
  
  reactValues <- reactiveValues(treatmentSamplesDF = NULL,
                                errorSeqOneLook = FALSE,
                                errorSeqTwoLooks = FALSE,
                                errorBayesian = FALSE,
                                iterationList = NULL,
                                TwoLooksSeq1 = NULL,
                                TwoLooksSeq2 = NULL)
  
  observeEvent(input$samplesFile, {
    inFile  <- input$samplesFile
    if (!is.null(inFile)) {
      # Read the uploaded file into treatmentSamplesDF
      reactValues$treatmentSamplesDF <-  readRDS(inFile$datapath)
      
      
      output$TDist <- renderPlot({
        SHELF::plotfit(reactValues$treatmentSamplesDF$fit1, d = reactValues$treatmentSamplesDF$d[1])
      })
      
      output$HRDist <- renderPlot({
        SHELF::plotfit(reactValues$treatmentSamplesDF$fit2, d = reactValues$treatmentSamplesDF$d[2])
      })
      
      output$P_S <- renderText({
        paste0("P_S is: ", reactValues$treatmentSamplesDF$P_S)
      })
      
      output$P_DTE <- renderText({
        paste0("P_DTE is: ", reactValues$treatmentSamplesDF$P_DTE)
      })
    }
  })
  
  
  # No Look Logic ---------------------------------
  
  observe({
    if (is.null(reactValues$treatmentSamplesDF)) {
      updateActionButton(session, "calcNoLook", disabled = TRUE)
    } else {
      updateActionButton(session, "calcNoLook", disabled = FALSE)
    }
  })
  
  observeEvent(input$calcNoLook, {
    
    NRep <- 500
    
    conc.probs <- matrix(0, 2, 2)
    conc.probs[1, 2] <- 0.5
    
    iterationList <- vector("list", length = 10)
    SSVec <- ceiling(seq(30, 1001, length = length(iterationList)))
    nEventsVec <- ceiling((input$numEvents/input$numPatients)*SSVec)
    
    for (i in 1:length(iterationList)){
      iterationList[[i]]$SampleSize <- SSVec[i]
      iterationList[[i]]$nEvents <- nEventsVec[i]
    }
    
    withProgress(message = 'Calculating', value = 0, {
    for (i in 1:length(iterationList)){
      treatmentSamplesDF <- SHELF::copulaSample(reactValues$treatmentSamplesDF$fit1, reactValues$treatmentSamplesDF$fit2,
                                                     cp = conc.probs, n = 1e4, d = reactValues$treatmentSamplesDF$d)
      
      for (j in 1:NRep){
        
          #Compute treatment times
          HRStar <- sample(treatmentSamplesDF[,2], 1)
          bigT <- sample(treatmentSamplesDF[,1], 1)

          #Simulate control and treatment data
          dataCombined <- SimDTEDataSet(iterationList[[i]]$SampleSize, input$lambdac, bigT, HRStar, input$recTime)

          #Perform looks at different Information Fractions
          finalDF <- CensFunc(dataCombined, iterationList[[i]]$nEvents)
          test <- survdiff(Surv(survival_time, status)~group, data = finalDF$dataCombined)
          coxmodel <- coxph(Surv(survival_time, status)~group, data = finalDF$dataCombined)
          deltad <- as.numeric(exp(coef(coxmodel)))
          
          iterationList[[i]]$PowerVec[j] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
          iterationList[[i]]$DurationVec[j] <- finalDF$censTime
          iterationList[[i]]$SSVec[j] <- finalDF$SS
        
      }
      incProgress(1/length(iterationList))
    }
    })
    
    
    for (i in 1:length(iterationList)){
      iterationList[[i]]$Power <- mean(iterationList[[i]]$PowerVec)
      iterationList[[i]]$Duration <- mean(iterationList[[i]]$DurationVec)
    }
    
    myDF <- data.frame(SampleSize = rep(NA, length(iterationList)), Power = rep(NA, length(iterationList)),
                       Duration = rep(NA, length(iterationList))) 
    
    for (i in 1:length(iterationList)){
      myDF[i,]$Power <- iterationList[[i]]$Power
      myDF[i,]$SampleSize <- iterationList[[i]]$SampleSize
      myDF[i,]$Duration <- iterationList[[i]]$Duration
    }
    
    smoothedPower <- loess(Power ~ SampleSize, data = myDF)
    smoothedDuration <- loess(Duration ~ SampleSize, data = myDF)
    
    updateNumericInput(session, "noLookAssuranceValue", value = input$numPatients*2)
    shinyjs::show("noLookAssuranceValue")
    
    output$noLookAssurancePlot <- renderPlot({

      plot(myDF$SampleSize*2, predict(smoothedPower), ylim = c(0,1), type = "l", xlab = "Sample Size", ylab = "Assurance")
      abline(v = input$noLookAssuranceValue, lty = 2)
      
    })
    
    
    output$noLookAssuranceText <- renderText({
      
      estAss <- predict(smoothedPower, newdata = input$noLookAssuranceValue/2)
      
      paste0("With ", input$noLookAssuranceValue, " patients altogether, the assurance is estimated to be: ", 
             round(estAss, 2))
             
    })
    
  #     }
  #   })
  #   
  #   
  #   output$finalAssTableNoLook <- renderTable({
  #     
  #     FinalAss <- mean(iterationArray[1, 1, ])
  #     FinalDuration <- mean(iterationArray[2, 1, ])
  #     FinalSS <- mean(iterationArray[3, 1, ])
  #     
  #     
  #     FinalAss <- data.frame(Assurance = FinalAss, Duration = FinalDuration, SS = FinalSS)
  #     
  #     colnames(FinalAss) <- c("Assurance", "Duration", "Sample Size")
  #     
  #     FinalAss
  #   }, digits = 3)
  #   
  })
  
  
  
  # One Look Logic ---------------------------------
  
  
  observe({
    if (!is.null(reactValues$treatmentSamplesDF)&reactValues$errorSeqOneLook==F) {
      updateActionButton(session, "calcOneLook", disabled = FALSE)
    } else {
      updateActionButton(session, "calcOneLook", disabled = TRUE)
    }
  })
  
  
  observe({
    # Check if any of the inputs are NA
    if (is.na(input$OneLookLB) ||
        is.na(input$OneLookUB) ||
        is.na(input$OneLookBy))  {
      reactValues$errorSeqOneLook <- TRUE
    } else {
      # Check validity of input ranges
      if (input$OneLookLB <= 0 | input$OneLookLB >= 1 ||
          input$OneLookUB <= 0 | input$OneLookUB >= 1 ||
          input$OneLookBy <= 0 | input$OneLookBy >= 1 ||
          input$OneLookLB > input$OneLookUB) {
        reactValues$errorSeqOneLook <- TRUE
      } else {
        reactValues$errorSeqOneLook <- FALSE
      }
    }
  })
  
  observe({
    if (reactValues$errorSeqOneLook==TRUE){
      shinyjs::show("oneLookErrorMessage")
    } else{
      shinyjs::hide("oneLookErrorMessage")
    }
  })
  
  output$oneLookErrorMessage <- renderText({
    if (reactValues$errorSeqOneLook==TRUE){
      return("Your inputs are incorrect!")
    } 
    
    return("")
    
  })
  
  
  
  output$OneLookText <- renderText({
    if (reactValues$errorSeqOneLook == TRUE) {
      return("")
    }
    
    OneLookSeq <- seq(input$OneLookLB, input$OneLookUB, by = input$OneLookBy)
    value <- paste0("We perform IA1 at: ", paste(OneLookSeq, collapse = ", "))
    return(value)
  })
  
  
  output$spendingOneLook <- renderRHandsontable({
    initial_data <- data.frame(
      Stage = 1:2,
      alphaspending = c("0.0125", "0.025"),
      betaspending = c("0.05", "0.1")
    )
    
    colnames(initial_data) <- c("Stage", "Alpha spending", "Beta spending")
    
    rhandsontable(initial_data, rowHeaders = FALSE)  %>%
      hot_col(col = "Stage", readOnly = TRUE)
  })
  
  observe({
    if (reactValues$errorSeqOneLook==F){
      updateSelectInput(session, "oneLookBoundaryIA", choices = seq(input$OneLookLB, input$OneLookUB, by = input$OneLookBy))
    }
  })
  
  observeEvent(c(input$spendingOneLook, input$oneLookBoundaryIA), {
    
    values <- hot_to_r(input$spendingOneLook)
    
    if (!is.null(values)){
      design <- getDesignGroupSequential(typeOfDesign = "asUser", 
                                         informationRates = c(as.numeric(input$oneLookBoundaryIA), 1),
                                         userAlphaSpending = as.numeric(values[,2]), 
                                         typeBetaSpending = "bsUser",
                                         userBetaSpending = as.numeric(values[,3]))
      
      output$oneLookBoundaries <- renderPlotly({
        plot(design)
      })
    } 
    
    
  })
  
  
  observeEvent(input$calcOneLook, {
    NRep <- 20
    
    IAVec <- seq(input$OneLookLB, input$OneLookUB, by = input$OneLookBy)
    
    iterationList <- vector("list", length = length(IAVec))
    
    noLooksDF <- data.frame(Power = numeric(NRep),
                            Duration = numeric(NRep),
                            SS = numeric(NRep))
    
    values <- hot_to_r(input$spendingOneLook)
    
    for (k in 1:length(IAVec)){
      
      iterationList[[k]]$IF <- c(IAVec[k], 1) 
      
      design <- getDesignGroupSequential(typeOfDesign = "asUser", 
                                         informationRates = iterationList[[k]]$IF,
                                         userAlphaSpending = as.numeric(values[,2]), 
                                         typeBetaSpending = "bsUser",
                                         userBetaSpending = as.numeric(values[,3]))
      
      iterationList[[k]]$EffBounds <- design$criticalValues
      iterationList[[k]]$futBounds <- design$futilityBounds
      
      
    }
    
    conc.probs <- matrix(0, 2, 2)
    conc.probs[1, 2] <- 0.5
    
    treatmentSamplesDF <- SHELF::copulaSample(reactValues$treatmentSamplesDF$fit1, reactValues$treatmentSamplesDF$fit2,
                                              cp = conc.probs, n = 1e4, d = reactValues$treatmentSamplesDF$d)
    
    withProgress(message = 'Calculating', value = 0, {
      for (i in 1:NRep){
        
        #Compute treatment times
        HRStar <- sample(treatmentSamplesDF[,2], 1)
        bigT <- sample(treatmentSamplesDF[,1], 1)
        
        #Simulate control and treatment data
        dataCombined <- SimDTEDataSet(input$numPatients, input$lambdac, bigT, HRStar, input$recTime)  
        
        #Perform looks at different Information Fractions
        finalDF <- CensFunc(dataCombined, input$numEvents)
        test <- survdiff(Surv(survival_time, status)~group, data = finalDF$dataCombined)
        coxmodel <- coxph(Surv(survival_time, status)~group, data = finalDF$dataCombined)
        deltad <- as.numeric(exp(coef(coxmodel)))
        
        noLooksDF$Power[i] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
        noLooksDF$Duration[i] <- finalDF$censTime
        noLooksDF$SS[i] <- finalDF$SS
        
        
        for (k in 1:length(IAVec)){
          
          GSDOC <- GSDOneIAFunc(dataCombined, iterationList[[k]]$futBounds, 
                                iterationList[[k]]$EffBounds, iterationList[[k]]$IF, input$numEvents) 
          
          iterationList[[k]]$Power[i] <- ifelse(GSDOC$Outcome %in% c("Efficacy", "Successful"), 1, 0)
          iterationList[[k]]$Duration[i] <- GSDOC$Duration
          iterationList[[k]]$SS[i] <- GSDOC$SS
          iterationList[[k]]$IA1Time[i] <- GSDOC$IA1Time
          iterationList[[k]]$Outcome[i] <- GSDOC$Outcome
          iterationList[[k]]$Stop[i] <- ifelse(GSDOC$Outcome %in% c("Efficacy", "Futility"), 1, 0)
          iterationList[[k]]$StopEff[i] <- ifelse(GSDOC$Outcome == "Efficacy", 1, 0)
          iterationList[[k]]$StopFut[i] <- ifelse(GSDOC$Outcome == "Futility", 1, 0)
          iterationList[[k]]$Truth[i] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
          
          
        }
        
        incProgress(1/NRep)
      }
    })
    
    
    correctlyStop <- rep(NA, length(IAVec))
    correctlyStopEff <- rep(NA, length(IAVec))
    correctlyStopFut <- rep(NA, length(IAVec))
    correctlyContinue <- rep(NA, length(IAVec))
    
    for (k in 1:length(IAVec)) {
      
      # Calculating Correctly Stop
      IndStopVec <- which(iterationList[[k]]$Stop == 1)
      if (length(IndStopVec) != 0) {
        OutcomeVec <- iterationList[[k]]$Outcome[IndStopVec]
        TruthVec <- iterationList[[k]]$Truth[IndStopVec]
        correctlyStopSum <- 0
        
        for (i in 1:length(IndStopVec)) {
          if ((OutcomeVec[i] == "Efficacy" & TruthVec[i] == 1) | (OutcomeVec[i] == "Futility" & TruthVec[i] == 0)) {
            correctlyStopSum <- correctlyStopSum + 1
          }
        }
        
        correctlyStop[k] <- correctlyStopSum / length(IndStopVec)
      } else {
        correctlyStop[k] <- NA
      }
      
      # Calculating Correctly Stopping for Efficacy
      IndStopEffVec <- which(iterationList[[k]]$StopEff == 1)
      if (length(IndStopEffVec) != 0) {
        TruthVec <- iterationList[[k]]$Truth[IndStopEffVec]
        correctlyStopEffSum <- 0
        
        for (i in 1:length(IndStopEffVec)) {
          if (TruthVec[i] == 1) {
            correctlyStopEffSum <- correctlyStopEffSum + 1
          }
        }
        
        correctlyStopEff[k] <- correctlyStopEffSum / length(IndStopEffVec)
      } else {
        correctlyStopEff[k] <- NA
      }
      
      # Calculating Correctly Stopping for Futility
      IndStopFutVec <- which(iterationList[[k]]$StopFut == 1)
      if (length(IndStopFutVec) != 0) {
        TruthVec <- iterationList[[k]]$Truth[IndStopFutVec]
        correctlyStopFutSum <- 0
        
        for (i in 1:length(IndStopFutVec)) {
          if (TruthVec[i] == 0) {
            correctlyStopFutSum <- correctlyStopFutSum + 1
          }
        }
        
        correctlyStopFut[k] <- correctlyStopFutSum / length(IndStopFutVec)
      } else {
        correctlyStopFut[k] <- NA
      }
      
      # Calculating Correctly Continuing
      IndContinueCorrectlyVec <- which(iterationList[[k]]$Stop == 0)
      if (length(IndContinueCorrectlyVec) != 0) {
        TruthVec <- iterationList[[k]]$Truth[IndContinueCorrectlyVec]
        correctlyContinueSum <- 0
        
        for (i in 1:length(IndContinueCorrectlyVec)) {
          if (TruthVec[i] == 1) {
            correctlyContinueSum <- correctlyContinueSum + 1
          }
        }
        
        correctlyContinue[k] <- correctlyContinueSum / length(IndContinueCorrectlyVec)
      } else {
        correctlyContinue[k] <- NA
      }
      
    }
    
    
    
    IADFOneLook <- data.frame(IF = IAVec,
                              IATime = unlist(lapply(iterationList, function(x) mean(x$IA1Time))),
                              Assurance = unlist(lapply(iterationList, function(x) mean(x$Power))),
                              Duration = unlist(lapply(iterationList, function(x) mean(x$Duration))),
                              SS = unlist(lapply(iterationList, function(x) mean(x$SS))),
                              Stop = unlist(lapply(iterationList, function(x) mean(x$Stop))),
                              StopEff = unlist(lapply(iterationList, function(x) mean(x$StopEff))),
                              StopFut = unlist(lapply(iterationList, function(x) mean(x$StopFut))),
                              correctlyStop = correctlyStop,
                              correctlyStopEff = correctlyStopEff,
                              correctlyStopFut = correctlyStopFut,
                              correctlyContinue = correctlyContinue)
    
    
    
    shinyjs::show("selectedOptionsIATableOneLook")
    
    
    #Continuing is positive!
    
    colnames(IADFOneLook) <- c("Information Fraction", "Interim Analysis Time", "Assurance", "Duration", "Sample Size",
                               "% Stop", "% Stop for Efficacy", "% Stop for Futility",
                               "% Correctly Stop", "% Correctly Stop for Efficacy", "% Correctly Stop for Futility",
                               "% Correctly Continue")
    
    FinalAss <- t(colMeans(noLooksDF))
    FinalAss <- as.data.frame(FinalAss)
    colnames(FinalAss) <- c("Assurance", "Duration", "Sample Size")
    
    #Here
    
    output$IATableOneLook <- renderDT({
      
      IADFOneLook <- subset(IADFOneLook, select = c("Information Fraction", input$selectedOptionsIATableOneLook))
      
      datatable(IADFOneLook, options = list(rowCallback = JS(rowCallback)),
                rownames = F) %>% formatStyle(
                  columns = colnames(IADFOneLook)
                  ) %>%
        formatSignif(
          columns = colnames(IADFOneLook),
          digits = 3
        )
    })
    
    
    output$noIATableOneLook <- renderTable({
      FinalAss
    }, digits = 3)
    
    
    output$oneLookPlotDuration <- renderPlotly({
      
      p <- plot_ly(IADFOneLook, x = ~Assurance, y = ~Duration,
                   text = ~ paste0("Information Fraction = ", `Information Fraction`), mode = "markers",
                   type = "scatter", marker = list(size = 10, color = "blue"), name = "Chosen Rules") %>%
        add_trace(x = ~FinalAss$Assurance, y = ~FinalAss$Duration, type = "scatter", mode = "markers", 
                  marker = list(size = 10, color = "red"), 
                  text = ~ "No Interim Analysis",
                  name = "No Interim Analysis") %>%
        layout(
          xaxis = list(title = list(text = "Assurance", font = list(color = "black", size = 14, family = "Arial", weight = "bold")), range = c(0, 1)),
          yaxis = list(title = list(text = "Duration", font = list(color = "black", size = 14, family = "Arial", weight = "bold"))),
          legend = list(orientation = "v", x = 1.05, y = 0.5),  # Position legend to the right
          title = "Assurance vs Duration for the different stopping rules"
        )
      
      p
      
      
      # plot(IADFOneLook$Assurance, IADFOneLook$Duration, xlab = "Assurance", xlim = c(0,1), ylab = "Duration",
      #      ylim = c(min(IADFOneLook$Duration), FinalAss$Duration),
      #      main = "Assurance vs Duration for the different stopping rules", pch = 19)
      # text(IADFOneLook$Assurance, IADFOneLook$Duration, IADFOneLook$`Information Fraction`, pos = 2)  
      # points(FinalAss$Assurance, FinalAss$Duration, col = "red", pch = 19)
      # text(FinalAss$Assurance, FinalAss$Duration, "No IA", col = "red", pos = 4)
      # 
      # legend("topleft", legend = c("Chosen Rules", "No IA"), col = c("black", "red"), pch = 19)
      
    })
    
    output$oneLookPlotSS <- renderPlotly({
      
      p <- plot_ly(IADFOneLook, x = ~Assurance, y = ~`Sample Size`,
                   text = ~ paste0("Information Fraction = ", `Information Fraction`), mode = "markers",
                   type = "scatter", marker = list(size = 10, color = "blue"), name = "Chosen Rules") %>%
        add_trace(x = ~FinalAss$Assurance, y = ~FinalAss$`Sample Size`, type = "scatter", mode = "markers", 
                  marker = list(size = 10, color = "red"), 
                  text = ~ "No Interim Analysis",
                  name = "No Interim Analysis") %>%
        layout(
          xaxis = list(title = list(text = "Assurance", font = list(color = "black", size = 14, family = "Arial", weight = "bold")), range = c(0, 1)),
          yaxis = list(title = list(text = "Sample Size", font = list(color = "black", size = 14, family = "Arial", weight = "bold"))),
          legend = list(orientation = "v", x = 1.05, y = 0.5),  # Position legend to the right
          title = "Assurance vs Sample Size for the different stopping rules"
        )
      
      p
      
      # plot(IADFOneLook$Assurance, IADFOneLook$`Sample Size`, xlab = "Assurance", xlim = c(0,1),
      #      ylab = "Sample size", ylim = c(min(IADFOneLook$`Sample Size`), FinalAss$`Sample Size`),
      #      main = "Assurance vs Sample Size for the different stopping rules", pch = 19)
      # text(IADFOneLook$Assurance, IADFOneLook$`Sample Size`, IADFOneLook$`Information Fraction`, pos = 2)  
      # points(FinalAss$Assurance, FinalAss$`Sample Size`, col = "red", pch = 19)
      # text(FinalAss$Assurance, FinalAss$`Sample Size`, "No IA", col = "red", pos = 4)
      # 
      # legend("topleft", legend = c("Chosen Rules", "No IA"), col = c("black", "red"), pch = 19)
      # 
      
    })
    
    output$finalAssTable1LookText <- renderUI({
      p(HTML("<b>No Interim Analysis</b>"))
    })
    
    shinyjs::show("finalAssTable1LookText")
    
    
  })
  
  # Two Looks Logic ---------------------------------
  
  observe({
    if (!is.null(reactValues$treatmentSamplesDF)&reactValues$errorSeqTwoLooks==F) {
      updateActionButton(session, "calcTwoLooks", disabled = FALSE)
    } else {
      updateActionButton(session, "calcTwoLooks", disabled = TRUE)
    }
  })
  
  
  observe({
    # Check if any of the inputs are NA
    if (is.na(input$TwoLooksLB1) ||
        is.na(input$TwoLooksUB1) ||
        is.na(input$TwoLooksBy1) ||
        is.na(input$TwoLooksLB2) ||
        is.na(input$TwoLooksUB2) ||
        is.na(input$TwoLooksBy2)) {
      reactValues$errorSeqTwoLooks <- TRUE
    } else {
      # Check validity of input ranges
      if (input$TwoLooksLB1 <= 0 | input$TwoLooksLB1 >= 1 ||
          input$TwoLooksUB1 <= 0 | input$TwoLooksUB1 >= 1 ||
          input$TwoLooksBy1 <= 0 | input$TwoLooksBy1 >= 1 ||
          input$TwoLooksLB2 <= 0 | input$TwoLooksLB2 >= 1 ||
          input$TwoLooksUB2 <= 0 | input$TwoLooksUB2 >= 1 ||
          input$TwoLooksBy2 <= 0 | input$TwoLooksBy2 >= 1 ||
          input$TwoLooksLB1 > input$TwoLooksUB1 ||
          input$TwoLooksLB2 > input$TwoLooksUB2) {
        reactValues$errorSeqTwoLooks <- TRUE
      } else {
        # Check the generated sequences
        TwoLooksSeq1 <- seq(input$TwoLooksLB1, input$TwoLooksUB1, by = input$TwoLooksBy1)
        TwoLooksSeq2 <- seq(input$TwoLooksLB2, input$TwoLooksUB2, by = input$TwoLooksBy2)
        
        if (TwoLooksSeq1[1] >= TwoLooksSeq2[1] ||
            TwoLooksSeq1[length(TwoLooksSeq1)] >= TwoLooksSeq2[length(TwoLooksSeq2)]) {
          reactValues$errorSeqTwoLooks <- TRUE
        } else {
          reactValues$errorSeqTwoLooks <- FALSE
        }
      }
    }
  })
  
  
  observe({
    if (reactValues$errorSeqTwoLooks) {
      shinyjs::show("twoLooksErrorMessage")
    } else {
      shinyjs::hide("twoLooksErrorMessage")
    }
  })
  
  output$twoLooksErrorMessage <- renderText({
    if (reactValues$errorSeqTwoLooks) {
      return("Your inputs are incorrect!")
    } else {
      return("")
    }
  })
  
  output$TwoLooksText1 <- renderText({
    
    if (reactValues$errorSeqTwoLooks==TRUE) {
      return("")
    }
    
    TwoLooksSeq1 <- seq(input$TwoLooksLB1, input$TwoLooksUB1, by = input$TwoLooksBy1)
    value <- paste0("We perform IA1 at: ", paste(TwoLooksSeq1, collapse = ", "))
    return(value)
  })
  
  output$TwoLooksText2 <- renderText({
    
    if (reactValues$errorSeqTwoLooks==TRUE) {
      return("")
    }
    
    TwoLooksSeq2 <- seq(input$TwoLooksLB2, input$TwoLooksUB2, by = input$TwoLooksBy2)
    value <- paste0("We perform IA2 at: ", paste(TwoLooksSeq2, collapse = ", "))
    return(value)
  })
  
  output$spendingTwoLooks <- renderRHandsontable({
    initial_data <- data.frame(
      Stage = 1:3,
      alphaspending = c("0.0083", "0.0167", "0.0250"),
      betaspending = c("0.0667", "0.1333", "0.2000")
    )
    
    colnames(initial_data) <- c("Stage", "Alpha spending", "Beta spending")
    
    rhandsontable(initial_data, rowHeaders = FALSE)  %>%
      hot_col(col = "Stage", readOnly = TRUE)
  })
  
  observe({
    
    if (reactValues$errorSeqTwoLooks==FALSE) {
      
      TwoLooksSeq1 <- seq(input$TwoLooksLB1, input$TwoLooksUB1, by = input$TwoLooksBy1)
      TwoLooksSeq2 <- seq(input$TwoLooksLB2, input$TwoLooksUB2, by = input$TwoLooksBy2)
      
      TwoLooksBoundaryIAChoices <- vector("list", length = sum(outer(TwoLooksSeq1, TwoLooksSeq2, "<")))
      
      myCount <- 1
      
      for (j in 1:length(TwoLooksSeq1)){
        for (k in 1:length(TwoLooksSeq2)){
          if (TwoLooksSeq1[j] < TwoLooksSeq2[k]){
            TwoLooksBoundaryIAChoices[[myCount]] <- c(TwoLooksSeq1[j], TwoLooksSeq2[k])
            myCount <- myCount + 1
          }
        }
      }
      
      TwoLooksBoundaryIAChoices <- sapply(TwoLooksBoundaryIAChoices, function(x) paste(x, collapse = ", "))
      TwoLooksBoundaryIAChoices <- setNames(TwoLooksBoundaryIAChoices, sapply(TwoLooksBoundaryIAChoices, function(x) paste(x, collapse = ", ")))
      updateSelectInput(session, "twoLooksBoundaryIA", choices = TwoLooksBoundaryIAChoices)
    
    }
  })
  
  
  observeEvent(c(input$spendingTwoLooks, input$twoLooksBoundaryIA), {
    values <- hot_to_r(input$spendingTwoLooks)
    
    if (!is.null(values)){
      design <- getDesignGroupSequential(typeOfDesign = "asUser", 
                                         informationRates = as.numeric(c(strsplit(input$twoLooksBoundaryIA, ", ")[[1]], 1)), 
                                         userAlphaSpending = as.numeric(values[,2]), 
                                         typeBetaSpending = "bsUser",
                                         userBetaSpending = as.numeric(values[,3]))
      
      output$twoLooksBoundaries <- renderPlotly({
        plot(design)
      })
    }
    
  })
  
  
  # Reactive expression to calculate metrics
  calculateTablemetrics <- function(iterationList, TwoLooksSeq1, TwoLooksSeq2) {
    
    
    PowerArray <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    DurationArray <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    SSArray <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    IATime1 <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    IATime2 <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    PercentStop <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    PercentStopLook1 <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    PercentStopLook2 <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    PercentStopLook1Fut <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    PercentStopLook2Fut <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    PercentStopLook1Eff <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    PercentStopLook2Eff <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    PercentStopEff <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    PercentStopFut <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    
    myCount <- 1

    for (j in 1:length(TwoLooksSeq1)){
      for (k in 1:length(TwoLooksSeq2)){
        if (TwoLooksSeq1[j]<TwoLooksSeq2[k]){

          PowerArray[k,j] <- mean(iterationList[[myCount]]$Power)
          DurationArray[k,j] <- mean(iterationList[[myCount]]$Duration)
          SSArray[k,j] <- mean(iterationList[[myCount]]$SS)
          IATime1[k,j] <- mean(iterationList[[myCount]]$IA1Time)
          IATime2[k,j] <- mean(iterationList[[myCount]]$IA2Time)
          PercentStop[k,j] <- mean(iterationList[[myCount]]$Stop)
          PercentStopLook1[k,j] <- mean(iterationList[[myCount]]$StopLook1)
          PercentStopLook2[k,j] <- mean(iterationList[[myCount]]$StopLook2)
          PercentStopLook1Fut[k,j] <- mean(iterationList[[myCount]]$Outcome=="Futility1")
          PercentStopLook2Fut[k,j] <- mean(iterationList[[myCount]]$Outcome=="Futility2")
          PercentStopLook1Eff[k,j] <- mean(iterationList[[myCount]]$Outcome=="Efficacy1")
          PercentStopLook2Eff[k,j] <- mean(iterationList[[myCount]]$Outcome=="Efficacy2")
          PercentStopEff[k,j] <- mean(iterationList[[myCount]]$Outcome %in% c("Efficacy1", "Efficacy2"))
          PercentStopFut[k,j] <- mean(iterationList[[myCount]]$Outcome %in% c("Futility1", "Futility2"))

          myCount <- myCount + 1

        }
      }
    }
    
    
    correctlyStopDF <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    correctlyStopEffDF <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    correctlyStopFutDF <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    correctlyStopLook1DF <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    correctlyStopLook2DF <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    correctlyStopEffLook1DF <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    correctlyStopEffLook2DF <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    correctlyStopFutLook1DF <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    correctlyStopFutLook2DF <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    
    correctlyContinueDF <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    correctlyContinueLook1DF <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    correctlyContinueLook2DF <- data.frame(matrix(NA, nrow = length(TwoLooksSeq2), ncol = length(TwoLooksSeq1)))
    
    myCount <- 1
    
    for (j in 1:length(TwoLooksSeq1)){
      for (k in 1:length(TwoLooksSeq2)){
        if (TwoLooksSeq1[j]<TwoLooksSeq2[k]){
          
      
      #Calculating Correctly Stop
      IndStopVec <- which(iterationList[[myCount]]$Stop==1)
      if (length(IndStopVec)!=0){
        
        OutcomeVec <- iterationList[[myCount]]$Outcome[IndStopVec]
        TruthVec <- iterationList[[myCount]]$Truth[IndStopVec]
        correctlyStopSum <- 0
        
        for (i in 1:length(IndStopVec)){
          if ((OutcomeVec[i]=="Efficacy1"&TruthVec[i]==1)|(OutcomeVec[i]=="Futility1"&TruthVec[i]==0)|
              (OutcomeVec[i]=="Efficacy2"&TruthVec[i]==1)|(OutcomeVec[i]=="Futility2"&TruthVec[i]==0)){
            correctlyStopSum <- correctlyStopSum + 1
          }
        }
        
        correctlyStopDF[k,j] <- correctlyStopSum/length(IndStopVec)
        
      } else {
        correctlyStopDF[k,j] <- NA
      }
      
      #Calculating Correctly Stopping for Efficacy
      IndStopEffVec <- which(iterationList[[myCount]]$StopEff==1)
      if (length(IndStopEffVec)!=0){
        
        OutcomeVec <- iterationList[[myCount]]$Outcome[IndStopEffVec]
        TruthVec <- iterationList[[myCount]]$Truth[IndStopEffVec]
        correctlyStopEffSum <- 0
        
        for (i in 1:length(IndStopEffVec)){
          if (TruthVec[i]==1){
            correctlyStopEffSum <- correctlyStopEffSum + 1
          }
        }
        
        correctlyStopEffDF[k,j] <- correctlyStopEffSum/length(IndStopEffVec)
        
      } else {
        correctlyStopEffDF[k,j] <- NA
      }
      
      #Calculating Correctly Stopping for Futility
      IndStopFutVec <- which(iterationList[[myCount]]$StopFut==1)
      if (length(IndStopFutVec)!=0){
        
        OutcomeVec <- iterationList[[myCount]]$Outcome[IndStopFutVec]
        TruthVec <- iterationList[[myCount]]$Truth[IndStopFutVec]
        correctlyStopFutSum <- 0
        
        for (i in 1:length(IndStopFutVec)){
          if (TruthVec[i]==0){
            correctlyStopFutSum <- correctlyStopFutSum + 1
          }
        }
        
        correctlyStopFutDF[k,j] <- correctlyStopFutSum/length(IndStopFutVec)
        
      } else {
        correctlyStopFutDF[k,j] <- NA
      }
      
      
      
      #Calculating Correctly Stopping at Look 1
      IndStopLook1Vec <- which(iterationList[[myCount]]$StopLook1==1)
      if (length(IndStopLook1Vec)!=0){
        
        OutcomeVec <- iterationList[[myCount]]$Outcome[IndStopLook1Vec]
        TruthVec <- iterationList[[myCount]]$Truth[IndStopLook1Vec]
        correctlyStopLook1Sum <- 0
        
        for (i in 1:length(IndStopLook1Vec)){
          if ((OutcomeVec[i]=="Efficacy1"&TruthVec[i]==1)|(OutcomeVec[i]=="Futility1"&TruthVec[i]==0)){
            correctlyStopLook1Sum <- correctlyStopLook1Sum + 1
          }
        }
        
        correctlyStopLook1DF[k,j] <- correctlyStopLook1Sum/length(IndStopLook1Vec)
        
      } else {
        correctlyStopLook1DF[k,j] <- NA
      }
      
      #Calculating Correctly Stopping at Look 2
      IndStopLook2Vec <- which(iterationList[[myCount]]$StopLook2==1)
      if (length(IndStopLook2Vec)!=0){
        
        OutcomeVec <- iterationList[[myCount]]$Outcome[IndStopLook2Vec]
        TruthVec <- iterationList[[myCount]]$Truth[IndStopLook2Vec]
        correctlyStopLook2Sum <- 0
        
        for (i in 1:length(IndStopLook2Vec)){
          if ((OutcomeVec[i]=="Efficacy2"&TruthVec[i]==1)|(OutcomeVec[i]=="Futility2"&TruthVec[i]==0)){
            correctlyStopLook2Sum <- correctlyStopLook2Sum + 1
          }
        }
        
        correctlyStopLook2DF[k,j] <- correctlyStopLook2Sum/length(IndStopLook2Vec)
        
      } else {
        correctlyStopLook2DF[k,j] <- NA
      }
      
      
      #Calculating Correctly Stopping for Efficacy at Look 1
      IndStopEffLook1Vec <- which(iterationList[[myCount]]$Outcome=="Efficacy1")
      if (length(IndStopEffLook1Vec)!=0){
        
        OutcomeVec <- iterationList[[myCount]]$Outcome[IndStopEffLook1Vec]
        TruthVec <- iterationList[[myCount]]$Truth[IndStopEffLook1Vec]
        correctlyStopEffLook1Sum <- 0
        
        for (i in 1:length(IndStopEffLook1Vec)){
          if (TruthVec[i]==1){
            correctlyStopEffLook1Sum <- correctlyStopEffLook1Sum + 1
          }
        }
        
        correctlyStopEffLook1DF[k,j] <- correctlyStopEffLook1Sum/length(IndStopEffLook1Vec)
        
      } else {
        correctlyStopEffLook1DF[k,j] <- NA
      }
      
      #Calculating Correctly Stopping for Efficacy at Look 2
      IndStopEffLook2Vec <- which(iterationList[[myCount]]$Outcome=="Efficacy2")
      if (length(IndStopEffLook2Vec)!=0){
        
        OutcomeVec <- iterationList[[myCount]]$Outcome[IndStopEffLook2Vec]
        TruthVec <- iterationList[[myCount]]$Truth[IndStopEffLook2Vec]
        correctlyStopEffLook2Sum <- 0
        
        for (i in 1:length(IndStopEffLook2Vec)){
          if (TruthVec[i]==1){
            correctlyStopEffLook2Sum <- correctlyStopEffLook2Sum + 1
          }
        }
        
        correctlyStopEffLook2DF[k,j] <- correctlyStopEffLook2Sum/length(IndStopEffLook2Vec)
        
      } else {
        correctlyStopEffLook2DF[k,j] <- NA
      }
      
      #Calculating Correctly Stopping for Futility at Look 1
      IndStopFutLook1Vec <- which(iterationList[[myCount]]$Outcome=="Futility1")
      if (length(IndStopFutLook1Vec)!=0){
        
        OutcomeVec <- iterationList[[myCount]]$Outcome[IndStopFutLook1Vec]
        TruthVec <- iterationList[[myCount]]$Truth[IndStopFutLook1Vec]
        correctlyStopFutLook1Sum <- 0
        
        for (i in 1:length(IndStopFutLook1Vec)){
          if (TruthVec[i]==0){
            correctlyStopFutLook1Sum <- correctlyStopFutLook1Sum + 1
          }
        }
        
        correctlyStopFutLook1DF[k,j] <- correctlyStopFutLook1Sum/length(IndStopFutLook1Vec)
        
      } else {
        correctlyStopFutLook1DF[k,j] <- NA
      }
      
      #Calculating Correctly Stopping for Futility at Look 2
      IndStopFutLook2Vec <- which(iterationList[[myCount]]$Outcome=="Futility2")
      if (length(IndStopFutLook2Vec)!=0){
        
        OutcomeVec <- iterationList[[myCount]]$Outcome[IndStopFutLook2Vec]
        TruthVec <- iterationList[[myCount]]$Truth[IndStopFutLook2Vec]
        correctlyStopFutLook2Sum <- 0
        
        for (i in 1:length(IndStopFutLook2Vec)){
          if (TruthVec[i]==0){
            correctlyStopFutLook2Sum <- correctlyStopFutLook2Sum + 1
          }
        }
        
        correctlyStopFutLook2DF[k,j] <- correctlyStopFutLook2Sum/length(IndStopFutLook2Vec)
        
      } else {
        correctlyStopFutLook2DF[k,j] <- NA
      }
      
      
      
      #Calculating Correctly continuing 
      IndContinueCorrectlyVec <- which(iterationList[[myCount]]$Stop==0)
      if (length(IndContinueCorrectlyVec)!=0){
        
        TruthVec <- iterationList[[myCount]]$Truth[IndContinueCorrectlyVec]
        correctlyContinueSum <- 0
        
        for (i in 1:length(IndContinueCorrectlyVec)){
          if (TruthVec[i]==1){
            correctlyContinueSum <- correctlyContinueSum + 1
          }
        }
        
        correctlyContinueDF[k,j] <- correctlyContinueSum/length(IndContinueCorrectlyVec)
        
      } else {
        correctlyContinueDF <- NA
      }
      
      #Calculating Correctly continuing at Look 1
      IndContinueCorrectlyLook1Vec <- which(iterationList[[myCount]]$StopLook1==0)
      if (length(IndContinueCorrectlyLook1Vec)!=0){
        
        TruthVec <- iterationList[[myCount]]$Truth[IndContinueCorrectlyLook1Vec]
        correctlyContinueLook1Sum <- 0
        
        for (i in 1:length(IndContinueCorrectlyLook1Vec)){
          if (TruthVec[i]==1){
            correctlyContinueLook1Sum <- correctlyContinueLook1Sum + 1
          }
        }
        
        correctlyContinueLook1DF[k,j] <- correctlyContinueLook1Sum/length(IndContinueCorrectlyLook1Vec)
        
      } else {
        correctlyContinueLook1DF <- NA
      }
      
      #Calculating Correctly continuing at Look 2
      IndContinueCorrectlyLook2Vec <- which(iterationList[[myCount]]$StopLook2==0)
      if (length(IndContinueCorrectlyLook2Vec)!=0){
        
        TruthVec <- iterationList[[myCount]]$Truth[IndContinueCorrectlyLook2Vec]
        correctlyContinueLook2Sum <- 0
        
        for (i in 1:length(IndContinueCorrectlyLook2Vec)){
          if (TruthVec[i]==1){
            correctlyContinueLook2Sum <- correctlyContinueLook2Sum + 1
          }
        }
        
        correctlyContinueLook2DF[k,j] <- correctlyContinueLook2Sum/length(IndContinueCorrectlyLook2Vec)
        
      } else {
        correctlyContinueLook2DF <- NA
      }
      
      
      myCount <- myCount + 1
      
        }
      }
      
      }
    
    
    return(list(PowerArray = PowerArray, DurationArray = DurationArray,
                SSArray = SSArray, IATime1 = IATime1, IATime2 = IATime2,
                PercentStop = PercentStop,
                PercentStopLook1 = PercentStopLook1, PercentStopLook2 = PercentStopLook2,
                PercentStopLook1Fut = PercentStopLook1Fut, PercentStopLook2Fut = PercentStopLook2Fut,
                PercentStopLook1Eff = PercentStopLook1Eff, PercentStopLook2Eff = PercentStopLook2Eff,
                PercentStopEff = PercentStopEff, PercentStopFut = PercentStopFut,
                correctlyStopDF = correctlyStopDF, 
                correctlyStopEffDF = correctlyStopEffDF, correctlyStopFutDF = correctlyStopFutDF,
                correctlyStopLook1DF = correctlyStopLook1DF, correctlyStopLook2DF = correctlyStopLook2DF,
                correctlyStopEffLook1DF = correctlyStopEffLook1DF, correctlyStopEffLook2DF = correctlyStopEffLook2DF,
                correctlyStopFutLook1DF = correctlyStopFutLook1DF, correctlyStopFutLook2DF = correctlyStopFutLook2DF,
                correctlyContinueDF = correctlyContinueDF, correctlyContinueLook1DF = correctlyContinueLook1DF,
                correctlyContinueLook2DF = correctlyContinueLook2DF))
                
    
  }
  
  output$metricsTableoutput <- renderUI({
    selected_metrics <- input$selectedOptionstwoLooks
    
    totalChoices <- c("Assurance", "Duration", "Sample Size", 
                      "Interim Analysis 1 Time", "Interim Analysis 2 Time",
                      "% Stop", "% Stop Look 1", "% Stop Look 2",
                      "% Stop Look 1 for Futility", "% Stop Look 2 for Futility",
                      "% Stop Look 1 for Efficacy", "% Stop Look 2 for Efficacy",
                      "% Stop for Efficacy", "% Stop for Futility",
                      "Correctly Stop", "Correctly Stop for Efficacy", "Correctly Stop for Futility",
                      "Correctly Stop at Look 1", "Correctly Stop at Look 2",
                      "Correctly Stop for Efficacy at Look 1", "Correctly Stop for Efficacy at Look 2",
                      "Correctly Stop for Futility at Look 1", "Correctly Stop for Futility at Look 2",
                      "Correctly Continue", "Correctly Continue at Look 1", "Correctly Continue at Look 2")
                      
    
    if (!is.null(selected_metrics) && length(selected_metrics) > 0) {
      chosenIndices <- which(totalChoices %in% selected_metrics)
      
      data <- calculateTablemetrics(reactValues$iterationList, 
                                    reactValues$TwoLooksSeq1,  
                                    reactValues$TwoLooksSeq2)
      
      num_tables <- length(chosenIndices)
      tables_list <- lapply(chosenIndices, function(i) {
        title <- totalChoices[i]  # Title for the table
        table_data <- as.data.frame(data[[i]])
        colnames(table_data) <- reactValues$TwoLooksSeq1
        rownames(table_data) <- reactValues$TwoLooksSeq2
        table_output <- renderDT({
          
          cutoffs <- quantile(na.omit(unlist(table_data)), probs = c(0.25, 0.75))
          
          datatable(table_data, options = list(rowCallback = JS(rowCallback)),
                    rownames = TRUE) %>% formatStyle(
            columns = colnames(table_data),
            backgroundColor = styleInterval(
              c(min(table_data, na.rm = T), cutoffs[1], cutoffs[2], max(table_data, na.rm = T)),
              c('#EDF8E9', '#BAE4B3', '#74C476', '#31A354', '#006D2C')
            )) %>%
            formatSignif(
              columns = colnames(table_data),
              digits = 3
            )
        })
        # Wrap the rendered table in a tagList with the title
        tagList(
          h3(title),  # Use h3 for the title (you can adjust this based on your preference)
          table_output
        )
      })
      
      # Return tables_list directly
      tables_list
    } else {
      
    }
  })
  
  
  
  
  
  
  
  output$corTable <- renderDataTable({
    
    data <-  doCorrelation()$corMatrix
    
   
      
      
      
  
    
    
  })
  
  
  
  
  
  
  

  observeEvent(input$calcTwoLooks, {
    NRep <- 20
    TwoLooksSeq1 <- seq(input$TwoLooksLB1, input$TwoLooksUB1, by = input$TwoLooksBy1)
    TwoLooksSeq2 <- seq(input$TwoLooksLB2, input$TwoLooksUB2, by = input$TwoLooksBy2)
    
    iterationList <- vector("list", length = sum(outer(TwoLooksSeq1, TwoLooksSeq2, "<")))
    
    noLooksDF <- data.frame(Power = numeric(NRep),
                            Duration = numeric(NRep),
                            SS = numeric(NRep))
    
    values <- hot_to_r(input$spendingTwoLooks)
    
    myCount <- 1
    
    for (j in 1:length(TwoLooksSeq1)){
      for (k in 1:length(TwoLooksSeq2)){
        if (TwoLooksSeq1[j]<TwoLooksSeq2[k]){
          
          iterationList[[myCount]]$IF <- c(TwoLooksSeq1[j], TwoLooksSeq2[k], 1) 
          
          design <- getDesignGroupSequential(typeOfDesign = "asUser", 
                                             informationRates = iterationList[[myCount]]$IF,
                                             userAlphaSpending = as.numeric(values[,2]), 
                                             typeBetaSpending = "bsUser",
                                             userBetaSpending = as.numeric(values[,3]))
          
          iterationList[[myCount]]$EffBounds <- design$criticalValues
          iterationList[[myCount]]$futBounds <- design$futilityBounds
          
          myCount <- myCount + 1
          
        }
      }
    }
    
    
    
    conc.probs <- matrix(0, 2, 2)
    conc.probs[1, 2] <- 0.5
    
    
    treatmentSamplesDF <- SHELF::copulaSample(reactValues$treatmentSamplesDF$fit1, reactValues$treatmentSamplesDF$fit2,
                                              cp = conc.probs, n = 1e4, d = reactValues$treatmentSamplesDF$d)
    
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
        
        
        #Perform looks at different Information Fractions
        finalDF <- CensFunc(dataCombined, input$numEvents)
        test <- survdiff(Surv(survival_time, status)~group, data = finalDF$dataCombined)
        coxmodel <- coxph(Surv(survival_time, status)~group, data = finalDF$dataCombined)
        deltad <- as.numeric(exp(coef(coxmodel)))
        
        noLooksDF$Power[i] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
        noLooksDF$Duration[i] <- finalDF$censTime
        noLooksDF$SS[i] <- finalDF$SS
        
        
        for (k in 1:length(iterationList)){
          
          GSDOC <- GSDTwoIAFunc(dataCombined, iterationList[[k]]$futBounds, 
                                iterationList[[k]]$EffBounds, iterationList[[k]]$IF, input$numEvents) 
          
          iterationList[[k]]$Power[i] <- ifelse(GSDOC$Outcome %in% c("Efficacy1", "Efficacy2", "Successful"), 1, 0)
          iterationList[[k]]$Duration[i] <- GSDOC$Duration
          iterationList[[k]]$SS[i] <- GSDOC$SS
          iterationList[[k]]$IA1Time[i] <- GSDOC$IA1Time
          iterationList[[k]]$IA2Time[i] <- GSDOC$IA2Time
          iterationList[[k]]$Outcome[i] <- GSDOC$Outcome
          iterationList[[k]]$Stop[i] <- ifelse(GSDOC$Outcome %in% c("Efficacy1", "Futility1", "Efficacy2", "Futility2"), 1, 0)
          iterationList[[k]]$StopLook1[i] <- ifelse(GSDOC$Outcome %in% c("Efficacy1", "Futility1"), 1, 0)
          iterationList[[k]]$StopLook2[i] <- ifelse(GSDOC$Outcome %in% c("Efficacy2", "Futility2"), 1, 0)
          iterationList[[k]]$StopEff[i] <- ifelse(GSDOC$Outcome %in% c("Efficacy1", "Efficacy2"), 1, 0)
          iterationList[[k]]$StopFut[i] <- ifelse(GSDOC$Outcome %in% c("Futility1", "Futility2"), 1, 0)
          iterationList[[k]]$Truth[i] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
        }
        
        incProgress(1/NRep)
      }
      
    })
    
    reactValues$iterationList <- iterationList
    reactValues$TwoLooksSeq1 <- TwoLooksSeq1
    reactValues$TwoLooksSeq2 <- TwoLooksSeq2
    
    
    #Do proposed rule logic
    proposedDF$power <- proposedDF$Look1Power*proposedDF$Look2Power*proposedDF$FinalLookPower
    proposedDF$SS <- ifelse(proposedDF$Look1Power==0, proposedDF$Look1SS, ifelse(proposedDF$Look2Power==0, 
                                                                                 proposedDF$Look2SS, proposedDF$FinalLookSS))
    proposedDF$Duration <- ifelse(proposedDF$Look1Power==0, proposedDF$Look1Duration, ifelse(proposedDF$Look2Power==0, 
                                                                                             proposedDF$Look2Duration, proposedDF$FinalLookDuration))

    shinyjs::show("selectedOptionstwoLooks")
    
    updateSelectizeInput(session, "selectedOptionstwoLooks", choices = c("Assurance", "Duration", "Sample Size", 
                                                                         "Interim Analysis 1 Time", "Interim Analysis 2 Time",
                                                                         "% Stop", "% Stop Look 1", "% Stop Look 2",
                                                                         "% Stop Look 1 for Futility", "% Stop Look 2 for Futility",
                                                                         "% Stop Look 1 for Efficacy", "% Stop Look 2 for Efficacy",
                                                                         "% Stop for Efficacy", "% Stop for Futility",
                                                                         "Correctly Stop", "Correctly Stop for Efficacy", "Correctly Stop for Futility",
                                                                         "Correctly Stop at Look 1", "Correctly Stop at Look 2",
                                                                         "Correctly Stop for Efficacy at Look 1", "Correctly Stop for Efficacy at Look 2",
                                                                         "Correctly Stop for Futility at Look 1", "Correctly Stop for Futility at Look 2",
                                                                         "Correctly Continue", "Correctly Continue at Look 1", "Correctly Continue at Look 2"), 
                         selected = c("Assurance", "Duration", "Sample Size"))
    

    
    #Making the proposed DF correctly
    FinalProposedDF <- data.frame(Assurance = mean(proposedDF$power),
                                  Duration = mean(proposedDF$Duration),
                                  SS = mean(proposedDF$SS))
    
    colnames(FinalProposedDF) <- c("Assurance", "Duration", "Sample Size")
    
    
    FinalAss <- t(colMeans(noLooksDF))
    FinalAss <- as.data.frame(FinalAss)
    colnames(FinalAss) <- c("Assurance", "Duration", "Sample Size")
    
    
    
    # output$finalAssTable2Looks <- renderTable({
    #   
    #   FinalAss
    #   
    # }, digits = 3)
    # 
    # output$proposedTable2Looks <- renderTable({
    #   
    #   FinalProposedDF
    #   
    # }, digits = 3)
    # 
    # 
    
    
    data <- calculateTablemetrics(reactValues$iterationList, 
                                  reactValues$TwoLooksSeq1,  
                                  reactValues$TwoLooksSeq2)
    
    myDFLength <- sum(outer(TwoLooksSeq1, TwoLooksSeq2, "<"))
    
    myDF <- data.frame(IF1 = rep(NA, myDFLength), IF2 = rep(NA, myDFLength), 
                       Power = rep(NA, myDFLength), SampleSize = rep(NA, myDFLength),
                       Duration = rep(NA, myDFLength))
    

    myCount <- 1
    
    for (j in 1:length(TwoLooksSeq1)){
      for (k in 1:length(TwoLooksSeq2)){
        if (TwoLooksSeq1[j]<TwoLooksSeq2[k]){
          
         myDF[myCount,]$IF1 <- TwoLooksSeq1[j]
         myDF[myCount,]$IF2 <- TwoLooksSeq2[k]
         myDF[myCount,]$Power <- data$PowerArray[k,j] 
         myDF[myCount,]$SampleSize <- data$SSArray[k,j] 
         myDF[myCount,]$Duration <- data$DurationArray[k,j] 

         myCount <- myCount + 1
          
        }
      }
    }
    

    output$twoLooksPlotDuration <- renderPlotly({
      
      
      p <- plot_ly(myDF, x = ~Power, y = ~Duration,
                   text = ~ paste0("Information Fraction = ", IF1, ", ", IF2), mode = "markers",
                   type = "scatter", marker = list(size = 10, color = "blue"), name = "Chosen Rules") %>%
        add_trace(x = ~FinalAss$Assurance, y = ~FinalAss$Duration, type = "scatter", mode = "markers",
                  marker = list(size = 10, color = "red"),
                  text = ~ "No Interim Analysis",
                  name = "No Interim Analysis") %>%
        add_trace(x = ~FinalProposedDF$Assurance, y = ~FinalProposedDF$Duration, type = "scatter", mode = "markers",
                  marker = list(size = 10, color = "green"),
                  text = ~ "Proposed Rule",
                  name = "Proposed Rule") %>%
        layout(
          xaxis = list(title = list(text = "Assurance", font = list(color = "black", size = 14, family = "Arial", weight = "bold")), range = c(0, 1)),
          yaxis = list(title = list(text = "Duration", font = list(color = "black", size = 14, family = "Arial", weight = "bold"))),
          legend = list(orientation = "v", x = 1.05, y = 0.5),  # Position legend to the right
          title = "Assurance vs Duration for the different stopping rules"
        )
      
      p
      

      # plot(0, 0, type = "n", ylim = c(min(data$DurationArray, na.rm = TRUE), FinalAss$Duration), xlim = c(0, 1),
      #      xlab = "Power", ylab = "Sample Size")
      # 
      # # Add points for each pair of corresponding elements
      # for (i in 1:nrow(data$PowerArray)) {
      #   for (j in 1:ncol(data$DurationArray)) {
      #     if (!is.na(data$PowerArray[i, j]) && !is.na(data$DurationArray[i, j])) {
      #       points(data$PowerArray[i, j], data$DurationArray[i, j], col = "black", pch = 19)
      #     }
      #   }
      # }
      # 
      # points(FinalProposedDF$Assurance, FinalProposedDF$Duration, col = "blue", pch = 19)
      # 
      # points(FinalAss$Assurance, FinalAss$Duration, col = "red", pch = 19)
      # 
      # legend("topleft", legend = c("Chosen Rules", "Proposed Rule", "No IA"), col = c("black",  "blue", "red"), pch = 19)

    })
    

    output$twoLooksPlotSS <- renderPlotly({
      
      
      p <- plot_ly(myDF, x = ~Power, y = ~SampleSize,
                   text = ~ paste0("Information Fraction = ", IF1, ", ", IF2), mode = "markers",
                   type = "scatter", marker = list(size = 10, color = "blue"), name = "Chosen Rules") %>%
        add_trace(x = ~FinalAss$Assurance, y = ~FinalAss$`Sample Size`, type = "scatter", mode = "markers",
                  marker = list(size = 10, color = "red"),
                  text = ~ "No Interim Analysis",
                  name = "No Interim Analysis") %>%
        add_trace(x = ~FinalProposedDF$Assurance, y = ~FinalProposedDF$`Sample Size`, type = "scatter", mode = "markers",
                  marker = list(size = 10, color = "green"),
                  text = ~ "Proposed Rule",
                  name = "Proposed Rule") %>%
        layout(
          xaxis = list(title = list(text = "Assurance", font = list(color = "black", size = 14, family = "Arial", weight = "bold")), range = c(0, 1)),
          yaxis = list(title = list(text = "Sample Size", font = list(color = "black", size = 14, family = "Arial", weight = "bold"))),
          legend = list(orientation = "v", x = 1.05, y = 0.5),  # Position legend to the right
          title = "Assurance vs Sample Size for the different stopping rules"
        )
      
      p
      

      # plot(0, 0, type = "n", ylim = c(min(data$SSArray, na.rm = TRUE), FinalAss$`Sample Size`), xlim = c(0, 1),
      #      xlab = "Power", ylab = "Sample Size")
      # 
      # # Add points for each pair of corresponding elements
      # for (i in 1:nrow(data$PowerArray)) {
      #   for (j in 1:ncol(data$SSArray)) {
      #     if (!is.na(data$PowerArray[i, j]) && !is.na(data$SSArray[i, j])) {
      #       points(data$PowerArray[i, j], data$SSArray[i, j], col = "black", pch = 19)
      #     }
      #   }
      # }
      # 
      # points(FinalProposedDF$Assurance, FinalProposedDF$`Sample Size`, col = "blue", pch = 19)
      # 
      # points(FinalAss$Assurance, FinalAss$`Sample Size`, col = "red", pch = 19)
      # 
      # legend("topleft", legend = c("Chosen Rules", "Proposed Rule", "No IA"), col = c("black",  "blue", "red"), pch = 19)


    })
   
    
  })
  
  
  # Bayesian Logic ---------------------------------
  
  observe({
    if (!is.null(reactValues$treatmentSamplesDF)&reactValues$errorBayesian==F) {
      updateActionButton(session, "calcBayesian", disabled = FALSE)
    } else {
      updateActionButton(session, "calcBayesian", disabled = TRUE)
    }
  })
  
  
  observe({
    # Check if any of the inputs are NA
    if (is.na(input$IFBayesian) ||
        is.na(input$tEffBayesian)){
      reactValues$errorBayesian <- TRUE
    } else {
      # Check validity of input ranges
      if (input$IFBayesian <= 0 | input$IFBayesian >= 1 ||
          input$tEffBayesian <= 0) {
        reactValues$errorBayesian <- TRUE
      } else {
        reactValues$errorBayesian <- FALSE
      }
    }
  })
  
  observe({
    if (reactValues$errorBayesian==TRUE){
      shinyjs::show("bayesianErrorMessage")
    } else{
      shinyjs::hide("bayesianErrorMessage")
    }
  })
  
  output$bayesianErrorMessage <- renderText({
    if (reactValues$errorBayesian==TRUE){
      return("Your inputs are incorrect!")
    } 
    
    return("")
    
  })
  
  observeEvent(input$calcBayesian, {
    
    #Parallel: # Set up parallel processing
    cl <- makeCluster(detectCores())  # Use all available cores
    registerDoParallel(cl)
    
    NRep <- 50
    
    conc.probs <- matrix(0, 2, 2)
    conc.probs[1, 2] <- 0.5
    
    # Extract required input values
    numPatients <- input$numPatients
    lambdac <- input$lambdac
    recTime <- input$recTime
    numEvents <- input$numEvents
    IFBayesian <- input$IFBayesian
    tEffBayesian <- input$tEffBayesian
    
    treatmentSamplesDF <- SHELF::copulaSample(reactValues$treatmentSamplesDF$fit1, reactValues$treatmentSamplesDF$fit2,
                                              cp = conc.probs, n = 1e4, d = reactValues$treatmentSamplesDF$d)
    
    BPPVec <- foreach(i = 1:NRep, .combine = c, .export = c("SimDTEDataSet", "CensFunc", "BPPFunc"),
                      .packages = c("survival", "rjags", "dplyr")) %dopar% {
                        
                        # Compute treatment times
                        HRStar <- sample(treatmentSamplesDF[,2], 1)
                        bigT <- sample(treatmentSamplesDF[,1], 1)
                        
                        # Simulate control and treatment data
                        dataCombined <- SimDTEDataSet(numPatients, lambdac, bigT, HRStar, recTime)  
                        
                        # Perform looks at different Information Fractions
                        finalDF <- CensFunc(dataCombined, numEvents)
                        
                        test <- survdiff(Surv(survival_time, status) ~ group, data = finalDF$dataCombined)
                        coxmodel <- coxph(Surv(survival_time, status) ~ group, data = finalDF$dataCombined)
                        deltad <- as.numeric(exp(coef(coxmodel)))
                        
                        BPPOutcome <- BPPFunc(dataCombined, numPatients, numEvents * IFBayesian, numEvents, recTime, tEffBayesian)
                        
                        Success <- (test$chisq > qchisq(0.95, 1) & deltad<1)
                        
                        return(list(BPP = BPPOutcome$BPP, Success = Success, propEffect = BPPOutcome$propEffect))
                      }
    
    stopCluster(cl)  # Stop the cluster
    
    BPPVec <- data.frame(BPP = unlist(BPPVec[seq(1, length(BPPVec), by = 3)]),
                         Success = unlist(BPPVec[seq(2, length(BPPVec), by = 3)]),
                         propEffect = unlist(BPPVec[seq(3, length(BPPVec), by = 3)]))
    
    
    
    output$BayesianPlot <- renderPlot({
      
      # Plotting histogram colored by ColorVar
      ggplot(BPPVec, aes(x = BPP, fill = Success)) +
        geom_histogram(position = "identity", alpha = 0.5) +
        scale_x_continuous(limits = c(0, 1)) + xlab("Bayesian Predictive Probability")
      
      
    })
    
    
    output$BayesianEffPlot <- renderPlot({
      
      # Plotting histogram colored by ColorVar
      ggplot(BPPVec, aes(x = propEffect, fill = Success)) +
        geom_histogram(position = "identity", alpha = 0.5) +
        scale_x_continuous(limits = c(0, 1)) + xlab("Proportion less than target effect")
      
    })
    
    output$BayesianBPPvTE <- renderPlot({
      
      plot(BPPVec$BPP, BPPVec$propEffect, xlim = c(0,1), ylim = c(0,1),
           xlab = "Bayesian Predictive Probability", ylab = "Proportion less than target effect")
      abline(a = 0, b = 1, lty = 2)
      
    })
    
    
  })
  
  
  # Report Logic ---------------------------------
  
  output$downloadHTML <- downloadHandler(
    filename = function() {
      paste("report-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      # Path to the Rmd file
      rmd_file <- "report.Rmd"
      
      # Render the Rmd file to HTML with parameters
      rmarkdown::render(
        input = rmd_file,
        output_file = file,
        params = list(dataset = cars),
        envir = new.env(parent = globalenv())
      )
    }
  )
  

}

# Run the Shiny app
shinyApp(ui, server)
