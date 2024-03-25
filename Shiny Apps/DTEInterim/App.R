library(shiny)
library(shinyjs)
library(ggplot2)
library(survival)
library(dplyr)

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
                   plotOutput("TSamples"),
                   plotOutput("HRStarSamples")
                   
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
                   # Generate plotOutputs dynamically based on the number of delay and HR combinations
                  tableOutput("futilityTable"),
                  tableOutput("finalAssTable")
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
                   hidden(uiOutput("AssTable")),
                   tableOutput("twoLooksAss"),
                   hidden(uiOutput("SSTable")),
                   tableOutput("twoLooksSS"),
                   hidden(uiOutput("DurationTable")),
                   tableOutput("twoLooksDuration")
                 )
               )
      ),
      
      # Bayesian UI ---------------------------------
      
      tabPanel("Bayesian", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   actionButton("calcFutilityBayesian", label  = "Calculate", disabled = T)
                 ), 
                 mainPanel = mainPanel(
                   # Generate plotOutputs dynamically based on the number of delay and HR combinations
                   tableOutput("BayesianAss"),
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
    file <- input$samplesFile
    if (!is.null(file)) {
      # Read the uploaded file into treatmentSamplesDF
      data$treatmentSamplesDF <- read.csv(file$datapath)

      # Enable the action button when a file is selected
      shinyjs::enable("calcFutilityOneLook")
      shinyjs::enable("calcFutilityTwoLooks")
      
      
      # output$TSamples <- renderPlot({
      #   hist(data$treatmentSamplesDF$bigT, freq = F, main  = "Histogram of T", xlab = "Time")
      # })
      # 
      # output$HRStarSamples <- renderPlot({
      #   hist(data$treatmentSamplesDF$HRStar, freq = F, main = "Histogram of post-delay HR", xlab = "Post-delay HR")
      # })
      # 
      # output$P_S <- renderText({
      # 
      #   paste0("P_S is estimated to be: ",
      #          nrow(subset(data$treatmentSamplesDF, HRStar != 1)) / nrow(data$treatmentSamplesDF))
      # 
      # 
      # })
      # 
      # output$P_DTE <- renderText({
      # 
      #   separateDF <- data$treatmentSamplesDF %>%
      #     filter(HRStar!=1)
      # 
      # 
      # 
      #   paste0("P_DTE is estimated to be: ",
      #          nrow(subset(separateDF, bigT != 0)) / nrow(separateDF))
      # 
      # 
      # })
      
      
      
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
    
    withProgress(message = 'Calculating', value = 0, {
      for (i in 1:NRep){
        
        treatmentSamplesDF <- data$treatmentSamplesDF
        
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
    
    FinalAss <- mean(iterationArray[1, length(futilityVec)+1, ])
    FinalDuration <- mean(iterationArray[2, length(futilityVec)+1, ])
    FinalSS <- mean(iterationArray[3, length(futilityVec)+1, ])
  
    
    FinalAss <- data.frame(Assurance = FinalAss, Duration = FinalDuration, SS = FinalSS)

    colnames(FinalAss) <- c("Assurance", "Duration", "Sample Size")

    output$finalAssTable <- renderTable({
      FinalAss
    }, digits = 3)
    
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
    NRep <- 3
    TwoLooksSeq1 <- seq(input$TwoLooksLB1, input$TwoLooksUB1, by = input$TwoLooksBy1)
    TwoLooksSeq2 <- seq(input$TwoLooksLB2, input$TwoLooksUB2, by = input$TwoLooksBy2)
    TwoLooksCombined <- unique(c(round(TwoLooksSeq1, 2), round(TwoLooksSeq2, 2)))
    
    #print(TwoLooksCombined)
    
    iterationArray <- array(NA,  dim = c(3, length(TwoLooksCombined)+1, NRep))
    
    dimnames(iterationArray) <- list(c("Power", "Duration", "Sample Size"), c(TwoLooksCombined, "Final"), 1:NRep)
    
    #ass3D <- SS3d <- Duration3d <- array(0, dim = c(length(TwoLooksSeq2), length(TwoLooksSeq1), NRep))
    
    treatmentSamplesDF <- data$treatmentSamplesDF
    
    withProgress(message = 'Calculating', value = 0, {
      for (i in 1:NRep){
        
        treatmentSamplesDF <- data$treatmentSamplesDF
        
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
      
      print(iterationArray)
      
      # Define a function to create the matrix
      create_matrix <- function(data3D, row_names, col_names) {
        mean_ij <- function(i, j) {
          mean(data3D[i, j, ])
        }
        
        # Create a matrix to store the means
        matrix_means <- matrix(0, nrow = length(row_names), ncol = length(col_names))
        
        # Apply the mean_ij function to each combination of i and j
        vector_means <- apply(expand.grid(i = 1:length(row_names), j = 1:length(col_names)), 1, function(idx) {
          mean_ij(idx[1], idx[2])
        })
        
        matrix_means[] <- vector_means
        
        colnames(matrix_means) <- col_names
        rownames(matrix_means) <- row_names
        
        return(matrix_means)
      }
      
      # Create matrices
      ass_matrix <- create_matrix(ass3D, TwoLooksSeq2, TwoLooksSeq1)
      ss_matrix <- create_matrix(SS3d, TwoLooksSeq2, TwoLooksSeq1)
      duration_matrix <- create_matrix(Duration3d, TwoLooksSeq2, TwoLooksSeq1)
      
      
      
      output$twoLooksAss <- renderTable({
        ass_matrix <- round(ass_matrix, 3)
        
        ass_matrix[upper.tri(ass_matrix)] <- ""
        
        ass_matrix
      }, rownames = TRUE)
      
      
      output$twoLooksSS <- renderTable({
        ss_matrix <- round(ss_matrix, 3)
        
        ss_matrix[upper.tri(ss_matrix)] <- ""
        
        ss_matrix
      }, rownames = T)
      
      output$twoLooksDuration <- renderTable({
        
        duration_matrix <- round(duration_matrix, 3)
        
        duration_matrix[upper.tri(duration_matrix)] <- ""
        
        
        duration_matrix
      }, rownames = T)
      
    })
    
    
    output$AssTable <- renderUI({
      p(HTML("<b>Assurance table</b>"))
    })
    
    output$SSTable <- renderUI({
      p(HTML("<b>Sample size table</b>"))
    })
    
    output$DurationTable <- renderUI({
      p(HTML("<b>Duration table</b>"))
    })
    
    shinyjs::show("AssTable")
    shinyjs::show("SSTable")
    shinyjs::show("DurationTable")
    
    
    
    
  })

}

# Run the Shiny app
shinyApp(ui, server)
