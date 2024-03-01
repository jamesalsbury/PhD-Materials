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
  titlePanel("Bayesian Interim Analyses: Delayed Treatment Effects"),
  
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
    
    
    # Create empty data frames
    assDF <- SSDF <- DurationDF <- IATimeDF <- data.frame(matrix(NA, nrow = NRep, ncol = length(futilityVec)))
    stopDF <- data.frame(matrix(0, nrow = NRep, ncol = length(futilityVec)))
    colnames(assDF) <- colnames(SSDF) <- colnames(DurationDF) <- paste0(futilityVec)
    
    FinalAss <- data.frame(Assurance = numeric(NRep), SS = numeric(NRep), Duration = numeric(NRep))
    
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
        assDF[i,] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
        SSDF[i,] <- finalDF$SS
        DurationDF[i,] <- finalDF$censTime
        
        #Input into the "final" assurance table
        FinalAss[i,1] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
        FinalAss[i,2] <- finalDF$SS
        FinalAss[i,3] <- finalDF$censTime
        
        
        for (k in 1:length(futilityVec)){
          futilityCens <- CensFunc(dataCombined, input$numEvents*futilityVec[k])
          futilityLook <- interimLookFunc(futilityCens$dataCombined, input$OneLookHR)
          IATimeDF[i,k] <- futilityCens$censTime
          if (futilityLook=="Stop"){
            stopDF[i,k] <- 1
            assDF[i,k] <- 0
            SSDF[i,k] <- futilityCens$SS
            DurationDF[i,k] <- futilityCens$censTime
          } 
        }
        incProgress(1/NRep)
      }
    })
    
    pHat <- colMeans(assDF)
    
    LB <- pHat - 1.96 * sqrt(pHat * (1 - pHat) / NRep)
    
    UB <- pHat + 1.96 * sqrt(pHat * (1 - pHat) / NRep)
    
    
    futilityDF <- data.frame(IF = futilityVec,
                             censTime = colMeans(IATimeDF),
                             ass = paste0(round(pHat, 3), " [", round(LB, 3), ",", round(UB, 3), "]"),
                             stop = colMeans(stopDF),
                             SS = colMeans(SSDF),
                             Duration = colMeans(DurationDF))
    
    colnames(futilityDF) <- c("Information Fraction", "IA Time", "Assurance", "% stop","Sample Size", "Duration")
    
    output$futilityTable <- renderTable({
      futilityDF
    }, digits = 3)
    
    FinalAss <- colMeans(FinalAss)
    
    FinalAss <- data.frame(Assurance = FinalAss[1], SS = FinalAss[2], Duration = FinalAss[3])
    
    colnames(FinalAss) <- c("Assurance", "Sample Size", "Duration")
    
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
    NRep <- 500
    TwoLooksSeq1 <- seq(input$TwoLooksLB1, input$TwoLooksUB1, by = input$TwoLooksBy1)
    TwoLooksSeq2 <- seq(input$TwoLooksLB2, input$TwoLooksUB2, by = input$TwoLooksBy2)
    TwoLooksCombined <- c(TwoLooksSeq1, TwoLooksSeq2)
    
    #print(TwoLooksCombined)
    
    ass3D <- SS3d <- Duration3d <- array(0, dim = c(length(TwoLooksSeq2), length(TwoLooksSeq1), NRep))
    
    treatmentSamplesDF <- data$treatmentSamplesDF
    
    #Compute treatment times
    HRStar <- sample(treatmentSamplesDF[,2], NRep, replace = T)
    bigT <- sample(treatmentSamplesDF[,1], NRep, replace = T)
    
    withProgress(message = 'Calculating', value = 0, {
      for (i in 1:NRep){
        
        
        futilityDF1 <- data.frame(futility = c(TwoLooksSeq1),
                                 outcome = numeric(length(TwoLooksSeq1)),
                                 SS = numeric(length(TwoLooksSeq1)),
                                 duration = numeric(length(TwoLooksSeq1)))
        
        futilityDF2 <- data.frame(futility = c(TwoLooksSeq2),
                                  outcome = numeric(length(TwoLooksSeq2)),
                                  SS = numeric(length(TwoLooksSeq2)),
                                  duration = numeric(length(TwoLooksSeq2)))
        
        
        #Simulate control and treatment data
        dataCombined <- SimDTEDataSet(input$numPatients, input$lambdac, bigT[i], HRStar[i], input$recTime)  
        
        #Perform futility look at different Information Fractions
        finalDF <- CensFunc(dataCombined, input$numEvents)
        test <- survdiff(Surv(survival_time, status)~group, data = finalDF$dataCombined)
        coxmodel <- coxph(Surv(survival_time, status)~group, data = finalDF$dataCombined)
        deltad <- as.numeric(exp(coef(coxmodel)))
        FinalOutcome <- (test$chisq > qchisq(0.95, 1) & deltad<1)
        FinalSS <- finalDF$SS
        FinalDuration <- finalDF$censTime
        
        for (k in 1:length(TwoLooksSeq1)){
          futilityCens <- CensFunc(dataCombined, input$numEvents*TwoLooksSeq1[k])
          futilityLook <- interimLookFunc(futilityCens$dataCombined, input$TwoLooksHR1)
          futilityDF1[k,3] <- futilityCens$SS
          futilityDF1[k,4] <- futilityCens$censTime
          if (futilityLook!="Stop"){ futilityDF1[k,2] <- 1} 
        }
        
        for (k in 1:length(TwoLooksSeq2)){
          futilityCens <- CensFunc(dataCombined, input$numEvents*TwoLooksSeq2[k])
          futilityLook <- interimLookFunc(futilityCens$dataCombined, input$TwoLooksHR2)
          futilityDF2[k,3] <- futilityCens$SS
          futilityDF2[k,4] <- futilityCens$censTime
          if (futilityLook!="Stop"){ futilityDF2[k,2] <- 1} 
        }
        

        
          for (j in 1:length(TwoLooksSeq1)){
            for (k in j:length(TwoLooksSeq2)){
              Look1Value <- TwoLooksSeq1[j]
              Look2Value <- TwoLooksSeq2[k]
              Look1Row <- match(round(as.numeric(Look1Value), 1), round(as.numeric(futilityDF1$futility), 1))
              Look2Row <- match(round(as.numeric(Look2Value), 1), round(as.numeric(futilityDF2$futility), 1))
              
              
              Look1Outcome <- futilityDF1[Look1Row,]$outcome
              Look2Outcome <- futilityDF2[Look2Row,]$outcome

              #Assurance logic
              if (Look1Outcome == 1 && Look2Outcome == 1 && FinalOutcome == 1) {
                ass3D[k, j, i] <- 1
              }
              
              #Sample size & duration logic
              if (Look1Outcome == 1 & Look2Outcome == 1) {
                SS3d[k,j,i] <- FinalSS
                Duration3d[k, j, i] <- FinalDuration
              } else if (Look1Outcome == 0) {
                SS3d[k,j,i] <- futilityDF1[Look1Row,]$SS
                Duration3d[k, j, i] <- futilityDF1[Look1Row,]$duration
              } else if (Look1Outcome == 1 & Look2Outcome == 0) {
                SS3d[k,j,i] <- futilityDF2[Look2Row,]$SS
                Duration3d[k, j, i] <- futilityDF2[Look2Row,]$duration
              }
              
            }
          }
        
        incProgress(1/NRep)
      }
      
      
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
