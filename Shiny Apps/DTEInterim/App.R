library(shiny)
library(shinyjs)
library(ggplot2)
library(survival)

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
                   
                 )
               )
      ),
      # One Look UI ---------------------------------
      
      tabPanel("One Look", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   actionButton("calcFutilityOneLook", label  = "Calculate", disabled = T)
                 ), 
                 mainPanel = mainPanel(
                   # Generate plotOutputs dynamically based on the number of delay and HR combinations
                  tableOutput("futilityTable")
                 )
               )
      ),
      
      # Two Looks UI ---------------------------------
      
      tabPanel("Two Looks", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   actionButton("calcFutilityTwoLooks", label  = "Calculate", disabled = T)
                 ), 
                 mainPanel = mainPanel(
                   # Generate plotOutputs dynamically based on the number of delay and HR combinations
                   tableOutput("twoLooksAss")
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
    #print(data$treatmentSamplesDF)
      # Enable the action button when a file is selected
      shinyjs::enable("calcFutilityOneLook")
      shinyjs::enable("calcFutilityTwoLooks")
    }
  })
  
  
  # One Look Logic ---------------------------------
  
  observeEvent(input$calcFutilityOneLook, {
    NRep <- 500
    futilityVec <- seq(0.2, 1, by = 0.1)
    
    # Create empty data frames
    assDF <- SSDF <- DurationDF <- IATimeDF <- data.frame(matrix(NA, nrow = NRep, ncol = length(futilityVec)))
    stopDF <- data.frame(matrix(0, nrow = NRep, ncol = length(futilityVec)))
    colnames(assDF) <- colnames(SSDF) <- colnames(DurationDF) <- paste0(futilityVec)
    
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
        
        for (k in 1:length(futilityVec)){
          futilityCens <- CensFunc(dataCombined, input$numEvents*futilityVec[k])
          futilityLook <- interimLookFunc(futilityCens$dataCombined)
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
  })
  
  
  # Two Looks Logic ---------------------------------
  
  
  observeEvent(input$calcFutilityTwoLooks, {
    NRep <- 500
    Look1Vec <- seq(0.2, 0.8, by = 0.1)
    Look2Vec <- seq(0.3, 0.9, by = 0.1)
    futilityVec <- seq(0.2, 0.9, by = 0.1)
    
    
    # # Create empty data frames
    # assDF <- SSDF <- DurationDF <- IATimeDF <- data.frame(matrix(NA, nrow = NRep, ncol = length(futilityVec)))
    # stopDF <- data.frame(matrix(0, nrow = NRep, ncol = length(futilityVec)))
    # colnames(assDF) <- colnames(SSDF) <- colnames(DurationDF) <- paste0(futilityVec)
    
    ass3D <- array(0, dim = c(length(Look2Vec), length(Look1Vec), NRep))
    
    withProgress(message = 'Calculating', value = 0, {
      for (i in 1:NRep){
        
        
        futilityDF <- data.frame(futility = c(futilityVec,1),
                                 outcome = numeric(length(futilityVec)+1),
                                 SS = numeric(length(futilityVec)+1),
                                 duration = numeric(length(futilityVec)+1))
        
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
        futilityDF[nrow(futilityDF),2] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
        futilityDF[nrow(futilityDF),3] <- finalDF$SS
        futilityDF[nrow(futilityDF),4] <- finalDF$censTime
        
        for (k in 1:length(futilityVec)){
          futilityCens <- CensFunc(dataCombined, input$numEvents*futilityVec[k])
          futilityLook <- interimLookFunc(futilityCens$dataCombined)
          futilityDF[k,3] <- futilityCens$SS
          futilityDF[k,4] <- futilityCens$censTime
          if (futilityLook!="Stop"){ futilityDF[k,2] <- 1} 
        }
        
        
        if (futilityDF[nrow(futilityDF),]$outcome!=0){
          for (j in 1:length(Look1Vec)){
            for (k in j:length(Look2Vec)){
              Look1Value <- Look1Vec[j]
              Look2Value <- Look2Vec[k]
              Look1Row <- match(round(as.numeric(Look1Value), 1), round(as.numeric(futilityDF$futility), 1))
              Look2Row <- match(round(as.numeric(Look2Value), 1), round(as.numeric(futilityDF$futility), 1))
              
              #Assurance logic
              Look1Outcome <- futilityDF[Look1Row,]$outcome
              Look2Outcome <- futilityDF[Look2Row,]$outcome
              FinalOutcome <- futilityDF[nrow(futilityDF),]$outcome
              
              if (Look1Outcome == 1 && Look2Outcome == 1 && FinalOutcome == 1) {
                ass3D[k, j, i] <- 1
              }
              
              #Sample size logic
              
            }
          }
        }
        
        
        incProgress(1/NRep)
      }
    
    

    mean_ij <- function(i, j) {
      mean(ass3D[i, j, ])
    }
    
    # Create a matrix to store the means
    means_matrix <- matrix(0, nrow = length(Look2Vec), ncol = length(Look1Vec))
    
    # Apply the mean_ij function to each combination of i and j
    means_vector <- apply(expand.grid(i = 1:length(Look2Vec), j = 1:length(Look1Vec)), 1, function(idx) {
      mean_ij(idx[1], idx[2])
    })
    
    means_matrix[] <- means_vector
    
    colnames(means_matrix) <- Look1Vec
    rownames(means_matrix) <- Look2Vec
    
    

    output$twoLooksAss <- renderTable({
      means_matrix
    }, digits = 3, rownames = T)
    
  })
  
  
  
  })

}

# Run the Shiny app
shinyApp(ui, server)
