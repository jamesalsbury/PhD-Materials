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
                   # Generate plotOutputs dynamically based on the number of delay and HR combinations
                   tableOutput("futilityTable")
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
      )
      
      
    ), style='width: 1000px; height: 600px'
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Output Logic ---------------------------------
  

  data <- reactiveValues(treatmentSamplesDF = NULL)
  
  observeEvent(input$samplesFile, {
    file <- input$samplesFile
    if (!is.null(file)) {
      # Read the uploaded file into treatmentSamplesDF
      data$treatmentSamplesDF <- read.csv(file$datapath)
      
      # Enable the action button when a file is selected
      shinyjs::enable("calcFutilityOneLook")
    }
  })
  
  
  observeEvent(input$calcFutilityOneLook, {
    NRep <- 50
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
}

# Run the Shiny app
shinyApp(ui, server)
