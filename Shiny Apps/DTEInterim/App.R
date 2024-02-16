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
                   
                   fileInput("samplesFile", "Upload samples"),
                   actionButton("calcAssurance", label  = "Calculate assurance", disabled = T)
                   
                 ), 
                 mainPanel = mainPanel(
                   # Generate plotOutputs dynamically based on the number of delay and HR combinations
                   textOutput("assuranceText")
                 )
               ),
      ),
      
      # Output UI ---------------------------------
      tabPanel("Output", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   
                   # Render the initial trial panel when the app starts
                   uiOutput("init_rule"),
                   
                   actionButton("calcFutility", "Calculate"),
                   
                   actionButton("addRuleButton", "Add rule"),
                   uiOutput("ruleInputs"),
                   actionButton("removeRule", "Remove last rule", disabled = TRUE)
                 ), 
                 mainPanel = mainPanel(
                   textOutput("futilityText")
                 )
               ),
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
      shinyjs::enable("calcAssurance")
    }
  })
  

  
  observeEvent(input$calcAssurance, {
   
    
    treatmentSamplesDF <- data$treatmentSamplesDF
    
    NRep <- 1000
    assvec <- rep(NA, NRep)
    
    
    withProgress(message = 'Calculating assurance', value = 0, {
      
      
    
    for (i in 1:NRep){
      
      #Compute treatment times
      HRStar <- sample(treatmentSamplesDF[,2], 1)
      bigT <- sample(treatmentSamplesDF[,1], 1)
      
      #Simulate control and treatment data
      dataCombined <- SimDTEDataSet(input$numPatients, input$lambdac, bigT, HRStar, input$recTime)  
        
      #Censor these data at numEvents
      dataCombined <- CensFunc(dataCombined, input$numEvents)$dataCombined

      #Do a log-rank test on these data
      test <- survdiff(Surv(survival_time, status)~group, data = dataCombined)
      coxmodel <- coxph(Surv(survival_time, status)~group, data = dataCombined)
      deltad <- as.numeric(exp(coef(coxmodel)))
      assvec[i] <- (test$chisq > qchisq(0.95, 1) & deltad<1)
      
      incProgress(1/NRep)
      
    }
      
    })
    
    
    output$assuranceText <- renderText({
      phat <- mean(assvec)
      LB <- phat - 1.96 * sqrt(phat * (1 - phat) / NRep)
      UB <- phat + 1.96 * sqrt(phat * (1 - phat) / NRep)
      result <- paste0("Assurance is ~ ", round(phat, 3), " (", round(LB, 3), ", ", round(UB, 3), ")")
      HTML(result)
    })
    

    
  })
  
  

  
  # Output Logic ---------------------------------
  
  
  output$init_rule <- renderUI({
    tags$div(
      id = "rule1",
      wellPanel(
        withMathJax(),
        title = "Set 1",
        selectInput("lookType1", 'Look type', choices = c("One look", "Two looks", "Proposed", "Bayesian"), selected = "One look")
      )
    )
  })
  
  
  observeEvent(input$calcFutility, {
    
    NRep <- 100
    assvec <- rep(NA, NRep)
    SSvec <- rep(NA, NRep)
    DurationVec <- rep(NA, NRep)
    Outcome <- ""
    
    for (i in 1:NRep){
      treatmentSamplesDF <- data$treatmentSamplesDF
      
      #Compute treatment times
      HRStar <- sample(treatmentSamplesDF[,2], 1)
      bigT <- sample(treatmentSamplesDF[,1], 1)
      
      #Simulate control and treatment data
      dataCombined <- SimDTEDataSet(input$numPatients, input$lambdac, bigT, HRStar, input$recTime)  
      
      #Perform a futility look at 50% Information Fraction
      
      #Censor these data at numEvents
      futilityDF <- CensFunc(dataCombined, input$numEvents*0.5)
      
      coxmodel <- coxph(Surv(survival_time, status)~group, data = futilityDF$dataCombined)
      deltad <- as.numeric(exp(coef(coxmodel)))
      if (deltad>1){
        Outcome <- "Stop at look 1"
        assvec[i] <- 0
        SSvec[i] <- futilityDF$SS
        DurationVec[i] <- futilityDF$censTime
      } else {
        
      }
      
    }
    
    
    
    print(SSvec)
    
    
    
    
    
    
    output$futilityText <- renderText({
      
      #futilityDF$
      
    })
    
  })
  

  
}

# Run the Shiny app
shinyApp(ui, server)




  
  
  
  