library(shiny)
library(shinyjs)
library(ggplot2)

# UI definition
ui <- fluidPage(
  withMathJax(),
  shinyjs::useShinyjs(),
  
  # Application title
  titlePanel("Bayesian Interim Analyses: Delayed Treatment Effects"),
  
  mainPanel(
    tabsetPanel(
      # Data Generating UI ---------------------------------
      tabPanel("Data Generating", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   numericInput("numPatients", 'Number of Patients (in each group, 1:1)', value=340, min=0),
                   numericInput("numEvents", 'Number of Events', value=512, min=0),
                   numericInput("lambdac", '\\( \\lambda_c \\)', value = 0.08),
                   
                   # Render the initial trial panel when the app starts
                   uiOutput("init_trial"),
                   
                   actionButton("addTrialButton", "Add trial"),
                   uiOutput("delayHRInputs"),
                   actionButton("removeTrial", "Remove last trial", disabled = TRUE)
                 ), 
                 mainPanel = mainPanel(
                   # Generate plotOutputs dynamically based on the number of delay and HR combinations
                   uiOutput("DGPlots")
                 )
               ),
      ),
      
      # Elicited Distributions UI ---------------------------------
      tabPanel("Elicited Distributions", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   wellPanel(
                     selectInput("TDist", 'T', choices = c("Normal", "Student-t", "Gamma", "Log normal", 
                                                           "Log student-t", "Beta", "Mirror gamma", "Mirror log normal",
                                                           "Mirror log Student-t"), selected = "Normal"),
                     uiOutput("TDistParams") # Output for TDist parameters
                   ),
                   wellPanel(
                     selectInput("HRStarDist", label = HTML(' \\(\\text{HR}^*\\)'), choices = c("Normal", "Student-t", "Gamma", "Log normal", 
                                                                                                "Log student-t", "Beta", "Mirror gamma", "Mirror log normal",
                                                                                                "Mirror log Student-t"), selected = "Normal"),
                     uiOutput("HRStarDistParams") # Output for HRStarDist parameters
                   )
                 ), 
                 mainPanel = mainPanel(
                   plotOutput("TDistPlot"),
                   plotOutput("HRStarDistPlot")
                 )
               )
      ),
      
      
      # Rules UI ---------------------------------
      tabPanel("Rules", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   
                   # Render the initial trial panel when the app starts
                   uiOutput("init_rule"),
                   
                   actionButton("addRuleButton", "Add rule"),
                   uiOutput("ruleInputs"),
                   actionButton("removeRule", "Remove last rule", disabled = TRUE)
                 ), 
                 mainPanel = mainPanel(
                   
                 )
               ),
      ),
      
      # Output UI ---------------------------------
      tabPanel("Output", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   
                    textOutput("numDG"),
                    textOutput("numRules"),
                   actionButton("Compute", "Compute", disabled = TRUE)
                 ), 
                 mainPanel = mainPanel(
                   
                 )
               ),
      )
      
    ), style='width: 1000px; height: 600px'
  )
)

# Initialize initial trial panel
init_trial <- tags$div(
  id = "trial1",
  wellPanel(
    title = "Set 1",
    numericInput(inputId = "delay1", label = "Delay", value = 0),
    numericInput(inputId = "HRStar1", label = HTML(' \\(\\text{HR}^*\\)'),  value = 1),
    numericInput(inputId = "Rec1", label = "Recruitment length", value = 34)
  )
)





# Server logic
server <- function(input, output, session) {
  
  # Initialize reactive values
  rv <- reactiveValues(trial_input_count = 1,
                       rule_input_count = 1)
  
  
  # Data Generating Logic ---------------------------------
  
  observeEvent(input$addTrialButton, {
    rv$trial_input_count <- rv$trial_input_count + 1
    insertUI(
      selector = "#addTrialButton",
      where = "beforeBegin",
      ui = tags$div(
        id = paste0("trial", rv$trial_input_count),
        wellPanel(
          withMathJax(),
          title = paste("Set", rv$trial_input_count),
          numericInput(inputId = paste0("delay", rv$trial_input_count), label = "Delay", value = 0),
          numericInput(inputId = paste0("HRStar", rv$trial_input_count), label = HTML(' \\(\\text{HR}^*\\)'),  value = 1),
          numericInput(inputId = paste0("Rec", rv$trial_input_count), label = "Recruitment length", value = 34)
        )
      )
    )
    if (rv$trial_input_count > 1) {
      shinyjs::enable("removeTrial")
    }
  })
  

  
  
  # Function to handle removing inputs
  observeEvent(input$removeTrial, {
    if (rv$trial_input_count > 1) {
      removeUI(selector = paste0("#trial", rv$trial_input_count))
      rv$trial_input_count <- rv$trial_input_count - 1
      if (rv$trial_input_count == 1) {
        shinyjs::disable("removeTrial")
      }
    }
  })
  
  # Render the initial trial panel when the app starts
  output$init_trial <- renderUI({
    tags$div(
      id = "trial1",
      wellPanel(
        withMathJax(),
        title = "Set 1",
        numericInput(inputId = "delay1", label = "Delay", value = 0),
        numericInput(inputId = "HRStar1", label = HTML(' \\(\\text{HR}^*\\)'), value = 1),
        numericInput(inputId = "Rec1", label = "Recruitment length", value = 34)
      )
    )
  })
  
  # Render the plots
  output$DGPlots <- renderUI({
    plot_output_list <- lapply(1:rv$trial_input_count, function(i) {
      plotOutput(outputId = paste0("plot", i))
    })
    do.call(tagList, plot_output_list)
  })
  
  # Reactive expression for rendering plots
  observe({
    lapply(1:rv$trial_input_count, function(i) {
      output[[paste0("plot", i)]] <- renderPlot({
        # Extract delay and HRStar inputs
        delay <- input[[paste0("delay", i)]]
        HRStar <- input[[paste0("HRStar", i)]]
        
        # Generate survival curves
        Time <- seq(0, 1000, by = 0.1)
        ControlSurv <- exp(-input$lambdac * Time)
        TreatmentSurv <- ifelse(Time <= delay, exp(-input$lambdac * Time), exp(-input$lambdac * delay - input$lambdac * HRStar * (Time - delay)))
        
        # Determine max time for plotting
        MaxTimeC <- which(ControlSurv < 0.001)[1]
        MaxTimeT <- which(TreatmentSurv < 0.001)[1]
        MaxTime <- max(MaxTimeC, MaxTimeT)
        
        # Plot survival curves
        plot(Time[1:MaxTime], ControlSurv[1:MaxTime], ylim = c(0, 1), type = "l", col = "blue", xlab = "Time (months)", ylab = "Survival")
        lines(Time[1:MaxTime], TreatmentSurv[1:MaxTime], col = "red")
        
        # Add legend
        legend("topright", legend = c("Control", "Treatment"), col = c("blue", "red"), lty = 1)
      })
    })
  })
  
  
  
  # Elicited Distributions Logic ---------------------------------
  
  # Define UI outputs for parameters
  output$TDistParams <- renderUI({
    dist <- input$TDist
    
    if (dist %in% c("Normal", "Log normal", "Mirror log normal")) {
      list(
        numericInput("TDistMean", label = "Mean:", value = 0),
        numericInput("TDistSD", label = "Standard Deviation:", value = 1)
      )
    } else if (dist %in% c("Student-t", "Log student-t", "Mirror log Student-t")) {
      list(
        numericInput("TDistMean", label = "Mean:", value = 0),
        numericInput("TDistDF", label = "Degrees of Freedom:", value = 1)
      )
    } else if (dist %in% c("Gamma", "Beta", "Mirror gamma")) {
      list(
        numericInput("TDistShape", label = "Shape:", value = 1),
        numericInput("TDistRate", label = "Rate:", value = 1)
      )
    } 
    # Add conditions for other distributions
  })
  
  output$HRStarDistParams <- renderUI({
    dist <- input$HRStarDist
    
    if (dist %in% c("Normal", "Log normal", "Mirror log normal")) {
      list(
        numericInput("HRStarDistMean", label = "Mean:", value = 0),
        numericInput("HRStarDistSD", label = "Standard Deviation:", value = 1)
      )
    } else if (dist %in% c("Student-t", "Log student-t", "Mirror log Student-t")) {
      list(
        numericInput("HRStarDistMean", label = "Mean:", value = 0),
        numericInput("HRStarDistDF", label = "Degrees of Freedom:", value = 1)
      )
    } else if (dist %in% c("Gamma", "Beta", "Mirror gamma")) {
      list(
        numericInput("HRStarDistShape", label = "Shape:", value = 1),
        numericInput("HRStarDistRate", label = "Rate:", value = 1)
      )
    } 
    # Add conditions for other distributions
  })

  
  # Define reactive values to store distribution parameters
  distParams <- reactiveValues()
  
  observe({
    # Store TDist parameters
    if (input$TDist == "Normal") {
      distParams$TDist <- list(mean = input$TDistMean, sd = input$TDistSD)
    } else if (input$TDist == "Student-t") {
      distParams$TDist <- list(mean = input$TDistMean, df = input$TDistDF)
    } else if (input$TDist %in% c("Gamma", "Beta")) {
      distParams$TDist <- list(shape = input$TDistShape, rate = input$TDistRate)
    }
    
    # Store HRStarDist parameters
    if (input$HRStarDist == "Normal") {
      distParams$HRStarDist <- list(mean = input$HRStarDistMean, sd = input$HRStarDistSD)
    } else if (input$HRStarDist == "Student-t") {
      distParams$HRStarDist <- list(mean = input$HRStarDistMean, df = input$HRStarDistDF)
    } else if (input$HRStarDist %in% c("Gamma", "Beta")) {
      distParams$HRStarDist <- list(shape = input$HRStarDistShape, rate = input$HRStarDistRate)
    }
  })
  
  output$TDistPlot <- renderPlot({
    dist_name <- input$TDist
    dist_params <- distParams$TDist
    
    # Define the x-axis range dynamically based on the distribution parameters
    switch(dist_name,
           "Normal" = {
             x_range <- c(dist_params$mean - 4 * dist_params$sd, dist_params$mean + 4 * dist_params$sd)
             x <- seq(x_range[1], x_range[2], length.out = 1000)
           },
           "Student-t" = {
             x_range <- c(-10, 10) # Adjust according to your requirement
             x <- seq(x_range[1], x_range[2], length.out = 1000)
           },
           "Gamma" = {
             x_range <- c(0, 10) # Adjust according to your requirement
             x <- seq(x_range[1], x_range[2], length.out = 1000)
           },
           "Log normal" = {
             x_range <- c(exp(dist_params$mean - 4 * dist_params$sd), exp(dist_params$mean + 4 * dist_params$sd))
             x <- seq(x_range[1], x_range[2], length.out = 1000)
           },
           "Log student-t" = {
             x_range <- c(exp(-10), exp(10)) # Adjust according to your requirement
             x <- seq(x_range[1], x_range[2], length.out = 1000)
           },
           "Beta" = {
             x_range <- c(0, 1) # Adjust according to your requirement
             x <- seq(x_range[1], x_range[2], length.out = 1000)
           },
           "Mirror gamma" = {
             x_range <- c(0, 10) # Adjust according to your requirement
             x <- seq(x_range[1], x_range[2], length.out = 1000)
           },
           "Mirror log normal" = {
             x_range <- c(exp(dist_params$mean - 4 * dist_params$sd), exp(dist_params$mean + 4 * dist_params$sd))
             x <- seq(x_range[1], x_range[2], length.out = 1000)
           },
           "Mirror log Student-t" = {
             x_range <- c(exp(-10), exp(10)) # Adjust according to your requirement
             x <- seq(x_range[1], x_range[2], length.out = 1000)
           }
    )
    
    # Calculate y values based on the selected distribution and parameters
    y <- switch(dist_name,
                "Normal" = dnorm(x, mean = dist_params$mean, sd = dist_params$sd),
                "Student-t" = dt(x, df = dist_params$df, ncp = dist_params$mean),
                "Gamma" = dgamma(x, shape = dist_params$shape, rate = dist_params$rate),
                "Log normal" = dlnorm(x, meanlog = dist_params$mean, sdlog = dist_params$sd),
                "Log student-t" = dt(exp(x), df = dist_params$df, ncp = dist_params$mean) * exp(x),
                "Beta" = dbeta(x, shape1 = dist_params$shape, shape2 = dist_params$rate),
                "Mirror gamma" = dgamma(abs(x), shape = dist_params$shape, rate = dist_params$rate),
                "Mirror log normal" = dlnorm(abs(x), meanlog = dist_params$mean, sdlog = dist_params$sd),
                "Mirror log Student-t" = dt(exp(abs(x)), df = dist_params$df, ncp = dist_params$mean) * exp(abs(x))
    )
    
    # Plot the density plot
    ggplot(data.frame(x = x, y = y), aes(x, y)) +
      geom_line(color = "blue") +
      labs(title = paste("Density Plot of", dist_name), x = "Value", y = "Density") +
      xlim(x_range) # Set x-axis limits dynamically
  })
  
  
  
  # Plot HRStarDist density
  output$HRStarDistPlot <- renderPlot({
    x <- seq(-10, 10, length.out = 1000) # Generate x values for the density plot
    y <- switch(
      input$HRStarDist,
      "Normal" = dnorm(x, mean = distParams$HRStarDist$mean, sd = distParams$HRStarDist$sd),
      "Student-t" = dt(x, df = distParams$HRStarDist$df, ncp = distParams$HRStarDist$mean),
      "Gamma" = dgamma(x, shape = distParams$HRStarDist$shape, rate = distParams$HRStarDist$rate),
      "Log normal" = dlnorm(x, meanlog = distParams$HRStarDist$mean, sdlog = distParams$HRStarDist$sd),
      "Log student-t" = dt(exp(x), df = distParams$HRStarDist$df, ncp = distParams$HRStarDist$mean) * exp(x),
      "Beta" = dbeta(x, shape1 = distParams$HRStarDist$shape, shape2 = distParams$HRStarDist$rate),
      "Mirror gamma" = dgamma(abs(x), shape = distParams$HRStarDist$shape, rate = distParams$HRStarDist$rate),
      "Mirror log normal" = dlnorm(abs(x), meanlog = distParams$HRStarDist$mean, sdlog = distParams$HRStarDist$sd),
      "Mirror log Student-t" = dt(exp(abs(x)), df = distParams$HRStarDist$df, ncp = distParams$HRStarDist$mean) * exp(abs(x))
    )
    
    ggplot(data.frame(x = x, y = y), aes(x, y)) +
      geom_line(color = "red") +
      labs(title = paste("Density Plot of", input$HRStarDist), x = "Value", y = "Density")
  })
  
  
  # Rules Logic ---------------------------------
  
  observeEvent(input$addRuleButton, {
    rv$rule_input_count <- rv$rule_input_count + 1
    insertUI(
      selector = "#addRuleButton",
      where = "beforeBegin",
      ui = tags$div(
        id = paste0("rule", rv$rule_input_count),
        wellPanel(
          withMathJax(),
          title = paste("Set", rv$rule_input_count),
          selectInput(inputId = paste0("lookType", rv$rule_input_count), 'Look type', choices = c("One look", "Two looks", "Proposed", "Bayesian"), selected = "One look"),
          numericInput(inputId = paste0("numInput1", rv$rule_input_count), label = "At Information Fraction (0-1)", value = 0),
          numericInput(inputId = paste0("numInput2", rv$rule_input_count), label = "Stop if observed HR >", value = 0),
          numericInput(inputId = paste0("numInput3", rv$rule_input_count), label = "At Information Fraction (0-1)", value = 0),
          numericInput(inputId = paste0("numInput4", rv$rule_input_count), label = "Stop if observed HR >", value = 0),
          numericInput(inputId = paste0("numInput5", rv$rule_input_count), label = "At Information Fraction (0-1)", value = 0),
          numericInput(inputId = paste0("numInput6", rv$rule_input_count), label = "Stop if observed HR >", value = 0),
          numericInput(inputId = paste0("numInput7", rv$rule_input_count), label = "At Information Fraction (0-1)", value = 0),
          numericInput(inputId = paste0("numInput8", rv$rule_input_count), label = "Stop if observed HR >", value = 0)
        )
      )
    )
    if (rv$rule_input_count > 1) {
      shinyjs::enable("removeRule")
    }
  })
  
  
  
  
  # Function to handle removing inputs
  observeEvent(input$removeRule, {
    if (rv$rule_input_count > 1) {
      removeUI(selector = paste0("#rule", rv$rule_input_count))
      rv$rule_input_count <- rv$rule_input_count - 1
      if (rv$rule_input_count == 1) {
        shinyjs::disable("removeRule")
      }
    }
  })
  
  # Render the initial rule panel when the app starts
  output$init_rule <- renderUI({
    tags$div(
      id = "rule1",
      wellPanel(
        withMathJax(),
        title = "Set 1",
        selectInput("lookType1", 'Look type', choices = c("One look", "Two looks", "Proposed", "Bayesian"), selected = "One look"),
        numericInput("numInput11", label = "At Information Fraction (0-1)", value = 0),
        numericInput("numInput21", label = "Stop if observed HR >", value = 0),
        numericInput("numInput31", label = "At Information Fraction (0-1)", value = 0),
        numericInput("numInput41", label = "Stop if observed HR >", value = 0),
        numericInput("numInput51", label = "At Information Fraction (0-1)", value = 0),
        numericInput("numInput61", label = "Stop if observed HR >", value = 0),
        numericInput("numInput71", label = "At Information Fraction (0-1)", value = 0),
        numericInput("numInput81", label = "Stop if observed HR >", value = 0)
      )
    )
  })
  
  
  # Observe changes in lookType and adapt UI
  observeEvent(seq_len(rv$rule_input_count), {
    lapply(1:rv$rule_input_count, function(i) {
      observeEvent(input[[paste0("lookType", i)]], {
        shinyjs::hide(paste0("numInput1", i))
        shinyjs::hide(paste0("numInput2", i))
        shinyjs::hide(paste0("numInput3", i))
        shinyjs::hide(paste0("numInput4", i))
        shinyjs::hide(paste0("numInput5", i))
        shinyjs::hide(paste0("numInput6", i))
        shinyjs::hide(paste0("numInput7", i))
        shinyjs::hide(paste0("numInput8", i))
        if (!is.null(input[[paste0("lookType", i)]])) {
          if (input[[paste0("lookType", i)]] == "One look") {
            shinyjs::show(paste0("numInput1", i))
            shinyjs::show(paste0("numInput2", i))
            updateNumericInput(session, paste0("numInput1", i), label = "At Information Fraction (0-1)")
            updateNumericInput(session, paste0("numInput2", i), label = "Stop if observed HR >")
          } else if (input[[paste0("lookType", i)]] == "Two looks") {
            shinyjs::show(paste0("numInput1", i))
            shinyjs::show(paste0("numInput2", i))
            shinyjs::show(paste0("numInput3", i))
            shinyjs::show(paste0("numInput4", i))
            updateNumericInput(session, paste0("numInput1", i), label = "At Information Fraction (0-1)")
            updateNumericInput(session, paste0("numInput2", i), label = "Stop if observed HR >")
            updateNumericInput(session, paste0("numInput3", i), label = "At Information Fraction (0-1)")
            updateNumericInput(session, paste0("numInput4", i), label = "Stop if observed HR >")
          } else if (input[[paste0("lookType", i)]] == "Proposed") {
            shinyjs::show(paste0("numInput1", i))
            shinyjs::show(paste0("numInput2", i))
            shinyjs::show(paste0("numInput3", i))
            shinyjs::show(paste0("numInput4", i))
            shinyjs::show(paste0("numInput5", i))
            shinyjs::show(paste0("numInput6", i))
            shinyjs::show(paste0("numInput7", i))
            shinyjs::show(paste0("numInput8", i))
            updateNumericInput(session, paste0("numInput1", i), label = "At Information Fraction (0-1)")
            updateNumericInput(session, paste0("numInput2", i), label = "AND more than ")
            updateNumericInput(session, paste0("numInput3", i), label = "events occur after")
            updateNumericInput(session, paste0("numInput4", i), label = "Stop if observed HR >")
            updateNumericInput(session, paste0("numInput5", i), label = "At Information Fraction (0-1)")
            updateNumericInput(session, paste0("numInput6", i), label = "AND more than ")
            updateNumericInput(session, paste0("numInput7", i), label = "events occur after")
            updateNumericInput(session, paste0("numInput8", i), label = "Stop if observed HR >")
          } else if (input[[paste0("lookType", i)]] == "Bayesian") {
            shinyjs::show(paste0("numInput1", i))
            shinyjs::show(paste0("numInput2", i))
            shinyjs::show(paste0("numInput3", i))
            shinyjs::show(paste0("numInput4", i))
            updateNumericInput(session, paste0("numInput1", i), label = "At Information Fraction (0-1)")
            updateNumericInput(session, paste0("numInput2", i), label = "Stop if BPP <")
            updateNumericInput(session, paste0("numInput3", i), label = "At Information Fraction (0-1)")
            updateNumericInput(session, paste0("numInput4", i), label = "Stop if BPP <")
          }
        }
      })
    })
  })
  

  # Output Logic ---------------------------------
  
  output$numDG <- renderText({
    if (rv$trial_input_count==1){
      paste0("You have chosen ", rv$trial_input_count, " data generating mechaniam")
    } else {
      paste0("You have chosen ", rv$trial_input_count, " data generating mechaniams")
    }
  })
  
  output$numRules <- renderText({
    if (rv$rule_input_count==1){
      paste0("You have chosen ", rv$rule_input_count, " rule")
    } else {
      paste0("You have chosen ", rv$rule_input_count, " rules")
    } 
  })
  
  
  
  
  
  
}

# Run the Shiny app
shinyApp(ui, server)




  
  
  
  