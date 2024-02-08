library(shiny)
library(shinyjs)

# UI definition
ui <- fluidPage(
  withMathJax(),
  shinyjs::useShinyjs(),
  
  # Application title
  titlePanel("Assurance: Delayed Treatment Effects"),
  
  mainPanel(
    
    tabsetPanel(
      # Control UI ---------------------------------
      
      
      tabPanel("Data Generating", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   numericInput("numPatients", 'Number of Patients (in each group, 1:1)', value=340, min=0),
                   numericInput("numEvents", 'Number of Events', value=512, min=0),
                   numericInput("lambdac", '\\( \\lambda_c \\)', value = 0.08),
                   actionButton("addButton", "Add trial"),
                   uiOutput("delayHRInputs"),
                   actionButton("removeInput", "Remove last trial", disabled = TRUE)
                 ), 
                 mainPanel = mainPanel(
                   # Generate plotOutputs dynamically based on the number of delay and HR combinations
                   uiOutput("plots")
                 )
               ),
      ),
      
      
      tabPanel("Futility rules", 
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   radioButtons("uploadSampleCheck", "Do you wish to upload a MCMC sample?", choices = c("Yes", "No"), selected = "No"),
                   numericInput("lambdacmean", 'mean (\\(\\text{scale} = \\lambda_c \\))', value=0.08, min=0),
                   numericInput("gammacmean", 'mean \\(\\text{shape} = \\gamma_c \\)', value=0.8, min=0)
                   
                   
                 ), 
                 mainPanel = mainPanel(
                   plotOutput("plotControl"),
                   htmlOutput("recommendedParams")
                 )
               ),
      ),
      
      
      
      #Help UI ---------------------------------
      #Need to say what files can be uploaded in the control sample
      #Link to SHELF for the elicitation
      tabPanel("Help",
               HTML("<p>This app implements the method as outlined in this <a href='https://jamesalsbury.github.io/'>paper</a>. For every tab, there is a brief summary below, for any other questions, comments or
                  feedback, please contact <a href='mailto:jsalsbury1@sheffield.ac.uk'>James Salsbury</a>.</p>"),
               HTML("<p><u>Control</u></p>"),
               HTML("<p>Here, the parameters for the control survival are specified, the parameterisation for the control survival curve is:</p>"),
               withMathJax(paste0(" $$S_c(t) = \\text{exp}\\{-(\\lambda_ct)^{\\gamma_c}\\}$$")),
               HTML("<p>The control tab also allows an upload of an Excel file containing survival data for the control sample. 
                  The Excel file needs to have two columns: the first column containing the survival time, and the second column containing the event status (1 for dead, 0 for alive).</p>"),
               HTML("<p><u>Eliciting the two parameters: T and post-delay HR</u></p>"),
               HTML("<p>These two tabs elicit beliefs from the user about two quantities: the length of delay, T, and the post-delay HR. The parameterisation for the treatment survival curve is:</p>"),
               withMathJax(paste0(" $$S_t(t) = \\text{exp}\\{-(\\lambda_ct)^{\\gamma_c}\\}, t \\leq T$$")),
               withMathJax(paste0(" $$S_t(t) = \\text{exp}\\{-(\\lambda_cT)^{\\gamma_c}-\\lambda_t^{\\gamma_t}(t^{\\gamma+t}-T^{\\gamma_t})\\}, t > T$$")),
               HTML("<p>The elicitation technique is based on SHELF, more guidance can be found <a href='https://shelf.sites.sheffield.ac.uk/'>here.</a></p>"),
               HTML("<p>There is also an option to include some mass at T = 0 and HR = 1.</p>"),
               HTML("<p><u>Feedback</u></p>"),
               HTML("<p>The plot shows the control survival curve (from the control tab), along with the median elicited treatment line (calculated from the previous two tabs). 
                  There are also three optional quantities to view: the first is the median survival time for both groups, the second is a 95% confidence interval for T and the
                  third shows a 80% confidence interval for the treatment curve. When the median survival time is added to the plot, a second plot is shown below. Initially, this plot is a histogram for the 
                  treatment median survival time. However, the user is able to change this median to any other quantile of interest.</p>"),
               HTML("<p><u>Calculating assurance</u></p>"),
               HTML("<p>This tab allows the user to calculate assurance, given the control parameters and elicited prior distributions  for T and post-delay HR. Some additional questions
                  about the trial are found on the left-hand panel. The app assumes uniform recruitment and uses a log-rank test for analysis of the simulated data. Depending on your processor speed, this 
                  calculation can take between 30-40 seconds. Once the calculation is complete, the app shows two plots. The top plot shows assurance (along with standard error curves) and target effect curve - which
                  shows the proportion of trials in which the target effect was observed. The bottom plot shows the average hazard ratio observed and the corresponding confidence intervals for this. </p>")
      ),
      
      
    ), style='width: 1000px; height: 600px',
    
  )
)


# Server logic
server <- function(input, output, session) {
  # Initialize reactive values
  rv <- reactiveValues(input_count = 1,
                       delays = c(3),
                       HRs = c(1),
                       recruit_lengths = c(34))
  
  # Function to render delay and HRstar inputs
  output$delayHRInputs <- renderUI({
    lapply(1:rv$input_count, function(i) {
      wellPanel(
        title = paste("Set", i),
        numericInput(inputId = paste0("delay", i), label = "Delay", value = rv$delays[i]),
        numericInput(inputId = paste0("HRStar", i), label = "HR", value = rv$HRs[i]),
        numericInput(inputId = paste0("Rec", i), label = "Recruitment length", value = rv$recruit_lengths[i])
      )
    })
  })
  
  # Function to handle adding inputs
  observeEvent(input$addButton, {
    rv$input_count <- rv$input_count + 1
    rv$delays <- c(rv$delays, 3)  # Initialize new trial with default values
    rv$HRs <- c(rv$HRs, 1)
    rv$recruit_lengths <- c(rv$recruit_lengths, 34)
    if (rv$input_count > 1) {
      shinyjs::enable("removeInput")
    }
  })
  
  # Function to handle removing inputs
  observeEvent(input$removeInput, {
    if (rv$input_count > 0) {
      rv$input_count <- rv$input_count - 1
      rv$delays <- head(rv$delays, -1)  # Remove last trial
      rv$HRs <- head(rv$HRs, -1)
      rv$recruit_lengths <- head(rv$recruit_lengths, -1)
    }
    if (rv$input_count == 1) {
      shinyjs::disable("removeInput")
    }
  })
  
  # Create separate reactive expressions for each plot
  plot_reactive <- reactiveValues()
  observe({
    lapply(1:rv$input_count, function(i) {
      plot_reactive[[paste0("plot", i)]] <- renderPlot({
        # Capture the delay and HR inputs
        delay <- rv$delays[i]
        HRStar <- rv$HRs[i]
        
        ControlTime <- seq(0, 1000, by = 0.1)
        ControlSurv <- exp(-input$lambdac*ControlTime)
        MaxTime <- which(ControlSurv < 0.001)[1]
        
        plot(ControlTime[1:MaxTime], ControlSurv[1:MaxTime], ylim = c(0,1), type = "l", col = "blue", xlab = "Time (months)", ylab = "Survival")
        
        TreatmentTime1 <- seq(0, delay, by = 0.1)
        TreatmentSurv1 <- exp(-input$lambdac*TreatmentTime1)
        
        TreatmentTime2 <- seq(delay, 1000, by = 0.1)
        TreatmentSurv2 <- exp(-delay*input$lambdac-input$lambdac*HRStar*(TreatmentTime2-delay))
        
        lines(TreatmentTime1, TreatmentSurv1, col = "red")
        lines(TreatmentTime2, TreatmentSurv2, col = "red")
      })
    })
  })
  
  # Render the plots dynamically
  output$plots <- renderUI({
    plot_output_list <- lapply(1:rv$input_count, function(i) {
      plot_reactive[[paste0("plot", i)]]
    })
    do.call(tagList, plot_output_list)
  })
}

# Run the Shiny app
shinyApp(ui, server)






  
  
  
  