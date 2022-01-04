
    ui <- shinyUI(fluidPage(
      
      # Application title
      titlePanel("SHELF: bivariate distribution"),
      
      # sidebarLayout(
      mainPanel(tags$style(type="text/css",
                           ".shiny-output-error { visibility: hidden; }",
                           ".shiny-output-error:before { visibility: hidden; }"
      ),
      
      tabsetPanel(
        tabPanel("Parameter 1",
                 fluidRow(
                   column(4, 
                          textInput("limits1", label = h5("Parameter 1 limits"), 
                                    value = "0, 100")
                   ),
                   column(4,
                          textInput("values1", label = h5("Parameter 1 values"), 
                                    value = "25, 50, 75")
                   ),
                   column(4,
                          textInput("probs1", label = h5("Cumulative probabilities"), 
                                    value = "0.25, 0.5, 0.75")
                   )
                 ),
                 fluidRow(
                   column(4, 
                          selectInput("dist1", label = h5("Distribution"), 
                                      choices =  list(Histogram = "hist",
                                                      Normal = "normal", 
                                                      'Student-t' = "t",
                                                      Gamma = "gamma",
                                                      'Log normal' = "lognormal",
                                                      'Log Student-t' = "logt",
                                                      Beta = "beta",
                                                      'Mirror gamma' = "mirrorgamma",
                                                      'Mirror log normal' = "mirrorlognormal",
                                                      'Mirror log Student-t' = "mirrorlogt",
                                                      'Best fitting' = "best"),
                                      #choiceValues = 1:8,
                                      selected = 1
                          )),
                   column(4,conditionalPanel(
                     condition = "input.dist1 == 't' || input.dist1 == 'logt' || input.dist1 == 'mirrorlogt'",
                     numericInput("tdf1", label = h5("Student-t degrees of freedom"),
                                  value = 3)
                   )
                   )
                   
                 ),
                 
                 
                 
                 
                 plotOutput("distPlot1")
                 #tableOutput("valuesPDF1")
        ),
        tabPanel("Parameter 2",
                 fluidRow(
                   column(4, 
                          textInput("limits2", label = h5("Parameter limits"), 
                                    value = "0, 200")
                   ),
                   column(4,
                          textInput("values2", label = h5("Parameter values"), 
                                    value = "30, 40, 60")
                   ),
                   column(4,
                          textInput("probs2", label = h5("Cumulative probabilities"), 
                                    value = "0.25, 0.5, 0.75")
                   )
                 ),
                 fluidRow(
                   column(4, 
                          selectInput("dist2", label = h5("Distribution"), 
                                      choices =  list(Histogram = "hist",
                                                      Normal = "normal", 
                                                      'Student-t' = "t",
                                                      Gamma = "gamma",
                                                      'Log normal' = "lognormal",
                                                      'Log Student-t' = "logt",
                                                      Beta = "beta",
                                                      'Mirror gamma' = "mirrorgamma",
                                                      'Mirror log normal' = "mirrorlognormal",
                                                      'Mirror log Student-t' = "mirrorlogt",
                                                      'Best fitting' = "best"),
                                      #choiceValues = 1:8,
                                      selected = 1
                          )),
                   column(4,
                          conditionalPanel(
                            condition = "input.dist2 == 't' || input.dist2 == 'logt' || input.dist1 == 'mirrorlogt'",
                            numericInput("tdf2", label = h5("degrees of freedom"),
                                         value = 3)
                            
                            
                          ))
                   
                 ),
                 
                 plotOutput("distPlot2")
                 # tableOutput("valuesPDF2")
        ),
        
        tabPanel("Joint distribution",
                 fluidRow(
                   column(4, 
                          numericInput("concProb", h5("Concordance probability"),
                                       value = 0.5,
                                       min = 0, max = 1)
                   ),
                   column(4,
                          numericInput("sampleSize", h5("Sample size"),
                                       value = 1000,
                                       min = 1)
                   ),
                   column(4, 
                          checkboxInput("showSample", "Show sample")
                   )
                 ),
                 plotOutput("bivariatePlot")
        ),
        tabPanel("Help",
                 includeHTML(system.file("shinyAppFiles", "helpBivariate.html",
                                         package="SHELF"))
        )
        
      ),
      wellPanel(
        fluidRow(
          column(3, selectInput("outFormat", label = "Report format",
                                choices = list('html' = "html_document",
                                               'pdf' = "pdf_document",
                                               'Word' = "word_document"))
          ),
          column(3, offset = 1, 
                 numericInput("fs", label = "Font size", value = 12)
          )),
        fluidRow(
          column(3, downloadButton("report", "Download report")
          ),
          column(3, downloadButton("downloadData", "Download sample")
          ),
          column(3, actionButton("exit", "Quit")
          )
        )
        
      )
      
      )
    )
    )
    
    
    server <- function(input, output, session) {
    }
    
    
    shinyApp(ui, server)
    
    