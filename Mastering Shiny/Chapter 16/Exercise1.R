library(shiny)
ui <- fluidPage(
  actionButton("rnorm", "Normal"),
  actionButton("runif", "Uniform"),
  plotOutput("plot")
)

server <- function(input, output, session) {
  
  chosenplot <- reactiveValues(plot = "runif")
  
  observeEvent(input$rnorm, {
    chosenplot$plot = "rnorm"
  })
  
  observeEvent(input$runif, {
    chosenplot$plot = "runif"
  })
  
  output$plot <- renderPlot({
    if (chosenplot$plot=="rnorm"){
      hist(rnorm(100))
    } else{
      hist(runif(100))
    }
  })
}

shinyApp(ui, server)