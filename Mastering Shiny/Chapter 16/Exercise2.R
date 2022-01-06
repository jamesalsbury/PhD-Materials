library(shiny)
ui <- fluidPage(
  selectInput("type", "type", c("Normal", "Uniform")),
  actionButton("go", "go"),
  plotOutput("plot")
)

server <- function(input, output, session) {
  chosenplot <- reactiveValues(plot = 0)
  
  observeEvent(input$go, {
    
    if (input$type=="Normal"){
      chosenplot$plot = "rnorm"
    } else{
      chosenplot$plot = "runif"
    }
    
  })
  
  
  output$plot <- renderPlot({
    if (chosenplot$plot=="rnorm"){
      hist(rnorm(100))
    } else if (chosenplot$plot=="runif"){
      hist(runif(100))
    }
  })
}

shinyApp(ui, server)