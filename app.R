library(tidyverse)
library(shiny)

ui <- fluidPage(
  
  titlePanel('Metabolomics results app'),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = 'response', label = 'Response', choices = c('AD', 'PD', 'Age')),
      selectInput(inputId = 'data', label = 'Dataset', choices = c('GOT', 'Lipids', 'Both'))
    ),
    mainPanel(
      imageOutput(outputId = 'plot'),
      textOutput(outputId = 'filename')
    )
  )
)
















server <- function(input, output){
  plot_path <- 'plots'
  output$filename <- reactive({
    name <- character(0)
    if(input$response == 'Age'){
      name <- 'age'
    }
    paste(input$response, input$data, sep = '_')
  })
  
  output$plot <- renderImage({
    list(src = 'plots/age_pred_resid_all_insample_4.png')
  })
  
  
  
}
shinyApp(ui = ui, server = server)
