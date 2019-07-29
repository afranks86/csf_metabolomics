library(tidyverse)
library(shiny)

ui <- fluidPage(
  
  titlePanel('Metabolomics results app'),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = 'response', label = 'Response', choices = c('AD', 'PD', 'Age')),
      selectInput(inputId = 'data', label = 'Dataset', choices = c('GOT', 'Lipids', 'Both'))
    )
  )
)
















server <- function(input, output){
  
  
  
  
  
  
}
shinyApp(ui = ui, server = server)
