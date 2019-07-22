# Module for the about us tab

referencesTabUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(title = "References", status = "warning", solidHeader = TRUE,
          width = 12,
          uiOutput(ns('citations'))
          
      )
    )
  )
}

referencesTab <- function(input, output, session) {
  output$citations <- renderUI({
    ns <- session$ns
    rmarkdown::render(input = "tab_text/references/citations.Rmd",
                      output_format = rmarkdown::html_fragment(self_contained = TRUE),
                      output_file = 'citations.html')  
    shiny::includeHTML('tab_text/references/citations.html') 
}) 
}