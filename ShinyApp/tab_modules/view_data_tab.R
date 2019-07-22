# Module for the View Data tab


viewDataUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(
      width = 12,
      uiOutput(ns("data_tabs"))
    )
  )
}

viewData <- function(input, output, session, supplied_data_sets) {
  
  
  observe({
    mapply(function(table_data, table_name){
      output[[table_name]] <- renderDataTable(table_data)
    }, table_data = reactiveValuesToList(supplied_data_sets), table_name = names(reactiveValuesToList(supplied_data_sets)))
  })
  
  
  output$data_tabs <- renderUI({
    ns <- session$ns
    tab_args <- lapply(names(reactiveValuesToList(supplied_data_sets)), function(table_name){
      tabPanel(
        title = table_name, 
        dataTableOutput(ns(table_name))
      )
    })
    
    tab_args$width <- 12
    
    
    do.call(tabBox, tab_args)
  })
}


viewDataTabUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      viewDataUI(id = ns("raw_data_view"))
    )
  )
}

viewDataTab <- function(input, output, session, supplied_data_sets) {
  callModule(viewData, "raw_data_view", supplied_data_sets)
}