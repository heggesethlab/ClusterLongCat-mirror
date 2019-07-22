# Module for the Upload Data tab

uploadDataTabUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(title = "Upload data", status = "warning", solidHeader = TRUE,
        width = 12,
        includeMarkdown("tab_text/upload_instructions/upload_instructions.md"),

        column(
          width = 12,
          fluidRow(

            align = "center",
            dataTableOutput(ns("example_table")),
            br(),
            br(),
            h3("Choose Dataset"),
            fileInput(
              ns("upload_input"),
              "",
              placeholder = "example.csv",
              width = "50%",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")
            ),

            uiOutput(ns("uploaded_data"))
          )
        )
      )
    )
  )
}

uploadDataTab <- function(input, output, session) {
  output$example_table <- renderDataTable({
    tibble(
      ID =         c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4),
      state =      c(1, 2, 1, 2, 3, 1, 3, 2, 1, 2, 3, 2, 3, 1),
      duration =   c(5, 1, 2, 3, 1, 6, 6, 2, 1, 2, 2, 3, 2, 2),
      start_time = c(0, 5, 6, 0, 3, 4, 0, 6, 8, 0, 2, 4, 7, 9)
    )
  })

  new_data <- reactive({
    in_file <- input$upload_input
    tried_to_upload <- TRUE
    if (is.null(in_file))
      return(NULL)

    new_table <- read.csv(in_file$datapath) %>%
      select(-X)

    #check its valid
    if(sum(!(c("state", "ID", "duration", "start_time") %in% names(new_table)))>0){
      return(NULL)
    }

    if(!is.null(new_table)){
      return(list(data = new_table, name = gsub("\\..{3}", "", in_file$name)))
    }

  })

  observeEvent(input$upload_input, {
    print(class(new_data()))
    if(!is.null(new_data()) ) {
      output$contents <- renderDataTable({
        new_data()$data
      })

      output$uploaded_data <- renderUI({
        ns <- session$ns

        return(tagList(
          h3("Uploaded Data"),
          dataTableOutput(ns("contents")))
        )
      })
    } else {
      output$uploaded_data <- renderUI({
        h3("Something went wrong, please check your file and try again.")
      })
    }
  })

  return(list(name = reactive(new_data()$name), data = reactive(new_data()$data)))
}
