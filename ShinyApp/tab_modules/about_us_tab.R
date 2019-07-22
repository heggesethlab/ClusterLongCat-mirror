# Module for the about us tab

aboutUsTabUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(title = "Dr. Brianna Heggeseth", status = "warning", solidHeader = TRUE,
          width = 12, column(2, tags$img(height = 150, width = 115, src = "brianna_pic.png")), column(1,""),column(9, BriBio)

      )
    ),
    fluidRow(
      box(title = "Ellen Graham", status = "danger", solidHeader = TRUE,
          width = 4,
          fluidRow(column(4, tags$img(height = 130, width = 100, src = "ellen_pic.jpg"))),
          fluidRow(br(),column(12, EllenBio))
      ),
      box(title = "Kieu-Giang Nguyen", status = "danger", solidHeader = TRUE,
          width = 4,
          fluidRow(column(4, tags$img(height = 130, width = 100, src = "kg_pic.jpg"))),
          fluidRow(br(),column(width = 12, KGBio))
      ),
      box(title = "Zuofu Huang", status = "danger", solidHeader = TRUE,
          width = 4,
          fluidRow(column(4, tags$img(height = 130, width = 100, src = "zuofu_pic.jpg"))),
          fluidRow(br(),column(width = 12, ZuofuBio))
      )
    )
  )
}

aboutUsTab <- function(input, output, session) {

}
