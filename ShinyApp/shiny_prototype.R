# the shiny app
# need to condition on whether sequences are same length for several methods, eg hamming distance
# need to implement different number of states

# some libraries
library(bayesMCClust)
library(rmarkdown)
library(shiny)
library(ggplot2)
library(shinydashboard)
library(shinyMatrix)
library(dplyr)
library(tidyr)
library(igraph)
library(dendextend)
library(cluster)
library(markdown)
library(htmlTable)
library(shinybusy)
library(seriation)
library(abind)
library(TraMineR)
library(scales)
library(seqHMM)
library(stringr)
library(partitionComparison)
library(data.table)
library(viridis)
library(purrr)


# loading helper functions
source('tab_text/about_us/shiny_bios.R')
source("helper_functions/shiny_viz.R")
source("helper_functions/shiny_functions.R")

# loading modules
source("tab_modules/clustering_helper_modules.R")
source("tab_modules/overview_tab.R")
source("tab_modules/one_clustering_tab.R")
source("tab_modules/two_clustering_tab.R")
source("tab_modules/view_data_tab.R")
source("tab_modules/upload_data_tab.R")
source("tab_modules/about_us_tab.R")
source("tab_modules/references.R")


# prepping the data for use
tempSleepData <- read.csv("data/sleep_data_subset.csv")
tempHealthData <- read.csv("data/health_data_subset.csv")

supplied_data_sets_nr <- list(health_data = tempHealthData, sleep_data = tempSleepData)



ui <- dashboardPage(

  dashboardHeader(title = "Clustering Tools"),
  dashboardSidebar(
    tags$head(tags$style(HTML('.content-wrapper { height: 6000px !important;}'))), # height needs to depend on num plots!
    sidebarMenu(
      menuItem(tabName = "overview", text = "Overview"),
      menuItem(tabName = "one_clustering", text = "Single Clustering"),
      menuItem(tabName = "compare_clusters", text = "Compare Clusterings"),
      menuItem(tabName = "raw_data", text = "View Data"),
      menuItem(tabName = "upload_data", text = "Upload Data"),
      menuItem(tabName = "about", text = "About Us"),
      menuItem(tabName = "references", text = "References")
    )
  ),
  dashboardBody(
    add_busy_spinner(height = "30px", width = "30px"),
    tabItems(

      tabItem(
        tabName = "overview",
        overviewTabUI(id = "overview_tab")
      ),

      tabItem(
        tabName = "one_clustering",
        oneClusteringTabUI(id = "one_clustering_tab")
      ),

      tabItem(
        tabName = "compare_clusters",
        twoClusteringTabUI(id = "two_clustering_tab")
      ),

      tabItem(
        tabName = "raw_data",
        viewDataTabUI(id = "view_data_tab")
      ),

      tabItem(
        tabName = "upload_data",
        uploadDataTabUI(id = "upload_data_tab")
      ),

      tabItem(
        tabName = "about",
        aboutUsTabUI(id = "about_us_tab")
      ),
      tabItem(
        tabName = "references",
        referencesTabUI(id = "references")
      )
    )
  )
)


server <- shinyServer(function(input, output, session){

  supplied_data_sets <- do.call(reactiveValues, supplied_data_sets_nr)

  # only need to call modules that actually have a server component here
  callModule(oneClusteringTab, "one_clustering_tab", supplied_data_sets)
  callModule(twoClusteringTab, "two_clustering_tab", supplied_data_sets)

  callModule(viewDataTab, "view_data_tab", supplied_data_sets)
  new_data <- callModule(uploadDataTab, "upload_data_tab")
  observeEvent(new_data$data(),{
    supplied_data_sets[[new_data$name()]] <- new_data$data()
  })

  callModule(referencesTab, "references")
})

shinyApp(ui, server)
