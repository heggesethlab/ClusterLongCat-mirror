# Module for one cluster stuff
source("tab_modules/clustering_helper_modules.R")
oneClusteringTabUI <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      clusterControlInput(id = ns("one_clustering_controls"), box_width = 12)
    ),
    
    clusteringVisesUI(id = ns("one_clustering_vises"))
  )
}

oneClusteringTab <- function(input, output, session, supplied_data_sets){
  one_clust_info <- callModule(clusterControl, "one_clustering_controls", supplied_data_sets = supplied_data_sets)
  callModule(clusteringVises, "one_clustering_vises", clustering_data = list(one_clust_info[[1]]), two_per_row = FALSE, run_now = list(one_clust_info[[2]]))
}