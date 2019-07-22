# Module for two clustering tab


source("tab_modules/clustering_helper_modules.R")

twoClusteringTabUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      clusterControlInput(id = ns("two_clustering_controls_1"), box_width = 6),
      clusterControlInput(id = ns("two_clustering_controls_2"), box_width = 6)
    ),
    
    clusteringVisesUI(id = ns("two_clustering_vises"))
  )
}

twoClusteringTab <- function(input, output, session, supplied_data_sets) {
  two_clust_info_1 <- callModule(clusterControl, "two_clustering_controls_1", supplied_data_sets = supplied_data_sets)
  two_clust_info_2 <- callModule(clusterControl, "two_clustering_controls_2", supplied_data_sets = supplied_data_sets)
  
  callModule(clusteringVises, "two_clustering_vises", 
             clustering_data = list(two_clust_info_1[[1]], two_clust_info_2[[1]]), 
             two_per_row = TRUE, 
             run_now = list(two_clust_info_1[[2]],two_clust_info_2[[2]]),
             same_dataset = reactive({two_clust_info_1[[3]]() == two_clust_info_2[[3]]()}) )
}
