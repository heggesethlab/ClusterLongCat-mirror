# A module for the overview tab

overviewTabUI <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      box(title = "A quick overview", status = "warning", solidHeader = TRUE, withMathJax(includeMarkdown("tab_text/overview_sections/a_quick_overview.md")),width = 12
      )
    ),
    
    fluidRow(
      box(
        title = "What do we do?", status = "danger", includeMarkdown("tab_text/overview_sections/what_do_we_do.md"),
        solidHeader = TRUE, width = 12
      )
    ),
    fluidRow(
      tabBox(title = "Data Descriptions",
        tabPanel("Healthcare Utilization",includeMarkdown("tab_text/data_descriptions/health_data_desc.md")),
        tabPanel("Sleep",includeMarkdown("tab_text/data_descriptions/sleep_data_desc.md")),
        width = 12
      )
    ),
    
    fluidRow(
      tabBox(title = "Dissimilarity Metrics",
             tabPanel("Pros and Cons",includeHTML("tab_text/diss_metrics/compare.html")),
             # Using HTML instead of Markdown changes the font of the entire Shiny app!
             # tabPanel("Pros and Cons",includeMarkdown("tab_text/diss_metrics/compare.md")),
             tabPanel("Distribution Distance",includeMarkdown("tab_text/diss_metrics/prob_distribution.md")),
             tabPanel("Feature Distance",includeMarkdown("tab_text/diss_metrics/counts_attributes.md")),
             tabPanel("Edit Distance",includeMarkdown("tab_text/diss_metrics/edit_distance.md")),
             width = 12
      )
    ),
    fluidRow(
      tabBox(title = "Clustering Methods",
             tabPanel("Pros and Cons",includeHTML("tab_text/clustering_methods/compare.html")),
             # tabPanel("Pros and Cons",includeMarkdown("tab_text/clustering_methods/compare.md")),
             tabPanel("Distance-based Clustering", includeMarkdown("tab_text/clustering_methods/distance_based.md")),
             tabPanel("Model-based Clustering", includeMarkdown("tab_text/clustering_methods/model_based.md")),
             width = 12
      )
    ),
    fluidRow(
      tabBox(title = "Clustering Comparisons",
             tabPanel("Two-way Table",includeMarkdown("tab_text/clustering_comparisons/table.md")),
             tabPanel("Adjusted Rand Index",includeMarkdown("tab_text/clustering_comparisons/adjusted_rand.md")),
             tabPanel("Jaccard Index",includeMarkdown("tab_text/clustering_comparisons/jaccard.md")),
             tabPanel("Normalized Mutual Information",includeMarkdown("tab_text/clustering_comparisons/normalized_mutual_info.md")),
             tabPanel("Normalized Variation of Information",includeMarkdown("tab_text/clustering_comparisons/var_info.md")),
             width = 12
      )
    )
  )
}

overviewTab <- function(input, output, session){

}