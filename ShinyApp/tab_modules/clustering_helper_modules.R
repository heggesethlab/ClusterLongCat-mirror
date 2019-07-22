
#' UI that lets users define the clustering
#'
#' @param id
#' @param box_width how wide should this box be
clusterControlInput <- function(id, box_width) {
  ns <- NS(id)

  box(
    width = box_width,
    title = "Clustering Controls",
    add_busy_spinner(),
    fluidRow(
      column(
        width = 6,
        uiOutput(ns("data_set_selector")),
        radioButtons(ns("format_model"), label = "Data Format for Modeling",
                     choices = list("Original Data (Warning: Run Time is Longer!)" = 'orig', "Compressed Data" = 'compress'),
                     selected = 'compress'),
        selectizeInput(
          inputId = ns("clust_method"),
          label = "Select Clustering Method",
          choices = list(
            `Distance-based` = c(`Partioning Around Medoids` = "PAM",
                                 `Hierarchical Clustering` = "HIER") ,
            `Model-based` = c(`Mixture of Markov Models (MLE)` = "MMM",
                              `Mixture of Markov Models (Bayesian Estimation)` = "MCC",
                              `Dirichlet Multinomial Clustering (Bayesian Estimation)` = "DMCC"
                              )
          )
        ),

        conditionalPanel(
          condition = "input.clust_method == 'HIER'",
          ns = ns,

          selectInput(
            ns("hier_controls"),
            label = "Select Linkage",
            choices = list(`Complete` = "complete", `Single` = "single",
                           `Average` = "average", `Ward` = "ward")
          )
        ),
        conditionalPanel(
          condition = "input.clust_method != 'MMM' & input.clust_method != 'DMCC' & input.clust_method != 'MCC'",
          ns = ns,
          selectInput(
            inputId = ns("dissim_metric"),
            label = "Select Dissimilarity Metric",
            choices = list(`Distribution Distance` = c(`Chi-Squared Distance` = "CHI2"),
                           `Feature Distance` = c(`Longest Common Subsequence` = "LCS"),
                           `Edit Distance` = c(`Gower Distance` = "Gower", `Generalized Hamming Distance*` = "HAM", `Optimal Matching` = "OM", `Optimal Matching Sensitive to Spell Length` = "OMslen")
                           )
          ),
          p("*Only works for datasets with equal length sequences. Will always use compressed data if sequences are not of equal length.", style = "font-size:12px")
        ),
        sliderInput(
          inputId = ns("clust_num"),
          label = "Input Number of Clusters",
          min = 1, max = 5, value = 2
        )
      ),

      column(
        width = 6,
        
        conditionalPanel(
          condition = "input.dissim_metric == 'CHI2'& input.clust_method != 'MMM' & input.clust_method != 'DMCC' & input.clust_method != 'MCC'",
          ns = ns,
          uiOutput(ns("step_selector")),
          radioButtons(ns("overlap_selector"), label = "Should the intervals overlap?",
                       choices = list("Yes" = 'true', "No" = 'false'),
                       selected = 'false'),
          p("*Bandwidth must be even when overlap = TRUE.", style = "font-size:12px")
        ),
        
        conditionalPanel(
          condition = "input.dissim_metric == 'OM' & input.clust_method != 'MMM' & input.clust_method != 'DMCC' & input.clust_method != 'MCC'",
          ns = ns,
          uiOutput(ns("indel_selector")),
          uiOutput(ns("indel_cost")),
          htmlTableWidgetOutput(ns("indel_panel"), height = "100px")
        ),
        
        conditionalPanel(
          condition =  "input.dissim_metric == 'OMslen' & input.clust_method != 'MMM' & input.clust_method != 'DMCC' & input.clust_method != 'MCC'",
          ns = ns,
          numericInput(
            ns("indel_selector_only_one"),
            "Input Indel Cost",
            min = 0, max = 100, value = 1
          )
        ),
        
        conditionalPanel(
          condition = "(input.dissim_metric == 'HAM' |  input.dissim_metric == 'OM' | input.dissim_metric == 'OMslen') & input.clust_method != 'MMM' & input.clust_method != 'DMCC' & input.clust_method != 'MCC'",
          ns = ns,

          uiOutput(ns("sm_selector")),
          uiOutput(ns("sm_cost")),
          htmlTableWidgetOutput(ns("sm_panel"), height='100%')
        ),

        conditionalPanel(
          condition = "input.dissim_metric == 'OMslen'",
          ns = ns, 
          HTML("<br/><br/>"),
          selectInput(
            ns("link_selector"),
            "Select Method to Compute Substitution Costs",
            choices = c(`Arithmetic Average` = "mean",  `Geometric Average` = "gmean")
          ),
          numericInput(
            ns("h_selector"),
            "Select Exponent Time Weight",
            min = 0, max = 1, value = .5
          )
        ),
        
        # MCMC
        conditionalPanel(
          condition =  "input.clust_method == 'MCC' | input.clust_method == 'DMCC'",
          ns = ns,
          h4(strong("MCMC Controls")),
          numericInput(
            ns("num_iters"),
            "Select Number of MCMC Iterations",
            min = 500, max = 50000, value = 2000
          ),
          numericInput(
            ns("num_burn_in"),
            "Select Number of Iterations for Burn In of MCMC",
            min = 0, max = 1000, value = 500
          ),
          h4(strong("Model Priors")),
          numericInput(
            ns("group_size_prior"),
            "Select Dirichlet-prior for group sizes",
            min = 0, max = 10, value = 4
          )
          
        ),
        # Priors
        conditionalPanel(
          condition = "input.clust_method == 'MCC'",
          ns = ns,
          numericInput(
            ns("mcc_diag_prior"),
            "Select Prior for Diagonal Elements of Transition Matrices",
            min = 0, max = 10, value = 1,
          ),
          numericInput(
            ns("mcc_off_diag_prior"),
            "Select Prior for Off Diagonal Elements of Transition Matrices",
            min = 0, max = 10, value = 1,
          )
        ),
        conditionalPanel(
          condition = "input.clust_method == 'DMCC'",
          ns = ns,
          numericInput(
            ns("dmcc_a0"),
            "Select Prior for Negative Multinomial beta Hyperparameter",
            min = 0, max = 10, value = 1,
          ),
          numericInput(
            ns("dmcc_alpha"),
            "Select Prior for Negative Multinomial alpha Hyperparameter",
            min = 0, max = 10, value = 1,
          ),
          numericInput(
            ns("dmcc_N0"),
            "Select Prior for Negative Multinomial N0 Hyperparameter",
            min = 0, max = 20, value = 5,
          )
        )
        
      )
    ),
    fluidRow(
      column(
        align = "center",
        width = 12,
        actionButton(ns("run"), "Run Clustering")
      )
    )
  )
}

#' Server for cluster controls that creates reactive data frame with clustered data
#'
#' @param input, output, session
clusterControl <- function(input, output, session, supplied_data_sets) {
  
  selected_data <- reactive({reactiveValuesToList(supplied_data_sets)[[(function() {
    if(is.null(input$selected_data)){
      return("health_data")
    } else {
      return(input$selected_data)
    }
  })()]]})
  
  default_controls <- reactiveValues(
    sm_matrix = c(1,1,1,1,1,1),
    indel_vec = c(1,1,1,1)
  )

  observeEvent({selected_data(); input$format_model}, {
    default_cost <- set_default_cost(selected_data(), format = input$format_model)
    default_controls$sm_matrix <- default_cost$sm_matrix
    default_controls$indel_vec <- default_cost$indel_vec
  })
  
  sub_mat_info <- reactiveValues(sub_mat = isolate({default_controls$sm_matrix}), # format it right.
                                 data_states = c("Home", "Hospice", "Hospital", "SNF"),
                                 pairs = apply(combn(c("Home", "Hospice", "Hospital", "SNF"), 2), 2, function(v) paste0(v, collapse=' <-> ')))
  
  indel_info <- reactiveValues(indel_vec = isolate({default_controls$indel_vec}),
                               data_states = c("Home", "Hospice", "Hospital", "SNF"))

  observeEvent({selected_data(); input$format_model},
               {
                 sub_mat_info$sub_mat <- default_controls$sm_matrix
                 sub_mat_info$data_states <- levels(as.factor(selected_data()$state))
                 sub_mat_info$pairs <- apply(combn(sub_mat_info$data_states, 2), 2, function(v) paste0(v, collapse=' <-> '))

                 indel_info$indel_vec <- default_controls$indel_vec
                 indel_info$data_states <- levels(as.factor(selected_data()$state))
               })

  observeEvent({input$selected_sub; input$input_sm_cost},
               {
                 sub_mat_info$sub_mat[which(sub_mat_info$pairs == input$selected_sub)] <- isolate({input$input_sm_cost})
               })
  
  observeEvent({input$selected_state; input$input_indel_cost},
               {
                 indel_info$indel_vec[which(indel_info$data_states == input$selected_state)] <- isolate({input$input_indel_cost})
               })

  output$data_set_selector <- renderUI({
    ns <- session$ns
    selectInput(
      ns("selected_data"),
      "Select Dataset",
      choices = names(reactiveValuesToList(supplied_data_sets)),
      selected = default_controls$selected_data)
  })
  
  output$step_selector <- renderUI({
    ns <- session$ns
    numericInput(
      inputId = ns("step_selector"),
      label = "Input Bandwidth",
      min = 0, max = ncol(sequence_data(data = selected_data(), format = input$format_model)), value = ncol(sequence_data(data = selected_data(), format = input$format_model))
    )
  })
  
  output$indel_selector <- renderUI({
    ns <- session$ns
    selectInput(
      ns("selected_state"),
      "Select Indel State",
      choices = sub_mat_info$data_states,
      selected = sub_mat_info$data_states[1])
  })
  
  output$indel_cost <- renderUI({
    ns <- session$ns
    numericInput(
      inputId = ns("input_indel_cost"),
      label = "Input Indel Cost",
      min = 0, max = 500, value = indel_info$indel_vec[[which(indel_info$data_states == input$selected_state)]]
    )
  })
  
  output$indel_panel <- renderHtmlTableWidget({
    indel_df <- t(data.frame(indel_info$indel_vec, row.names = indel_info$data_states))
    htmlTableWidget(indel_df, number_of_entries = 1, rnames = FALSE, height = "50px")
  })

  output$sm_selector <- renderUI({
    ns <- session$ns
    selectInput(
      ns("selected_sub"),
      "Select Substitution",
      choices = sub_mat_info$pairs,
      selected = sub_mat_info$pairs[1])
  })

  output$sm_cost <- renderUI({
    ns <- session$ns
    numericInput(
      inputId = ns("input_sm_cost"),
      label = "Input Substitution Cost",
      min = 0, max = 100, value = sub_mat_info$sub_mat[[which(sub_mat_info$pairs == input$selected_sub)]]
    )
  })

  output$sm_panel <- renderHtmlTableWidget({
    mat <- matrix(0, nrow = length(sub_mat_info$data_states), ncol = length(sub_mat_info$data_states))
    mat[lower.tri(mat, diag=FALSE)] <- sub_mat_info$sub_mat
    colnames(mat) <- rownames(mat) <- sub_mat_info$data_states
    M <- makeSymm(mat)
    htmlTableWidget(as.data.frame(M), number_of_entries = 4)
  })

  clustering_info <- eventReactive(input$run, {
    if(input$run > 0){
      data <- isolate(reactiveValuesToList(supplied_data_sets)[[(function() {
        if(is.null(isolate(input$selected_data))){
          return(default_controls$selected_data)
        } else {
          return(isolate(input$selected_data))
        }
      })()]])}
      
      mat <- matrix(0, nrow = length(sub_mat_info$data_states), ncol = length(sub_mat_info$data_states))
      mat[lower.tri(mat, diag=FALSE)] <- sub_mat_info$sub_mat
      colnames(mat) <- rownames(mat) <- sub_mat_info$data_states
      sub_cost_matrix <- makeSymm(mat)
      
      indel_cost_vector <- indel_info$indel_vec
      
      data %>%
        cluster_data(cluster_method = isolate(input$clust_method),
                     cluster_params = list(hier_method = isolate(input$hier_controls),
                                           num_iters = input$num_iters,
                                           num_burn_in = input$num_burn_in,
                                           group_size_prior = input$group_size_prior,
                                           mcc_diag_prior = input$mcc_diag_prior,
                                           mcc_off_diag_prior = input$mcc_off_diag_prior,
                                           dmcc_a0 = input$dmcc_a0,
                                           dmcc_alpha = input$dmcc_alpha,
                                           dmcc_N0 = input$dmcc_N0),
                     num_clusts = isolate(input$clust_num),
                     distance_method = isolate(input$dissim_metric),
                     distance_params = list(
                       step = input$step_selector,
                       overlap = (function() {
                         if(input$overlap_selector == 'true'){
                           return(TRUE)
                         } else {
                           return(FALSE)
                         }
                       })(),
                       sm_matrix = (function() {
                         if(is.null(sub_cost_matrix)){
                           return(default_controls$sm_matrix)
                         } else {
                           return(sub_cost_matrix)
                         }
                       })(),
                       indel = (function() {
                         if(input$dissim_metric == "OMslen"){
                           return(input$indel_selector_only_one)
                         }
                         if(is.null(indel_cost_vector)) {
                           return(default_controls$indel_vec)
                         } else {
                           return(indel_cost_vector)
                         }
                       })(),
                       link = input$link_selector,
                       h = input$h_selector
                     ),
                     data_format_model = (function(data){
                       if(isolate(input$format_model) == "compress") {
                         return("compress") # if format is compressed we're fine
                       } else if ( is.null(input$dissim_metric) | !(input$dissim_metric %in% c("HAM")) ) {
                         return(isolate({input$format_model})) # if dissim metric isn't HAM we're fine
                       } else if(
                         (data %>% 
                           group_by(ID) %>% 
                           summarise(seq_length = sum(duration)))$seq_length %>% 
                           (function(x) ((x[1] == x))) %>% 
                           all()
                       ) {
                         return("orig")
                       } else {
                         return("compress")
                       }
                     })(.)
        )
      
    })


  run_now <- reactive({input$run})


  return(list(clustering_info, run_now, reactive({input$selected_data})))
}

#' UI for module that makes single visualizations
#'
#' @param id
#' @param width int the width of the visualization
clusteringVisUI <- function(id, width) {
  ns <- NS(id)

  column(
    width = width,
    plotOutput(
      ns("clusteringVis")
    )
  )
}

#' Server for module that makes a single visualization
#'
#' @param input, output, session
#' @param clustering_info list (reactive) of the data required to plot this visualization
#' @param vis_type string the type of visualization to display
clusteringVis <- function(input, output, session, clustering_info, vis_type, data_format_vis) {

  output$clusteringVis <- renderPlot({
    model_type <- clustering_info()$model_type
    switch (vis_type,
            "num_ids" = plot_num_ids(clustering_info()),
            "state_dist" = plot_state_dist(clustering_info()),
            "total_dur" = plot_total_duration(clustering_info()),
            "num_trans" = plot_avg_num_spells(clustering_info()),
            "samp_traj" = plot_sample_traj(clustering_info(),data_format_vis),
            "avg_dur" = plot_avg_state_duration(clustering_info(),data_format_vis),
            "state_dist_time" = plot_prop_time(clustering_info(),data_format_vis),
            "mark_chain" = if(model_type %in% c("MMM","DMCC", "MCC")){
              plot_mark_chain(clustering_info())
            },
            "mark_heat_map" = if(model_type %in% c("MMM","DMCC", "MCC")){
              if(model_type == "MMM"){
                plot_mark_heat_map(clustering_info())
              } else {
                plot_mcc_trans_heat_map(clustering_info())
              }
              
            },
            "mark_heat_map_sd" = if(model_type %in% c("DMCC", "MCC")){
              plot_mcc_sd_heat_map(clustering_info())
            },
            "mark_post_probs" = if(model_type %in% c("MMM","DMCC", "MCC")){
              plot_mark_post_probs(clustering_info())
            },
            "dend" = if(model_type == "HIER"){
              plot_dendro(clustering_info())
            }, 
            "dissim_mat" = if(!(model_type %in% c("MMM", "MCC", "DMCC"))) {
              plot_dissim_matrix(clustering_info())
            },
            "var_trans_prob" = if(model_type %in% c("DMCC")){
              plot_variance(clustering_info())
            },
            "unob_het" = if(model_type %in% c("DMCC")){
              plot_unob_het(clustering_info())
            }, 
            "group_trace" = if(model_type %in% c("DMCC", "MCC")){
              plot_group_size_trace(clustering_info())
            }, 
            "diag_trace" = if(model_type %in% c("DMCC", "MCC")) {
              plot_diag_trace(clustering_info())
            },
            "plot_group_het_trace" = if(model_type %in% c("DMCC")) {
              plot_group_het_trace(clustering_info())
            }

    )

    }, height = ifelse(vis_type == "mark_chain", ceiling(clustering_info()$num_clusts/2)*400, "auto")
  )
}

#' UI for module that displays visualizations
#'
#' @param id, \code{shiny} id
#' @return taglist with a selector for visualizations to display and a set of visualizations
clusteringVisesUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("download_data_button")),
    uiOutput(ns("vis_extras")),
    fluidRow(
      uiOutput(ns("vises"))
    )
  )
}

#' Server for module that displays visualizations
#'
#' @param input, output, session standard \code{shiny} boilerplate
#' @param clustering_data list that contains reactive data frames which hold the data required to plot
#' @param two_per_row boolean true if we plot two visualizations per row
clusteringVises <- function(input, output, session, clustering_data = NULL, two_per_row, run_now, same_dataset = NULL) {

  output$vis_extras <- renderUI({
    if((length(run_now) == 1 & sum(sapply(run_now, function(v) v()))  > 0) | (length(run_now) > 1 & sum(sapply(run_now, function(v) v())>0)  > 1)){
      ns <- session$ns
      do.call(tagList, list(
        fluidRow(
          (function(){if(two_per_row){
            if(isolate({same_dataset()})) {
              box(
                width = 12,
                title = "Clustering Comparison",
                height = "450px",
                tags$b('Two-Way Table of Cluster Labels'),
                tags$br(),tags$br(),
                htmlTableWidgetOutput(ns("clustering_table"), height='100%'),
                tags$br(),tags$br(),
                tags$b('Comparison Indices'),
                tags$br(),tags$br(),
                htmlTableWidgetOutput(ns("clustering_comparison"))
              )
            } else {
              box(
                width = 12,
                title = "Clustering Comparison",
                strong("The same datasets must be used in order to calculate comparison statistics.", style = "color:red")
              )
            }
          }
          })(),
          box(
            title = "Clustering Visualizations",
            width = 12,
            uiOutput(ns("vis_types"))
          )
          )

        )

      )

    }
  })

  output$vis_types <- renderUI({
    if((length(run_now) == 1 & sum(sapply(run_now, function(v) v()))  > 0) | (length(run_now) > 1 & sum(sapply(run_now, function(v) v()))  > 1)){
    ns <- session$ns
    vis_choices <- c(`Cluster Distribution` = "num_ids",
                     `State Distribution` = "state_dist",
                     `State Distribution over Time` = "state_dist_time",
                     `Total Duration Distribution` = "total_dur",
                     `Average Duration by State` = "avg_dur",
                     `Average Count of Spells by State` = "num_trans",
                     `Sample Trajectories` = "samp_traj")
    default = "num_ids"

    temp_list <- clustering_data[[1]]
    temp_val <- temp_list()$model_type

    if(length(clustering_data) > 1) {
      temp_list <- clustering_data[[2]]
      temp_val <- union(temp_val, temp_list()$model_type)
    }

    if(!is.null(temp_val)){
      if("MMM" %in% temp_val){
        vis_choices <- c(vis_choices, 
                         `Posterior Cluster Probabilities` = "mark_post_probs",
                         `Markov Chains` = "mark_chain",
                         `Transition Probabilities` = "mark_heat_map")
                         
        default <- union(default, c("mark_post_probs", "mark_chain", "mark_heat_map"))
      }
      if("DMCC" %in% temp_val) {
        vis_choices <- c(vis_choices, 
                         `Posterior Cluster Probabilities` = "mark_post_probs",
                         `Markov Chains` = "mark_chain",
                         `Transition Probabilities` = "mark_heat_map",
                         `SD of Transition Probabilities` = "mark_heat_map_sd",
                         `Variance of Transition Probabilities` = "var_trans_prob",
                         `Unobservered Heterogeneity` = "unob_het",
                         `Group Size Trace Plot` = "group_trace",
                         `Diagonal Transitions Trace Plot` = "diag_trace",
                         `Unobsereved Within Group Heterogeneity Trace Plot` = "plot_group_het_trace")
        
        default <- union(default, c("mark_post_probs", "mark_chain", "mark_heat_map"))
      }
      if("MCC" %in% temp_val) {
        vis_choices <- c(vis_choices, 
                         `Posterior Cluster Probabilities` = "mark_post_probs",
                         `Markov Chains` = "mark_chain",
                         `Transition Probabilities` = "mark_heat_map",
                         `SD of Transition Probabilities` = "mark_heat_map_sd",
                         `Group Size Trace Plot` = "group_trace",
                         `Diagonal Transitions Trace Plot` = "diag_trace")
        
        default <- union(default, c("mark_post_probs", "mark_chain", "mark_heat_map"))
      } 
      if(sum( c("HIER", "PAM")%in% temp_val) > 0 ){ #if one is dissim based
        vis_choices <- c(vis_choices, `Dissimilarity Matrix` = "dissim_mat")
        default <- union(default, "dissim_mat")
      }
      if("HIER" %in% temp_val) {
        vis_choices <- c(vis_choices, `Dendrogram` = "dend")
        default <- union(default, "dend")
      }

    }

    tagList(
      div(selectInput(
        ns("vis_types"),
        label = "Select Visualizations",
        choices = vis_choices,
        multiple = TRUE,
        selected = default
      ),
      radioButtons(ns("format_viz"), label = "Data Format for Visualization",
        choices = list("Original Data" = 'orig', "Compressed Data" = 'compress'),
        selected = 'compress')
      ),
      column(
        width = 12,
        align = "center",
        actionButton(ns("run_visualizations"),
                     "Visualize")
      ))
    }
  })

  if(two_per_row){
    observeEvent({run_now[[1]](); run_now[[2]]()},{
      if((length(run_now) == 1 & sum(sapply(run_now, function(v) v()))  > 0) | (length(run_now) > 1 & sum(sapply(run_now, function(v) v()))  > 1)){
        output$clustering_comparison <- renderHtmlTableWidget({htmlTableWidget(compare_clustering(isolate(clustering_data[[1]]()), isolate(clustering_data[[2]]())),
                                                                             number_of_entries = 1,
                                                                             rnames = FALSE)})
        #Add matrix of cluster labels as output: clust_table(isolate(clustering_data[[1]]()), isolate(clustering_data[[2]]()))
        output$clustering_table <- renderHtmlTableWidget({htmlTableWidget(clust_table(isolate(clustering_data[[1]]()), isolate(clustering_data[[2]]())), number_of_entries = 6)})
      }
    })
  }
  
  if(!two_per_row){
    observeEvent({run_now[[1]]()},{
      if((length(run_now) == 1 & sum(sapply(run_now, function(v) v()))  > 0) | (length(run_now) > 1 & sum(sapply(run_now, function(v) v()))  > 1)){
        output$download_data_button <- renderUI({
          ns <- session$ns
          
          fluidRow(
            # column(
            #   width = 4,
            #   offset = 4,
              box(
                width = 12,
                title = "Download Data for Further Analysis",
                column(
                  width = 12,
                  align= "center",
                  "Download the original supplied data set with cluster membership for this clustering adjoined to it.",
                  HTML("<br/><br/>"),
                  downloadButton(
                    ns("download_data"),
                    label = "Download"
                  )
                )
              )
            # )
            
          )
        })
        
        
        output$download_data <- downloadHandler(
          filename = "clustered_data.csv",
          content = function(file) {
            write.csv(clustering_data[[1]]()$clustered_data, file, row.names = FALSE)
          }
        )
      }
    })
  }
  
  output$vises <- renderUI({
    if((length(run_now) == 1 & sum(sapply(run_now, function(v) v()))  > 0) | (length(run_now) > 1 & sum(sapply(run_now, function(v) v()))  > 1)){
    ns <- session$ns


    num_selections <- isolate({length(input$vis_types)})




    ui_list <- eventReactive(input$run_visualizations,{

      lapply(min(1, num_selections):num_selections, function(i){
        if(num_selections > 0){
          if(!two_per_row) {
            box(
              width = 12,
              column(
                width = 12,
                callModule(clusteringVis,
                           ns("one_cluster_vis"),
                           clustering_info = clustering_data[[1]],
                           vis_type = input$vis_types[i],
                           data_format_vis = isolate(input$format_viz))
              )
            )
          } else {
            box(
              width = 12,
              column(
                width = 6,
                callModule(clusteringVis,
                           ns("two_cluster_vis_1"),
                           clustering_info = clustering_data[[1]],
                           vis_type = input$vis_types[i],
                           data_format_vis = isolate(input$format_viz))
              ),
              column(
                width = 6,
                callModule(clusteringVis,
                           ns("two_cluster_vis_2"),
                           clustering_info = clustering_data[[2]],
                           vis_type = input$vis_types[i],
                           data_format_vis = isolate(input$format_viz))
              )
            )
          }
        }
      })
    })

    do.call(tagList, ui_list())
    }
  })

}
