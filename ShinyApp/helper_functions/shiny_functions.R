library(tidyr)
library(dplyr)
library(cluster)
library(TraMineR)
library(scales)
library(seqHMM)
library(stringr)
library(partitionComparison)
library(data.table)

###### REFERENCES ######
#A function to make a symmetric function out of a lower triangular matrix
makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

# A function to make each row a (Day, ID) pair
data_to_daily <- function(data){
  # grabs a random sample (without replacement) and makes it so each row is a (day, id) pair
  duration_joiner <- tibble(
    order = rep(data$order, data$Duration),
    ID = rep(data$ID, data$Duration)
  )

  data_long <- data %>%
    full_join(duration.joiner, by = c("ID", "order")) %>%
    group_by(ID) %>%
    mutate(time = 1:n()) %>%
    ungroup()

  return(data_long)
}
#######################

#' Formats the healthcare utilization data 
#'
#' @param dataframe, a dataset containing ID, order, setting, Duration, FromDate, ThruDate, DurationRel, FromDateRel, ThruDateRel variables
#' @return a dataframe containing ID, order, state, start_time, duration, start_time_comp, duration_comp variables
format_health_data <- function(data) {
  data <- data %>%
    rename(state = setting, duration = Duration, start_time = FromDate, start_time_rel = FromDateRel, end_time_rel = ThruDateRel) %>%
    mutate(state = recode(state, "NoCoveredCare" = "Home")) %>%
    group_by(ID) %>%
    mutate(total_duration = sum(duration)) %>%
    filter(total_duration < 1000) %>%
    ungroup %>%
    sample_data(1000) %>%
    data_to_time_prop(100) %>%
    group_by(ID, order) %>%
    mutate(duration_comp = round((max(time)-min(time))*100)+1, start_time_comp = round(min(time)*100))

  
  data <- data %>%
    distinct(ID, order, state, start_time, duration, start_time_comp, duration_comp)
  
  return(data)
}

#' Formats the Fitbit sleep data 
#'
#' @param dataframe, a dataset containing DateTime, Stage, Time, Day variables
#' @return a dataframe containing ID, order, state, start_time, duration, start_time_comp, duration_comp variables
format_sleep_data <- function(data) {
  data <- data %>%
    rename(state = Stage, duration = Time) %>%
    mutate(state = recode(state, "rem" = "REM"))
  
  data$ID <- cumsum(c(0, diff(data$DateTime)) < 0) + 1
  data$duration <- data$duration / 60
  data$start_time <- unlist(tapply(data$duration, data$ID, function(v) c(0, cumsum(v)[-length(v)])))
  
  data <- data %>%
    group_by(ID) %>%
    mutate(duration_rel = duration/sum(duration))
  
  data$start_time_rel <- unlist(tapply(data$duration_rel, data$ID, function(v) c(0, cumsum(v)[-length(v)])))
  data$end_time_rel <- unlist(tapply(data$duration_rel, data$ID, function(v) c(cumsum(v))))
    
  data <- data %>%
    group_by(ID) %>%
    mutate(order = 1:n()) %>%
    ungroup()
  
  data <- data %>%
    data_to_time_prop(100) %>%
    group_by(ID, order) %>%
    mutate(duration_comp = round((max(time)-min(time))*100)+1, 
           start_time_comp = round(min(time)*100)) %>%
    ungroup()
  
  data <- data %>%
    distinct(ID, order, state, start_time, duration, start_time_comp, duration_comp)
  
  return(data)
}

#' Takes a random sample of the data
#'
#' @param dataframe, a dataset containing an ID variable
#' @param size, an integer indicating the sample size 
#' @return a dataframe of smaller size containing ID and the remaining original variables
sample_data <- function(data, size){
  data_sample <- sample(unique(data$ID), size = min(size, length(unique(data$ID))))
  
  data_small <- data %>%
    filter(ID  %in% data_sample)
  
  return(data_small)
}

#' Repeats things as many times as it appears
#'
#' @param dataframe, a dataset containing ID, order, state, duration, start_time, duration_rel, start_time_rel, end_time_rel variables
#' @param timechunks, an integer indicating the number of timechunks
#' @return a dataframe contatining ID, order, state, duration, start_time, duration_rel, start_time_rel, end_time_rel, time variables  
data_to_time_prop <- function(data, timechunks){
  foo <- trunc(data$end_time_rel*timechunks+10^(-13)) - trunc(data$start_time_rel*timechunks+10^(-13))  #Fix this again
  grid_joiner <- tibble(
    order = rep(data$order, foo),
    ID = rep(data$ID, foo )
  )
  
  grid_joiner <- grid_joiner %>%
    group_by(ID) %>%
    mutate(time = 0:(timechunks-1) / timechunks)
  
  
  data_long <- data %>%
    full_join(grid_joiner, by = c("ID", "order")) %>%
    filter(!is.na(time)) # gets rid of visits that lasted no time
  
  
  return(data_long)
}

#' Turns data to long format 
#'
#' @param dataframe, a dataset containing ID, state, duration, start_time variables
#' @return a dataframe containing ID, state, time variables
data_to_long <- function(data, format = NULL){
  if(!is.null(format))
    if(format == 'compress'){
    data <- data %>%
      dplyr::select(-start_time,-duration) %>%
      rename(start_time = start_time_comp,duration = duration_comp)
  }
  
  data <- data %>%
    group_by(ID) %>% 
    mutate(order = 1:n()) %>% 
    ungroup()
  
  time_expander <- tibble(ID = rep(data$ID, data$duration), start_time = rep(data$start_time, data$duration))  
  data %>% 
    #select(ID, state, duration, start_time) %>% 
    full_join(time_expander, by = c("ID", "start_time")) %>% 
    group_by(ID) %>% 
    mutate(time = 1:n()) %>% 
    ungroup()
}

#' Turns data to wide format  
#'
#' @param dataframe, a dataset containing ID, state, time variables
#' @return a dataframe containing ID and time_n variables with state values
data_to_wide <- function(data){
  data_wide <- data %>%
    dplyr::select(ID, state, time) %>%
    spread(key = time, value = state)
  
  return(data_wide)
}

#' Turns data into a state sequence object
#'
#' @param dataframe, a dataset containing ID, state, duration, start_time variables
#' @return a dataframe containing an ID and its corresponding sequence
sequence_data <- function(data, format = NULL) {
  sequences <- data %>%
    data_to_long(format) %>% 
    data_to_wide() %>%
    seqdef(var = 2:ncol(.), id = .$ID)

  return(sequences)
}

# A function to create transition matrix for MMM
create_trans_matrix<- function(n){
  # creates an nxn matrix with a higher weight towards staying in the same state
  trans_matrix <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    trans_matrix[i,] = runif(n, .05, 1/(2*n))
    trans_matrix[i,i] = 1-sum(trans_matrix[i,])+trans_matrix[i,i]
  }

  return(trans_matrix)
}

# A function to create initial probabilities for MMM
create_initial_probs <- function(n){
  initial_probs <- runif(n, 0.05, 1/n)
  initial_probs[1] <- 1-sum(initial_probs)+initial_probs[1]
  return(initial_probs)
}

# A function to set default indel and substitution costs
set_default_cost <- function(data, format = NULL){
  
  if(!is.null(format))
    if(format == 'compress'){
      data <- data %>%
        dplyr::select(-start_time,-duration) %>%
        rename(start_time = start_time_comp,duration = duration_comp)
    }
  
  # indel
  state_frequency <- data %>%
    group_by(state) %>%
    summarise(state_duration = sum(duration))
  indel_cost <- as.vector(((state_frequency$state_duration))^(-0.5))
  scale_times_indel <- 1/(min(indel_cost))
  indel_vec <- round(scale_times_indel * indel_cost, 3)
  
  # substitution
  sequences <- sequence_data(data)
  sub_cost <- seqtrate(sequences)
  sub_cost_trans <- t(sub_cost)
  sub_upper <- sub_cost_trans[!upper.tri(sub_cost_trans,diag=TRUE)]
  sub_lower <- sub_cost[!upper.tri(sub_cost,diag=TRUE)]
  sub_long <- 2 * (1 - pmax(sub_upper, sub_lower))
  
  scale_times_substitution <- 2/(max(sub_long[sub_long > 0]))
  sm_matrix <- round(sub_long * scale_times_substitution, 3)
  
  return(list(indel_vec = indel_vec, sm_matrix = sm_matrix))
}

# A function to compute the selected dissimilarity/distance matrix
create_dist_matrix <- function(sequences, data = NULL, distance_method, distance_params) {

  switch(distance_method,
         "CHI2" = {
                  dist_matrix <- seqdist(sequences, method = "CHI2", breaks = NULL, step = distance_params$step, overlap = distance_params$overlap)
                  },
         "HAM" = {
                  dist_matrix <- seqdist(sequences, method = "HAM", sm = distance_params$sm_matrix)
                 },
         "OM" = {
                  dist_matrix <- seqdist(sequences, method = "OM", sm = distance_params$sm_matrix, indel = distance_params$indel)
                },
         "OMslen" = {
                  dist_matrix <- seqdist(sequences, 
                                         method = "OMslen", 
                                         sm = distance_params$sm_matrix,
                                         indel = as.numeric(distance_params$indel), # some set number,
                                         h = distance_params$h,
                                         link = distance_params$link)
                    },
         "LCS" = {
                   dist_matrix <- seqdist(sequences, method = "LCS")
                 },
         "Gower" = {
                   dist_matrix <- data %>% 
                            data_to_long(format = distance_params$data_format_model) %>%
                            data_to_wide() %>%
                            dplyr::select(-ID) %>% 
                            apply(.,2,function(v) as.factor(v)) %>% as.data.frame() %>%
                            daisy(., metric = "gower") %>% as.matrix()
         }
  )

  dist_matrix <- dist_matrix / max(dist_matrix) * 100
  dist_matrix
}


# A function to perform the selected clustering method
perform_clustering <- function(cluster_method, cluster_params) {

  switch(cluster_method,
         "PAM"  = {
                    clust_result <- pam(cluster_params$dist_matrix, diss = TRUE, k = cluster_params$num_clusts)
                  },
         "HIER" = {
                    clust_result <- agnes(cluster_params$dist_matrix, diss = TRUE, method = cluster_params$hier_method)
                  },
         "MMM"  = {
                    trans_probs <- lapply(rep(length(alphabet(cluster_params$sequences)), cluster_params$num_clusts), create_trans_matrix)
                    initial_probs <- lapply(rep(length(alphabet(cluster_params$sequences)), cluster_params$num_clusts), create_initial_probs)
                    
                    num_cores <- parallel::detectCores()
                    if(is.na(num_cores)){
                      num_cores <- 1
                    }

                    clust_result <- build_mmm(cluster_params$sequences, transition_probs = trans_probs, initial_probs = initial_probs) %>%
                      fit_model(threads = num_cores, log_space = FALSE, control_em = list(restart = list(times = 10, emission = FALSE)))
                  },
         "MCC"  = {
                    data <- cluster_params$data
                    formatted_data <- data %>% 
                    data_to_long(format = cluster_params$data_format_model) %>% 
                    dplyr::select(ID, state, time)
                     
                    formatted_array <- formatted_data %>%
                      group_by(ID) %>% 
                      transmute(from = state, to = lead(state)) %>% 
                      filter(!is.na(to)) %>% 
                      ungroup() %>%
                      count(from, to, ID) %>% 
                      spread(from, n, fill = 0) %>% 
                      gather(from, n, 3:6) %>% 
                      spread(to, n, fill = 0) %>% 
                      select(ID, from, sort(unique(.$from))) %>% 
                      arrange(ID, from) %>% 
                      split(.$ID) %>% 
                      lapply(function(curr_matrix){
                        state_names <- curr_matrix[[2]]
                        curr_matrix <- curr_matrix[,3:6]
                        curr_matrix <- as.matrix(curr_matrix)
                        rownames(curr_matrix) <- state_names
                        curr_matrix
                      }) %>% 
                      abind(along =  3)
                     
                      num_states <- length(unique(formatted_data$state))
                     
                      dimnames(formatted_array) <- list(1:num_states, 1:num_states, NULL)
                      
                     # running the model
                      tryCatch(clust_result <- mcClust(
                        Data = list(dataFile = formatted_array,
                                    storeDir = "temp_outfiles"
                        ),
                        Prior = list( H = cluster_params$num_clusts, 
                                      e0 = cluster_params$group_size_prior,
                                      c = cluster_params$mcc_diag_prior,
                                      cOff = cluster_params$mcc_off_diag_prior,
                                      usePriorFile = FALSE
                        ),
                        Initial = list( xi.start.ind = 2),
                        Mcmc = list( M = cluster_params$num_iters,
                                     M0 = cluster_params$num_burn_in,
                                     mOut = 100,
                                     mSave = 50,
                                     monitor = FALSE, 
                                     seed = sample(1:10000, 1))
                      ),
                               error = function(c) stop("The clustering method failed")
                      )
                      
                  },
         "DMCC" = {
                    # data prep
                    data <- cluster_params$data
                    formatted_data <- data %>% 
                      data_to_long(format = cluster_params$data_format_model) %>% 
                      dplyr::select(ID, state, time)

                    formatted_array <- formatted_data %>%
                      group_by(ID) %>% 
                      transmute(from = state, to = lead(state)) %>% 
                      filter(!is.na(to)) %>% 
                      ungroup() %>%
                      count(from, to, ID) %>% 
                      spread(from, n, fill = 0) %>% 
                      gather(from, n, 3:6) %>% 
                      spread(to, n, fill = 0) %>% 
                      select(ID, from, sort(unique(.$from))) %>% 
                      arrange(ID, from) %>% 
                      split(.$ID) %>% 
                      lapply(function(curr_matrix){
                        state_names <- curr_matrix[[2]]
                        curr_matrix <- curr_matrix[,3:6]
                        curr_matrix <- as.matrix(curr_matrix)
                        rownames(curr_matrix) <- state_names
                        curr_matrix
                      }) %>% 
                      abind(along =  3)
                    
                    num_states <- length(unique(formatted_data$state))
                    
                    dimnames(formatted_array) <- list(1:num_states, 1:num_states, NULL)
                    
                    # running the model
                    tryCatch(clust_result <- dmClust(
                      Data = list(dataFile = formatted_array,
                                  storeDir = "temp_outfiles", 
                                  mccFile = NULL),
                      Prior = list( H = cluster_params$num_clusts, 
                                    alpha0 = cluster_params$group_size_prior,
                                    a0 = cluster_params$dmcc_a0,
                                    alpha = cluster_params$dmcc_alpha,
                                    N0 = cluster_params$dmcc_N0,
                                    isPriorNegBin = FALSE,
                                    mccAsPrior = FALSE,
                                    persPrior = 0.7),
                      Initial = list( mccUse = FALSE,
                                      pers = 0.7),
                      Mcmc = list( kNo = 2,
                                   M = cluster_params$num_iters,
                                   M0 = cluster_params$num_burn_in,
                                   mOut = 100,
                                   mSave = 50,
                                   showAcc = TRUE,
                                   monitor = FALSE, 
                                   seed = sample(1:10000, 1))
                    ),
                      error = function(c) stop("The clustering method failed")
                    )
                  }
  )


  return(clust_result)
}

#  distance_params can include sm_matrix or indel
cluster_data <- function(data, cluster_method, num_clusts, cluster_params = NULL, distance_method = NULL, distance_params = NULL, data_format_model = NULL) {
  clustering_info <- list()
  
  if( !(cluster_method %in% c("MCC", "DMCC")) ){
    sequences <- sequence_data(data, format = data_format_model)
  } else {
    cluster_params$data <- data
    cluster_params$data_format_model <- data_format_model
    sequences <- NULL
  }
  
  if(is.null(cluster_params)){
    cluster_params <- list()
  }
  
  if(distance_method == "Gower"){
    distance_params$data_format_model <- data_format_model
  }
  
  if( !(cluster_method %in% c("MMM", "MCC", "DMCC")) ){
    cluster_params$dist_matrix <- create_dist_matrix(sequences, data = data, distance_method = distance_method, distance_params = distance_params)
    clustering_info$dist_matrix <- cluster_params$dist_matrix
  } else {
    cluster_params$sequences <- sequences
  }

  cluster_params$num_clusts <- num_clusts

  clustering <- perform_clustering(cluster_method = cluster_method, cluster_params)


  switch(cluster_method,
         "PAM"  = {
                    clustering_info$clustered_data <- tibble(ID = unique(data$ID), cluster = clustering$clustering) %>%
                      right_join(data, by = "ID") %>%
                      dplyr::select(-cluster, cluster)
                  },
         "HIER" = {
                    hier_cut <- cutree(clustering, k = num_clusts)
                    clustering_info$clustered_data <- tibble(ID = unique(data$ID), cluster = hier_cut) %>%
                      right_join(data, by = "ID") %>%
                      dplyr::select(-cluster, cluster)
                    clustering_info$cluster_output <- clustering
                  },
         "MMM"  = {
                    mmm_summary <- summary(clustering$model)
                    clustering_info$clustered_data <- tibble(ID = unique(data$ID), cluster = mmm_summary$most_probable_cluster) %>%
                      right_join(data, by = "ID") %>%
                      dplyr::select(-cluster, cluster) %>%
                      mutate(cluster = str_sub(cluster, 9))

                    post_clust_probs <- mmm_summary$posterior_cluster_probabilities
                    post_clust_probs <- as.data.frame(post_clust_probs)
                    names(post_clust_probs) <- paste0("post_prob_cluster_", 1:num_clusts)
                    post_clust_probs$ID <- unique(data$ID)

                    clustering_info$clustered_data <- clustering_info$clustered_data %>%
                      left_join(post_clust_probs, by = "ID")

                    clustering_info$model <- clustering$model
                  },
         "MCC"  = {
                    allocList <- calcAllocationsMCC(clustering, thin=1, maxi=round((cluster_params$num_iters-cluster_params$num_burn_in) / 10), plotPathsForEta = FALSE)
                     
                    clustering_info$clustered_data <- tibble(ID = unique(data$ID), cluster = allocList$class) %>% 
                      right_join(data, by = "ID") 
                     
                    post_clust_probs <- as.data.frame(allocList$classProbs)
                  
                    names(post_clust_probs) <- paste0("post_prob_cluster_", 1:num_clusts)
                    post_clust_probs$ID <- unique(data$ID)
                     
                    clustering_info$clustered_data <- clustering_info$clustered_data %>% 
                      left_join(post_clust_probs, by = "ID")
                     
                    clustering_info$transition_data <- calcTransProbs(clustering, estGroupSize=allocList$estGroupSize, thin=1, 
                                                                      printXtable=FALSE, printSd=FALSE, printTogether=FALSE,
                                                                      plotPaths = FALSE, 
                                                                      plotPathsForE = FALSE)
                    
                    clustering_info$model <- clustering
                  },
         "DMCC" = {
                    allocList <- calcAllocationsDMC(clustering, thin=1, maxi=round((cluster_params$num_iters-cluster_params$num_burn_in) / 10), plotPathsForEta = FALSE)
           
                    clustering_info$clustered_data <- tibble(ID = unique(data$ID), cluster = allocList$class) %>% 
                      right_join(data, by = "ID") 

                    post_clust_probs <- as.data.frame(allocList$classProbs)
                    names(post_clust_probs) <- paste0("post_prob_cluster_", 1:num_clusts)
                    post_clust_probs$ID <- unique(data$ID)
                    
                    clustering_info$clustered_data <- clustering_info$clustered_data %>% 
                      left_join(post_clust_probs, by = "ID")
                    
                    clustering_info$transition_data <- calcTransProbs(clustering, estGroupSize=allocList$estGroupSize, thin=1, 
                                   printXtable=FALSE, printSd=FALSE, printTogether=FALSE,
                                   plotPaths = FALSE, 
                                   plotPathsForE = FALSE)
                    
                    clustering_info$hetero <- calcVariationDMC(clustering, thin = 1, maxi = 50, printAllTogether = FALSE)
                    
                    clustering_info$model <- clustering
                  }
  )

  clustering_info$clustered_data <- clustering_info$clustered_data %>%
    mutate(cluster = as.factor(cluster))

  clustering_info$model_type = cluster_method

  
 
  
  clustering_info$num_clusts <- num_clusts
  
  clustering_info <- fix_cluster_orders(clustering_info)
  
  clustering_info
}

# A function to compare two clusters: Adjusted Rand Index, Jaccard Coefficient, Entropy, Minkowski
compare_clustering <- function(clustering_info_1, clustering_info_2) {
  clustering_1_memberships <- clustering_info_1$clustered_data %>% 
    group_by(ID) %>% 
    sample_n(1)
  
  clustering_1_memberships <- clustering_1_memberships$cluster
  
  clustering_2_memberships <- clustering_info_2$clustered_data %>% 
    group_by(ID) %>% 
    sample_n(1)
  
  clustering_2_memberships <- clustering_2_memberships$cluster
  partition_1 <- new("Partition", clustering_1_memberships)
  partition_2 <- new("Partition", clustering_2_memberships)
  rand_index <- adjustedRandIndex(partition_1, partition_2)
  jacc_coef <- jaccardCoefficient(partition_1, partition_2)
  norm_mutual_info <- normalizedMutualInformation(partition_1 ,partition_2)
  var_info <- variationOfInformation(partition_1, partition_2)
  normalized_var_info <- var_info/(min(log(dim(clustering_info_1$clustered_data)[1]),2*log(as.numeric(clustering_info_1$clust_num))))
  
  comparison_table <- data.table(`Adjusted Rand Index` = round(rand_index, 3), 
                                 `Jaccard Index` = round(jacc_coef, 3), 
                                 `Normalized Mutual Information` = round(norm_mutual_info, 3),
                                 `Normalized Variation of Information` = round(normalized_var_info, 3))
  rownames(comparison_table) <- NULL
  
  return(comparison_table)
  
}

#Two-way table for cluster labels
clust_table <- function(clustering_info_1, clustering_info_2) {
  clustering_1_memberships <- clustering_info_1$clustered_data %>% 
    group_by(ID) %>% 
    sample_n(1)
  
  first_method <- clustering_1_memberships$cluster
  
  clustering_2_memberships <- clustering_info_2$clustered_data %>% 
    group_by(ID) %>% 
    sample_n(1)
  
  second_method <- clustering_2_memberships$cluster
  
  cluster_table <- as.matrix(table(`First Method` = first_method, `Second Method` = second_method))
  
  return(cluster_table)
}


relevel_clusters <- function(vec, ordering){
  vapply(vec, function(ele){which(ele==ordering)}, FUN.VALUE = numeric(1))
}

fix_cluster_orders <- function(input){
  #counts <- input$clustered_data %>% 
  #  group_by(ID) %>% 
  #  sample_n(1) %>% 
  #  ungroup() %>% 
  #  group_by(cluster) %>% 
  #  count() %>% 
  #  arrange(desc(n))
  input <- input
  maxstate <- input$clustered_data %>%
    group_by(state) %>%
    summarize(duration = sum(duration)) %>%
    arrange(desc(duration)) %>%
    dplyr::select(state) %>% 
    slice(1) %>%
    pull()
  
  counts <- input$clustered_data %>%
    group_by(cluster,state) %>%
    summarize(duration = sum(duration)) %>%
    ungroup() %>%
    group_by(cluster) %>%
    mutate(total = sum(duration)) %>%
    mutate(prop = duration/total) %>%
    filter(state == maxstate) %>%
    arrange(desc(prop))
    
  
  
  correct_order <- counts$cluster
  if(length(correct_order) < length(unique(input$clustered_data$cluster))) {
    correct_order <- c(as.numeric(as.character(correct_order)), setdiff(1:length(unique(input$clustered_data$cluster)), as.numeric(levels(correct_order))[correct_order]))
  }
  fixed_cluster_assigments <- relevel_clusters(input$clustered_data$cluster, correct_order)

  cluster_data_fixed <- input$clustered_data %>% 
    mutate(cluster = as.factor(fixed_cluster_assigments))
  cluster_data_fixed <- cluster_data_fixed
  if(input$model_type %in% c("MMM", "DMCC", "MCC")){
    names(cluster_data_fixed)[(ncol(cluster_data_fixed)-input$num_clusts + 1):ncol(cluster_data_fixed)] <- paste0("post_prob_cluster_", 
                                                                relevel_clusters(1:input$num_clusts, c(as.numeric(as.character(correct_order)), setdiff(1:input$num_clusts, as.numeric(levels(correct_order))[correct_order]))))
    
    if(input$model_type %in% c("MMM")) {
      fixed_trans_probs <- list()
      fixed_initial_probs <- list()
      for(i in 1:length(correct_order)){
        fixed_trans_probs[[paste0("Cluster ", i)]] <- input$model$transition_probs[[correct_order[i]]]
        fixed_initial_probs[[paste0("Cluster ", i)]] <- input$model$initial_probs[[correct_order[i]]]
      }
  
      input$model$transition_probs <- fixed_trans_probs
      input$model$initial_probs <- fixed_initial_probs
    } else {
      input$transition_data$estTransProb <- input$transition_data$estTransProb[,,correct_order]
      input$transition_data$estTransProbSd <- input$transition_data$estTransProbSd[,,correct_order]
      
      if(input$model_type == "DMCC") {
        # names not swapped but also never used
        input$hetero$var_e <- input$hetero$var_e[,,correct_order]
        input$hetero$het <- input$hetero$het[,correct_order]
        input$hetero$hetsd <- input$hetero$hetsd[,correct_order]
        input$model$eta_m <- input$model$eta_m[,correct_order]
        input$model$xi_h_m <- input$model$xi_h_m[,,correct_order,]
        input$model$e_h_m <- input$model$e_h_m[,,correct_order,]
      } else {
        input$model$eta.m <- input$model$eta.m[correct_order,]
        input$model$xi.m <- input$model$xi.m[,,,correct_order]
      }
    }
  }
  
  
  input$clustered_data <- cluster_data_fixed
  input
}