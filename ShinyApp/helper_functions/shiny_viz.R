library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(dendextend)
library(igraph)
library(purrr)

# A function to visualize number of IDs by cluster (proportion or stacked bar)
plot_num_ids <- function(input) {
  input <- input
  num_ids_data <- input$clustered_data %>%
    group_by(cluster) %>%
    summarize(num_ids = n_distinct(ID), prop_num_ids = num_ids / length(unique(input$clustered_data$ID)))

  num_ids_plot <- num_ids_data %>%
    ggplot() +
    geom_col(aes(x = cluster, y = prop_num_ids), show.legend = FALSE) +
    ylim(0, 1) +
    labs(title = "Cluster Distribution", x = "Cluster", y = "Proportion") +
    #scale_fill_viridis_d(option = "C") +
    theme_minimal()

  return(num_ids_plot)
}

# A function to visualize state distribution by cluster
plot_state_dist <- function(input) { 
  state_dist_plot <- input$clustered_data %>%
    group_by(cluster,state) %>%
    summarize(duration = sum(duration)) %>%
    ungroup() %>%
    group_by(cluster) %>%
    mutate(total = sum(duration)) %>%
    mutate(prop = duration/total) %>%
    ggplot() +
    geom_bar(aes(x = cluster, fill = state,y = prop), stat='identity') +
    labs(title = "State Distribution", x = "Cluster", y = "Proportion") +
    scale_fill_manual(name = "State", values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")) +
    theme_minimal()

  return(state_dist_plot)
}




# A function to visualize total duration distribution by cluster (boxplot)
plot_total_duration <- function(input) {
  total_duration_data <- input$clustered_data %>%
    group_by(cluster, ID) %>%
    summarize(total_duration = sum(duration))

  total_duration_plot <- total_duration_data %>%
    ggplot() +
    geom_boxplot(aes(x = cluster, y = total_duration, fill = cluster), show.legend = FALSE) +
    labs(title = "Total Duration Distribution", x = "Cluster", y = "Duration") +
    scale_fill_viridis_d(option = "C") +
    theme_minimal()

  return(total_duration_plot)
}

# A function to visualize average duration within a state
plot_avg_state_duration <- function(input,format=NULL) {
  data <- input$clustered_data
  if(!is.null(format))
  if(format == 'compress'){
    data <- data %>%
      select(-start_time,-duration) %>%
      rename(start_time = start_time_comp,duration = duration_comp)
  }
  
  avg_state_duration_data <- data %>%
    group_by(cluster, state) %>%
    summarize(avg_state_duration = mean(duration))
  
  avg_state_duration_plot <- avg_state_duration_data %>%
    ggplot() +
    geom_col(aes(x = cluster, y = avg_state_duration, fill = state), position = "dodge") +
    labs(title = "Average Duration By State", x = "Cluster", y = "Duration") +
    scale_fill_manual(name = "State", values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")) +
    theme_minimal()
  
  return(avg_state_duration_plot)
} 

# A function to visualize average number of transitions
plot_avg_num_spells <- function(input) {
  avg_num_trans_data <- input$clustered_data %>%
    add_count(cluster, state) %>%
    group_by(cluster) %>%
    mutate(num_ids = n_distinct(ID), avg_num_trans = n/num_ids)

  avg_num_trans_plot <- avg_num_trans_data %>%
    ggplot() +
    geom_col(aes(x = cluster, y = avg_num_trans, fill = state), position = "dodge") +
    labs(title = "Average Count of Spells By State", x = "Cluster", y = "Count of Spells") +
    scale_fill_manual(name = "State", values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")) +
    theme_minimal()

  return(avg_num_trans_plot)
}


# A function to visualize sample trajectories
plot_sample_traj <- function(input,format=NULL) {
  data <- input$clustered_data
  if(!is.null(format))
    if(format == 'compress'){
    data <- data %>%
      select(-start_time,-duration) %>%
      rename(start_time = start_time_comp,duration = duration_comp)
  }
  
  sample_by_clust <- data %>% 
    select(ID, cluster) %>% 
    unique() %>% 
    count(cluster) %>% 
    mutate(n = pmin(n, 5))
  

  K <- length(levels(data$cluster))
  data %>%
    select(ID, cluster) %>%
    unique() %>%
    nest(-cluster) %>% 
    left_join(sample_by_clust, by = "cluster") %>% 
    mutate(sample = map2(data, n, sample_n)) %>% 
    unnest(sample) %>% 
    group_by(cluster) %>% 
    mutate(y = 1:n()) %>%
    inner_join(data) %>%
    mutate(cluster2 = paste('Cluster ',cluster)) %>%
    ggplot(aes(xmin = start_time, xmax = start_time + duration, ymin = y, ymax = y + 1)) + 
    geom_rect(aes(fill = state),color='darkgrey') + 
    facet_grid(~cluster2) +
    labs(title = "Sample Sequence Trajectories", x = "Time") +
    scale_fill_manual(name = "State", values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")) +
    theme_minimal()+
    theme(axis.text.y=element_blank(),
    axis.ticks=element_blank()) 
} 





# A function to visualize state proportions over time
plot_prop_time <- function(input,format=NULL) {
    #incorporate format into this...  
    longdat <- data_to_long(input$clustered_data,format)
    tmp_full <- as.data.frame(expand.grid(unique(longdat$time), unique(longdat$state),unique(longdat$cluster)))
    names(tmp_full) <- c('time','state','cluster')
  
    plot_prop_time <- longdat %>%
    count(cluster,time,state) %>%
    left_join(tmp_full,.) %>%
    mutate(n = ifelse(is.na(n), 0, n)) %>%
    group_by(cluster,time) %>%
    mutate(relfreq = n/sum(n)) %>%
    ggplot(aes(x = factor(time), y = relfreq, fill = state)) + 
    geom_bar(stat='identity',width=1) + 
    facet_grid(~cluster) +
    labs(title = "State Distribution over Time", x = "Time",y = 'Proportion') +
    scale_fill_manual(name = "State", values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")) +
    theme_minimal()+
    theme(axis.text.x=element_blank(),
      axis.ticks=element_blank()) 
    
    return(plot_prop_time)
} 






# A function to visualize dendrogram for hierarchical clustering
plot_dendro <- function(input) {
  dend <- as.dendrogram(input$cluster_output) %>% dendextend::set("branches_k_color", k = input$num_clusts)
  dend <- plot(dend, leaflab = "none")
  
  return(dend)
}


plot_mark_heat_map <- function(input){
  # par(mfrow = c(ceiling(input$model$n_clusters/2), 2))
  # par(mar = c(1,0,1,0) + 0.1)
  
  clust_info <- bind_rows(lapply(1:length(input$model$transition_probs), 
                                 function(index) {
                                   trans_data_frame <- melt(input$model$transition_probs[[index]])
                                   trans_data_frame[["Cluster"]] <- paste0("Cluster ", index)
                                   trans_data_frame
                                 }))
  
  ggplot(clust_info, aes(to, from)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_viridis(name = "Transition \nProbabilites",begin=0,end=1) +
    labs(x = "To", y = "From") +
    facet_wrap(~Cluster) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Estimated Transition Probabilities")
  
}

plot_mark_post_probs <- function(input){
  input$clustered_data %>% 
    group_by(ID) %>% 
    sample_n(1) %>% 
    ungroup() %>% 
    gather((ncol(input$clustered_data) - input$num_clusts + 1):ncol(input$clustered_data), key = "cluster_prob", value = "prob") %>% 
    select(-cluster_prob) %>% 
    mutate(cluster = paste0("Cluster ", cluster)) %>% 
    group_by(cluster, ID) %>% 
    summarize(posterior_prob = max(prob)) %>% 
    ggplot()+
    geom_histogram(aes(x = posterior_prob, y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), binwidth = .05, boundary = 0)+
    facet_wrap(~cluster, nrow = ceiling(input$num_clusts/3))+
    theme_minimal()+
    labs(x = "Posterior Cluster Probability", y= "Proportion") +
    ggtitle("Posterior Probability of Belonging to Most Probable Cluster")
}


plot_dissim_matrix <- function(input) {
  just_clust <- (input$clustered_data %>% 
    group_by(ID) %>% 
    sample_n(1) %>% 
    ungroup())$cluster
  just_clust <- as.numeric(just_clust)
  
  
  dissplot(input$dist_matrix, 
           labels = just_clust, 
           options = list(col = viridis::viridis(100),
                          main = "Dissimilarity Between Observations and Clusters",
                          silhouettes = TRUE))
}

plot_mcc_trans_heat_map <- function(input){
  # par(mfrow = c(ceiling(input$model$n_clusters/2), 2))
  # par(mar = c(1,0,1,0) + 0.1)
  
  clust_info <- bind_rows(lapply(1:(dim(input$transition_data$estTransProb)[3]), 
                                 function(index) {
                                   trans_mat <- input$transition_data$estTransProb[,,index]
                                   dimnames(trans_mat) <- list(from = sort(unique(input$clustered_data$state)), to = sort(unique(input$clustered_data$state)))
                                   trans_data_frame <- melt(trans_mat)
                                   trans_data_frame[["Cluster"]] <- paste0("Cluster ", index)
                                   trans_data_frame
                                 }))
  
  ggplot(clust_info, aes(to, from)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_viridis(name = "Transition \nProbabilites",begin=0,end=1) +
    labs(x = "To", y = "From") +
    facet_wrap(~Cluster) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Average Transition Probabilities across MCMC Iterations")
  
}

plot_mcc_sd_heat_map <- function(input){
  # par(mfrow = c(ceiling(input$model$n_clusters/2), 2))
  # par(mar = c(1,0,1,0) + 0.1)

  clust_info <- bind_rows(lapply(1:(dim(input$transition_data$estTransProbSd)[3]), 
                                 function(index) {
                                   trans_mat <- input$transition_data$estTransProbSd[,,index]
                                   dimnames(trans_mat) <- list(from = sort(unique(input$clustered_data$state)), to = sort(unique(input$clustered_data$state)))
                                   trans_data_frame <- melt(trans_mat)
                                   trans_data_frame[["Cluster"]] <- paste0("Cluster ", index)
                                   trans_data_frame
                                 }))
  
  ggplot(clust_info, aes(to, from)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_viridis(name = "SD of Transition \nProbabilites",begin=0,end=1) +
    labs(x = "To", y = "From") +
    facet_wrap(~Cluster) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Standard Deviation of Transition Probabilities across MCMC Iterations")
  
}

plot_mark_chain <- function(input) {
  model_type <- input$model_type
  par(mfrow = c(ceiling(input$num_clusts/2), 2))
  par(mar = c(1,0,1,0) + 0.1)
  par(oma = c(0,0,2,0))
  lapply(1:input$num_clusts, function(i){
    if(model_type == "MMM") {
      transM <- input$model$transition_probs[[i]]
    } else {
      transM <- input$transition_data$estTransProb[,,i]
    }
    transM[transM < .0001] <- 0
    diag(transM) <- 0
    edges <- transM
    edges[edges > 0] <- 1
    
    transitions <- transM
    transitions <- t(transitions)[t(transitions) > 0]
    if(model_type == "MMM") {
      edge.width <- transitions * (7 / max(sapply(input$model$transition_probs, function(trans_probs){
        transM <- trans_probs
        transM[transM < .0001] <- 0
        diag(transM) <- 0
        max(transM)
      })))
    } else {
      edge.width <- transitions * (7 / max(sapply(1:input$num_clusts, function(index){
        transM <- input$transition_data$estTransProb[,,index]
        transM[transM < .0001] <- 0
        diag(transM) <- 0
        max(transM)
      })))
    }
    
    
    g1 <<- graph.adjacency(transM > 0, mode = "directed")
    glayout <- layout_in_circle(g1)
    plot.igraph(g1,  layout = glayout,
                edge.width = edge.width,
                edge.curved = TRUE,
                vertex.label = (function(){
                  if(model_type == "MMM") {c(names(transM))} else {sort(unique(input$clustered_data$state))}})(),
                edge.label = round(transitions, 4),
                vertex.size = 50,
                edge.arrow.size = .9,
                main = paste("Cluster", i))
  })
  mtext("Estimated Markov Chains", outer = TRUE, cex = 1.5)
}


plot_variance <- function(input) {
  clust_info <- bind_rows(lapply(1:(dim(input$hetero$var_e)[3]), 
                                 function(index) {
                                   var_mat <- input$hetero$var_e[,,index]
                                   dimnames(var_mat) <- list(from = sort(unique(input$clustered_data$state)), to = sort(unique(input$clustered_data$state)))
                                   var_data_frame <- melt(var_mat)
                                   var_data_frame[["Cluster"]] <- paste0("Cluster ", index)
                                   var_data_frame
                                 }))
  
  ggplot(clust_info, aes(to, from)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_viridis(name = "Variance of \nTransition Probabilites",begin=0,end=1) +
    labs(x = "To", y = "From") +
    facet_wrap(~Cluster) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Posterior Expectation of Variance of Transition Probabilities")
}

plot_unob_het <- function(input) {
  het_df <- as.data.frame(input$hetero$het) %>% 
    mutate(State = sort(unique(input$clustered_data$state))) %>% 
    gather(key = "Cluster", value = "Het", 1:input$num_clusts) %>% 
    mutate(Cluster = str_sub(Cluster, 7))
  hetsd_df <- as.data.frame(input$hetero$hetsd) %>% 
    mutate(State = sort(unique(input$clustered_data$state))) %>% 
    gather(key = "Cluster", value = "SD", 1:input$num_clusts) %>% 
    mutate(Cluster = str_sub(Cluster, 7))
  
  het_df %>% 
    left_join(hetsd_df, by = c("State", "Cluster")) %>% 
    mutate(Cluster = paste0("Cluster ", Cluster)) %>% 
    ggplot() + 
    geom_col(aes(x = State, y = Het, fill = SD), position = "dodge") +
    scale_fill_viridis(name = "Standard \nDeviation") +
    ylab("Within Cluster Heterogeneity") +
    theme_minimal() +
    coord_flip() +
    facet_wrap(~Cluster) +
    ggtitle("Posterior Exepectation of Unobserved Within Cluster Heterogeneity")
}

plot_group_size_trace <- function(input) {
  if(input$model_type == "DMCC"){
    eta_m <- input$model$eta_m
    dimnames(eta_m) <- list(iteration = 1:dim(eta_m)[1],
                            cluster = as.character(1:dim(eta_m)[2]))
  } else {
    eta_m <- input$model$eta.m
    dimnames(eta_m) <- list(cluster = as.character(1:dim(eta_m)[1]),
                            iteration = 1:dim(eta_m)[2])
    
  }
  
  melt(eta_m) %>%  
    ggplot() +
    geom_line(aes(x = iteration, y = value, color = as.character(cluster))) +
    ylab("Group Sizes") +
    xlab("Iteration") +
    ggtitle("Trace Plot of Estimated Cluster Sizes") +
    theme_minimal() +
    scale_color_discrete(name= "Cluster")

}

plot_diag_trace <- function(input) {
  if(input$model_type == "DMCC") {
    ksi_arr <- input$model$xi_h_m
    dimnames(ksi_arr) <- list(from = sort(unique(input$clustered_data$state)),
                              to = sort(unique(input$clustered_data$state)),
                              cluster = 1:dim(ksi_arr)[3],
                              iteration = 1:dim(ksi_arr)[4])
  } else {
    ksi_arr <- input$model$xi.m
    dimnames(ksi_arr) <- list(iteration = 1:dim(ksi_arr)[1],
                              from = sort(unique(input$clustered_data$state)),
                              to = sort(unique(input$clustered_data$state)),
                              cluster = 1:dim(ksi_arr)[4])
  }
  
  melt(ksi_arr) %>% 
    filter(from == to) %>% 
    ggplot(aes(x = iteration, y = value, color = from)) +
    geom_line() +
    facet_wrap(~cluster) +
    scale_color_discrete(name = "Diagonal") +
    theme_minimal() +
    ggtitle("Trace Plot of Estimated Diagonal Entries of Transition Matrices")
}

plot_group_het_trace <- function(input) {
  e_array <- input$model$e_h_m
  dimnames(e_array) <- list(from = sort(unique(input$clustered_data$state)),
                            to = sort(unique(input$clustered_data$state)),
                            cluster = 1:dim(e_array)[3],
                            iteration = 1:dim(e_array)[4])
  
  melt(e_array) %>% 
    filter(from == to) %>% 
    ggplot(aes(x = iteration, y = value, color = from)) +
    geom_line() +
    facet_wrap(~cluster) +
    scale_color_discrete(name = "Diagonal") +
    theme_minimal() +
    ggtitle("Trace Plot of Estimated Within Group Heterogeneity on Diagonals")
}