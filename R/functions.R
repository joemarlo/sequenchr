plot_sequence_index <- function(seq_def_tidy, color_mapping, cluster_assignments = NULL, n_col_facets = NULL){
  
  if (is.null(cluster_assignments)){
    
    # plot the regular sequences without clustering
    p <- seq_def_tidy %>% 
      group_by(sequenchr_seq_id) %>% 
      mutate(entropy = DescTools::Entropy(table(value))) %>%
      ungroup() %>% 
      ggplot(aes(x = period, y = reorder(sequenchr_seq_id, entropy), fill = value)) +
      geom_tile() +
      scale_fill_manual(values = color_mapping) +
      scale_y_discrete(labels = NULL, breaks = NULL) +
      labs(title = "All sequences sorted by entropy",
           x = 'Period',
           y = 'Sequence',
           fill = NULL)
    
  } else {
    
    # plot the sequences with clusters
    p <- tibble(cluster = cluster_assignments,
                sequenchr_seq_id = 1:length(cluster_assignments)) %>%
      right_join(seq_def_tidy, by = 'sequenchr_seq_id') %>%
      group_by(sequenchr_seq_id) %>%
      mutate(entropy = DescTools::Entropy(table(value))) %>%
      ungroup() %>%
      ggplot(aes(x = period, y = reorder(sequenchr_seq_id, entropy), fill = value)) +
      geom_tile() +
      scale_fill_manual(values = color_mapping) +
      scale_y_discrete(labels = NULL, breaks = NULL) +
      facet_wrap(~cluster, scales = 'free_y', ncol = n_col_facets) +
      labs(title = "All sequences by cluster sorted by entropy",
           x = 'Period',
           y = 'Sequence',
           fill = NULL)
  }
  
  return(p)
}


plot_state <- function(seq_def_tidy, color_mapping, cluster_assignments = NULL, n_col_facets = NULL){
  
  if (is.null(cluster_assignments)){
    
    # plot without clustering
    p <- seq_def_tidy %>% 
      ggplot(aes(x = period, fill = value)) +
      geom_bar(width = 1) +
      scale_fill_manual(values = color_mapping) +
      labs(title = "State distributions",
           x = 'Period',
           y = 'Frequency',
           fill = NULL)
    
    } else {
    
    # plot with clustering
    p <- tibble(cluster = cluster_assignments,
                sequenchr_seq_id = 1:length(cluster_assignments)) %>%
      right_join(seq_def_tidy, by = 'sequenchr_seq_id') %>% 
      ggplot(aes(x = period, fill = value)) +
      geom_bar(width = 1) +
      scale_fill_manual(values = color_mapping) +
      facet_wrap(~cluster, scales = 'free_y', ncol = n_col_facets) +
      labs(title = "State distributions by cluster",
           x = 'Period',
           y = 'Frequency',
           fill = NULL)
  }
  
  return(p)
}


plot_modal <- function(seq_def_tidy, color_mapping, cluster_assignments = NULL, n_col_facets = NULL){
  
  if (is.null(cluster_assignments)){
    
    # plot without clustering
    p <- seq_def_tidy %>% 
      count(value, period) %>% 
      group_by(period) %>% 
      filter(n == max(n)) %>% 
      ggplot(aes(x = period, y = n, fill = value)) +
      geom_col() +
      scale_fill_manual(values = color_mapping) +
      labs(title = "Modal activity per period",
           caption = "Ties are shown as stacked bars",
           x = "Period",
           y = 'Frequency',
           fill = NULL)
  } else {
    # plot with cluster
    p <- tibble(cluster = cluster_assignments,
                sequenchr_seq_id = 1:length(cluster_assignments)) %>%
      right_join(seq_def_tidy, by = 'sequenchr_seq_id') %>% 
      count(cluster, value, period) %>% 
      group_by(cluster, period) %>% 
      filter(n == max(n)) %>% 
      ggplot(aes(x = period, y = n, fill = value)) +
      geom_col() +
      scale_fill_manual(values = color_mapping) +
      facet_wrap(~cluster, scales = 'free_y', ncol = n_col_facets) +
      labs(title = "Modal activity per period by cluster",
           caption = "Ties are shown as stacked bars",
           x = "Period",
           y = 'Frequency',
           fill = NULL)
  }
  
  return(p)
}


plot_legend <- function(color_mapping){
  p <- dplyr::tibble(value = names(color_mapping)) %>%
    mutate(index = row_number()) %>% 
    ggplot(aes(x=1, y = reorder(value, -index), fill = value)) + 
    geom_tile(color = 'white', size = 3) + 
    scale_fill_manual(values = color_mapping) +
    scale_x_continuous(labels = NULL) + 
    labs(x = NULL, y = NULL) + 
    theme(legend.position = 'none')
  
  return(p)
}


plot_dendrogram <- function(cluster_model, k, h = 100){

  # build base dendrogram
  dend <- as.dendrogram(cluster_model) %>% 
    dendextend::set("branches_k_color", k = k) %>% 
    dendextend::set("labels_colors")
  
  # cut off bottom of dendogram for computation performance
  dend <- cut(dend, h = h)$upper
  ggd1 <- dendextend::as.ggdend(dend)
  
  # set dashed line for non-cluster segments
  ggd1$segments$linetype <- 'solid'
  ggd1$segments$linetype[which(is.na(ggd1$segments$col))] <- 'dashed'
  
  # set connecting lines to grey
  ggd1$segments$col[is.na(ggd1$segments$col)] <- 'grey50'
  
  # set the label positions
  cluster_labels <- ggd1$segments %>% 
    filter(col != 'grey50') %>% 
    group_by(col) %>% 
    summarize(x = mean(x), .groups = 'drop') %>% 
    arrange(x) %>% 
    mutate(label = paste0("Cluster ", 1:k))
  
  # plot the dendrograms
  p <- ggd1$segments %>% 
    ggplot() + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
                 color = ggd1$segments$col, linetype = ggd1$segments$linetype,
                 lwd = 0.9, alpha = 0.7) +
    scale_x_continuous(labels = cluster_labels$label,
                       breaks = cluster_labels$x) +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(title = "Dendrogram",
         subtitle = 'Helpful subtitle goes here',
         x = NULL,
         y = NULL) +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = 'none')
  
  return(p)
}


plot_transition_matrix <- function(transition_matrix){

  # TODO: issue here that labels should be comprehensive regardless of period
  # TODO: add clustering
  p <- transition_matrix %>%
    pivot_longer(cols = -current, names_to = "previous", values_to = "n") %>% 
    ggplot(aes(x = previous, y = current, fill = n, label = round(n, 3))) +
    geom_tile() +
    geom_text(color = 'grey90') +
    scale_fill_viridis_c() +
    labs(title = "Transition matrix",
         subtitle = "A helpful subtitle",
         x = "\nFrom state",
         y = 'To state',
         fill = 'Transition rate') +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  
  return(p)
}


cluster_stats <- function(dist_matrix, cluster_model, k_min, k_max){
  all_stats <- lapply(k_min:k_max, function(k){
    c_stats <- fpc::cluster.stats(
      d = dist_matrix,
      clustering = stats::cutree(cluster_model, k = k),
      silhouette = TRUE
    )
    return(dplyr::tibble(k = k, ch = c_stats$ch, silhouette = c_stats$avg.silwidth))
  })
  
  all_stats <- dplyr::bind_rows(all_stats)
  scale_01 <- function(x) (x - min(x)) / diff(range(x))
  all_stats$ch_norm <- scale_01(all_stats$ch)
  all_stats$silhouette_norm <- scale_01(all_stats$silhouette)
  
  return(all_stats)
}


tidy_sequence_data <- function(sequence_data){
  sequence_data %>%
    as_tibble() %>% 
    setNames(1:ncol(sequence_data)) %>% 
    mutate(sequenchr_seq_id = row_number()) %>%
    pivot_longer(cols = setdiff(colnames(.), "sequenchr_seq_id")) %>% 
    mutate(period = as.numeric(name)) %>% 
    dplyr::select(-name)
}
