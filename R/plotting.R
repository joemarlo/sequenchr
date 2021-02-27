#' Generates a sequence index plot
#'
#' @param seq_def_tidy a tidy tibble generated from sequenchr::tidy_sequence_data
#' @param color_mapping a list of named colors where the names match the alphabet of the original sequence data
#' @param cluster_assignments optional. A vector of cluster assignments
#' @param n_col_facets optional. If cluster_assignments is provided then the number of facet columns
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' library(TraMineR)
#' data(mvad)
#' seqstatl(mvad[, 17:86])
#' mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school",
#'                    "training")
#' mvad.labels <- c("employment", "further education", "higher education",
#'                  "joblessness", "school", "training")
#' mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet, # states = mvad.scodes,
#'                    labels = mvad.labels, xtstep = 6)
#' seq_def_tidy <- tidy_sequence_data(mvad.seq)
#' color_mapping <- viridis::viridis_pal()(length(alphabet(mvad.seq)))
#' names(color_mapping) <- alphabet(mvad.seq)
#' plot_sequence_index(seq_def_tidy, color_mapping)
#'
#' dist_matrix <- TraMineR::seqdist(seqdata = mvad.seq, method = "DHD")
#' cluster_model <- fastcluster::hclust(d = as.dist(dist_matrix), method = 'ward.D2')
#' cluster_assignments <- stats::cutree(cluster_model, k = 5)
#' plot_sequence_index(seq_def_tidy, color_mapping, cluster_assignments = cluster_assignments)
plot_sequence_index <- function(seq_def_tidy, color_mapping, cluster_assignments = NULL, n_col_facets = 1){

  # TODO: write error handling for providing cluster_assignemnts but not n_col_facets
  # TODO: allow renaming of title and xy labels?

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

#' Generates a sequence state plot
#'
#' @param seq_def_tidy a tidy tibble generated from sequenchr::tidy_sequence_data
#' @param color_mapping a list of named colors where the names match the alphabet of the original sequence data
#' @param cluster_assignments optional. A vector of cluster assignments
#' @param n_col_facets optional. If cluster_assignments is provided then the number of facet columns
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' library(TraMineR)
#' data(mvad)
#' seqstatl(mvad[, 17:86])
#' mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school",
#'                    "training")
#' mvad.labels <- c("employment", "further education", "higher education",
#'                  "joblessness", "school", "training")
#' mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet, # states = mvad.scodes,
#'                    labels = mvad.labels, xtstep = 6)
#' seq_def_tidy <- tidy_sequence_data(mvad.seq)
#' color_mapping <- viridis::viridis_pal()(length(alphabet(mvad.seq)))
#' names(color_mapping) <- alphabet(mvad.seq)
#' plot_state(seq_def_tidy, color_mapping)
#'
#' dist_matrix <- TraMineR::seqdist(seqdata = mvad.seq, method = "DHD")
#' cluster_model <- fastcluster::hclust(d = as.dist(dist_matrix), method = 'ward.D2')
#' cluster_assignments <- stats::cutree(cluster_model, k = 5)
#' plot_state(seq_def_tidy, color_mapping, cluster_assignments = cluster_assignments)
plot_state <- function(seq_def_tidy, color_mapping, cluster_assignments = NULL, n_col_facets = 1){

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

#' Generates a plot of the modal states
#'
#' @param seq_def_tidy a tidy tibble generated from sequenchr::tidy_sequence_data
#' @param color_mapping a list of named colors where the names match the alphabet of the original sequence data
#' @param cluster_assignments optional. A vector of cluster assignments
#' @param n_col_facets optional. If cluster_assignments is provided then the number of facet columns
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' library(TraMineR)
#' data(mvad)
#' seqstatl(mvad[, 17:86])
#' mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school",
#'                    "training")
#' mvad.labels <- c("employment", "further education", "higher education",
#'                  "joblessness", "school", "training")
#' mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet, # states = mvad.scodes,
#'                    labels = mvad.labels, xtstep = 6)
#' seq_def_tidy <- tidy_sequence_data(mvad.seq)
#' color_mapping <- viridis::viridis_pal()(length(alphabet(mvad.seq)))
#' names(color_mapping) <- alphabet(mvad.seq)
#' plot_modal(seq_def_tidy, color_mapping)
#'
#' dist_matrix <- TraMineR::seqdist(seqdata = mvad.seq, method = "DHD")
#' cluster_model <- fastcluster::hclust(d = as.dist(dist_matrix), method = 'ward.D2')
#' cluster_assignments <- stats::cutree(cluster_model, k = 5)
#' plot_modal(seq_def_tidy, color_mapping, cluster_assignments = cluster_assignments)
plot_modal <- function(seq_def_tidy, color_mapping, cluster_assignments = NULL, n_col_facets = 1){

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


#' Plot the legend
#'
#' Plots just the legend given a list of named oclors
#'
#' @param color_mapping a list of named colors
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' library(TraMineR)
#' data(mvad)
#' seqstatl(mvad[, 17:86])
#' mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school",
#'                    "training")
#' mvad.labels <- c("employment", "further education", "higher education",
#'                  "joblessness", "school", "training")
#' mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet, # states = mvad.scodes,
#'                    labels = mvad.labels, xtstep = 6)
#' color_mapping <- viridis::viridis_pal()(length(alphabet(mvad.seq)))
#' names(color_mapping) <- alphabet(mvad.seq)
#' plot_legend(color_mapping)
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


#' Plot a dendrogram colored by cluster
#'
#' Plots a dedrogram where the colors of the segments represent cluster membership
#'
#' @param cluster_model a clustering model such as the output from fastcluster::hclust
#' @param k the number of clusters
#' @param h the minimum height to plot the segments. A lower height results in decreased performance
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' library(TraMineR)
#' data(mvad)
#' seqstatl(mvad[, 17:86])
#' mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school",
#'                    "training")
#' mvad.labels <- c("employment", "further education", "higher education",
#'                  "joblessness", "school", "training")
#' mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet, # states = mvad.scodes,
#'                    labels = mvad.labels, xtstep = 6)
#' dist_matrix <- TraMineR::seqdist(seqdata = mvad.seq, method = "DHD")
#' cluster_model <- fastcluster::hclust(d = as.dist(dist_matrix), method = 'ward.D2')
#' plot_dendrogram(cluster_model, 5)
plot_dendrogram <- function(cluster_model, k, h = 100){

  # build base dendrogram
  dend <- as.dendrogram(cluster_model) %>%
    dendextend::set("branches_k_color", k = k) %>%
    dendextend::set("labels_colors")

  # cut off bottom of dendogram for computation performance
  dend <- base::cut(dend, h = h)$upper
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


#' Plot a transition matrix
#'
#' Plots a 'heatmap' of a transition matrix
#'
#' @param transition_matrix a transition matrix
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' TODO
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
