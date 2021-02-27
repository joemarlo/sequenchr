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

  # TODO: write error handling for providing cluster_assignments but not n_col_facets
  # TODO: allow renaming of title and xy labels?

  if (is.null(cluster_assignments)){

    # plot the regular sequences without clustering
    p <- seq_def_tidy %>%
      dplyr::group_by(sequenchr_seq_id) %>%
      dplyr::mutate(entropy = DescTools::Entropy(table(value))) %>%
      dplyr::ungroup() %>%
      ggplot2::ggplot(ggplot2::aes(x = period, y = stats::reorder(sequenchr_seq_id, entropy), fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_manual(values = color_mapping) +
      ggplot2::scale_y_discrete(labels = NULL, breaks = NULL) +
      ggplot2::labs(title = "All sequences sorted by entropy",
           x = 'Period',
           y = 'Sequence',
           fill = NULL)

  } else {

    # plot the sequences with clusters
    p <- dplyr::tibble(cluster = cluster_assignments,
                sequenchr_seq_id = 1:length(cluster_assignments)) %>%
      dplyr::right_join(seq_def_tidy, by = 'sequenchr_seq_id') %>%
      dplyr::group_by(sequenchr_seq_id) %>%
      dplyr::mutate(entropy = DescTools::Entropy(table(value))) %>%
      dplyr::ungroup() %>%
      ggplot2::ggplot(ggplot2::aes(x = period, y = stats::reorder(sequenchr_seq_id, entropy), fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_manual(values = color_mapping) +
      ggplot2::scale_y_discrete(labels = NULL, breaks = NULL) +
      ggplot2::facet_wrap(~cluster, scales = 'free_y', ncol = n_col_facets) +
      ggplot2::labs(title = "All sequences by cluster sorted by entropy",
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
      ggplot2::ggplot(ggplot2::aes(x = period, fill = value)) +
      ggplot2::geom_bar(width = 1) +
      ggplot2::scale_fill_manual(values = color_mapping) +
      ggplot2::labs(title = "State distributions",
           x = 'Period',
           y = 'Frequency',
           fill = NULL)

  } else {

    # plot with clustering
    p <- dplyr::tibble(cluster = cluster_assignments,
                sequenchr_seq_id = 1:length(cluster_assignments)) %>%
      dplyr::right_join(seq_def_tidy, by = 'sequenchr_seq_id') %>%
      ggplot2::ggplot(ggplot2::aes(x = period, fill = value)) +
      ggplot2::geom_bar(width = 1) +
      ggplot2::scale_fill_manual(values = color_mapping) +
      ggplot2::facet_wrap(~cluster, scales = 'free_y', ncol = n_col_facets) +
      ggplot2::labs(title = "State distributions by cluster",
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
      dplyr::count(value, period) %>%
      dplyr::group_by(period) %>%
      dplyr::filter(n == max(n)) %>%
      ggplot2::ggplot(ggplot2::aes(x = period, y = n, fill = value)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = color_mapping) +
      ggplot2::labs(title = "Modal activity per period",
           caption = "Ties are shown as stacked bars",
           x = "Period",
           y = 'Frequency',
           fill = NULL)
  } else {
    # plot with cluster
    p <- dplyr::tibble(cluster = cluster_assignments,
                sequenchr_seq_id = 1:length(cluster_assignments)) %>%
      dplyr::right_join(seq_def_tidy, by = 'sequenchr_seq_id') %>%
      dplyr::count(cluster, value, period) %>%
      dplyr::group_by(cluster, period) %>%
      dplyr::filter(n == max(n)) %>%
      ggplot2::ggplot(ggplot2::aes(x = period, y = n, fill = value)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = color_mapping) +
      ggplot2::facet_wrap(~cluster, scales = 'free_y', ncol = n_col_facets) +
      ggplot2::labs(title = "Modal activity per period by cluster",
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
    dplyr::mutate(index = dplyr::row_number()) %>%
    ggplot2::ggplot(ggplot2::aes(x=1, y = stats::reorder(value, -index), fill = value)) +
    ggplot2::geom_tile(color = 'white', size = 3) +
    ggplot2::scale_fill_manual(values = color_mapping) +
    ggplot2::scale_x_continuous(labels = NULL) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(legend.position = 'none')

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
  dend <- stats::as.dendrogram(cluster_model) %>%
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
    dplyr::filter(col != 'grey50') %>%
    dplyr::group_by(col) %>%
    dplyr::summarize(x = mean(x), .groups = 'drop') %>%
    dplyr::arrange(x) %>%
    dplyr::mutate(label = paste0("Cluster ", 1:k))

  # plot the dendrograms
  p <- ggd1$segments %>%
    ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                 color = ggd1$segments$col, linetype = ggd1$segments$linetype,
                 lwd = 0.9, alpha = 0.7) +
    ggplot2::scale_x_continuous(labels = cluster_labels$label,
                       breaks = cluster_labels$x) +
    ggplot2::scale_y_continuous(labels = scales::comma_format()) +
    ggplot2::labs(title = "Dendrogram",
         subtitle = 'Helpful subtitle goes here',
         x = NULL,
         y = NULL) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank(),
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
#' #TODO
plot_transition_matrix <- function(transition_matrix){

  # TODO: issue here that labels should be comprehensive regardless of period
  # TODO: add clustering
  p <- transition_matrix %>%
    tidyr::pivot_longer(cols = -current, names_to = "previous", values_to = "n") %>%
    ggplot2::ggplot(ggplot2::aes(x = previous, y = current, fill = n, label = round(n, 3))) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(color = 'grey90') +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(title = "Transition matrix",
         subtitle = "A helpful subtitle",
         x = "\nFrom state",
         y = 'To state',
         fill = 'Transition rate') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  return(p)
}
