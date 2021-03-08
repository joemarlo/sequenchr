
#' Compute CH index and silhouette width
#'
#' Computes the raw and normalized Calinski-Harabasz index and silhouette width for various number of clusters.
#'
#' @param dist_matrix a distance matrix
#' @param cluster_model a clustering model such as the output from hclust
#' @param k_min the minimum number of clusters to test
#' @param k_max the maximum number of clusters to test
#'
#' @return tibble
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
#' cluster_model <- hclust(d = as.dist(dist_matrix), method = 'ward.D2')
#'
#' cluster_stats(
#'  dist_matrix = as.dist(dist_matrix),
#'  cluster_model = cluster_model,
#'  k_min = 2,
#'  k_max = 5
#' )
#'
cluster_stats <- function(dist_matrix, cluster_model, k_min, k_max){
  all_stats <- base::lapply(k_min:k_max, function(k){
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


#' Convert output of TraMineR::seqdef to tidy dataframe
#'
#' Converts the output matrix of TraMineR::seqdef into a tidy tibble
#'
#' @param sequence_data an object created from TraMineR::seqdef
#'
#' @return tibble
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
#' tidy_sequence_data(mvad.seq)
tidy_sequence_data <- function(sequence_data){
  tidy_df <- sequence_data %>%
    dplyr::as_tibble() %>%
    stats::setNames(1:ncol(sequence_data)) %>%
    dplyr::mutate(sequenchr_seq_id = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = base::setdiff(colnames(.), "sequenchr_seq_id"),
                        names_to = 'period', values_to = 'state') %>%
    dplyr::mutate(period = as.numeric(period))

  return(tidy_df)
}


#' Shannon entropy of a vector
#'
#' Calculates the Shannon entropy of a categorical vector
#'
#' @param x a categorical vector
#'
#' @return numeric
#' @export
#'
#' @examples
#' shannon_entropy(c('A', 'A', 'B', 'C', 'A', 'C'))
shannon_entropy <- function(x){

  tab <- table(x)
  prop <- tab / sum(tab)
  entropy <- -sum(ifelse(prop > 0, prop * log(prop, base = 2), 0))

  return(entropy)
}


#' Label clusters consistently
#'
#' Returns cluster labels that match the clusters in sequenchr::plot_dendrogram branches left-to-right.
#'
#' @param .model an hclust model
#' @param k the number of clusters
#'
#' @return a factor
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
#' cluster_model <- hclust(d = as.dist(dist_matrix), method = 'ward.D2')
#'
#' label_clusters(cluster_model, k = 5)
label_clusters <- function(.model, k){

  if (isFALSE(inherits(.model, "hclust"))) stop('.model must be a hclust model produced by stats::hclust or fastcluster::hclust')

  # raw cluster labels
  cluster_labels_raw <- stats::cutree(.model, k = k)

  # reorder clusters to match dendrogram left to right
  cluster_to_dend_mapping <- dplyr::tibble(cluster = cluster_labels_raw[.model$order]) %>%
    tidyr::nest(data = -cluster) %>%
    dplyr::mutate(cluster_dend = dplyr::row_number()) %>%
    tidyr::unnest(data) %>%
    dplyr::distinct()
  cluster_sorted <- dplyr::tibble(cluster = cluster_labels_raw) %>%
    dplyr::left_join(cluster_to_dend_mapping, by = 'cluster') %>%
    dplyr::pull(cluster_dend)

  # add label
  cluster_ns <- base::table(cluster_sorted)
  cluster_labels <- factor(
    cluster_sorted,
    labels = paste("Cluster", 1:k, " | n = ", cluster_ns)
  )

  return(cluster_labels)
}


#' Create a bootstrap html table
#'
#' Similar to knitr::kable and kableExtra::kable_styling without the dependencies.
#'
#' @param .rownames vector of length three
#' @param .values vector of length three
#'
#' @return html code for a 2x3 table
#' @export
#'
#' @examples
#' my_names <- c('n sequences', 'n unique seqences', 'n periods')
#' my_values <- c(50, 20, 10)
#' bootstrap_table(my_names, my_values)
bootstrap_table <- function(.rownames, .values){

  if (length(.rownames) != 3 | length(.values) !=3) stop(".rownames and .values must be length 3")

  # create html table
  html_table <- paste0(
    '<table class="table table-hover table-condensed" style="margin-left: auto; margin-right: auto;">
    <tbody>
    <tr>
    <td style="text-align:left;">',
    .rownames[1],
    '</td>
    <td style="text-align:right;">',
    .values[1],
    '</td>
    </tr>
    <tr>
    <td style="text-align:left;">',
    .rownames[2],
    '</td>
    <td style="text-align:right;">',
    .values[2],
    '</td>
    </tr>
    <tr>
    <td style="text-align:left;">',
    .rownames[3],
    '</td>
    <td style="text-align:right;">',
    .values[3],
    '</td>
    </tr>
    </tbody>
    </table>
  '
  )

  return(html_table)
}

#' Calculate the transition matrix of states
#'
#' Calculate the transition rates between states. If cluster_labels are provided then transition rates are calculated by cluster.
#'
#' @param seq_def_tidy a tidy tibble generated from sequenchr::tidy_sequence_data
#' @param period_min the minimum period to include in the transition calculation
#' @param period_max the maxmium period to include in the transition calculation
#' @param cluster_labels optional. A vector of cluster assignments
#'
#' @return a tidy tibble containing the transition rates (per cluster)
#' @export
#'
#' @seealso \code{\link{plot_transition_matrix}}
#'
#' @examples
#' library(TraMineR)
#' data(mvad)
#' seqstatl(mvad[, 17:86])
#' mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school",
#'                    "training")
#' mvad.labels <- c("employment", "further education", "higher education",
#'                  "joblessness", "school", "training")
#' mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet,
#'                    labels = mvad.labels, xtstep = 6)
#' seq_def_tidy <- tidy_sequence_data(mvad.seq)
#'
#' trans_tidy <- transition_matrix(seq_def_tidy)
#'
#' dist_matrix <- TraMineR::seqdist(seqdata = mvad.seq, method = "DHD")
#' cluster_model <- hclust(d = as.dist(dist_matrix), method = 'ward.D2')
#' cluster_labels <- stats::cutree(cluster_model, k = 5)
#'
#' trans_tidy <- transition_matrix(seq_def_tidy, cluster_labels = cluster_labels)
transition_matrix <- function(seq_def_tidy, period_min = min(seq_def_tidy$sequenchr_seq_id), period_max = max(seq_def_tidy$sequenchr_seq_id), cluster_labels = NULL){

  # TODO: issue here that labels should be comprehensive regardless of period

  # add cluster labels, if exists
  if (!is.null(cluster_labels)) {
    seq_def_tidy <- dplyr::left_join(
      x = seq_def_tidy,
      y = data.frame(sequenchr_seq_id = 1:length(cluster_labels),
                     cluster = cluster_labels),
      by = 'sequenchr_seq_id'
    )
  } else {
    seq_def_tidy$cluster <- 'dummy'
  }

  # filter to just the periods of interest
  filtered_data <- seq_def_tidy %>%
    dplyr::filter(period >= period_min,
                  period <= period_max) %>%
    dplyr::mutate(state = as.character(state))

  # calculate the transition rate per cluster
  cluster_rates <- filtered_data %>%
    dplyr::group_by(cluster) %>%
    dplyr::group_split() %>%
    lapply(X = ., FUN = function(df){

      # add NA filler rows after each group before calculating transition matrix
      # this prevents end of day looping back to beginning of day for next group
      df_NA <- df %>%
        dplyr::group_by(sequenchr_seq_id) %>%
        dplyr::group_split() %>%
        lapply(X = ., FUN = function(df){
          dplyr::add_row(df)
        }) %>%
        dplyr::bind_rows()

      # calculate transition matrix
      n <- nrow(df_NA)
      TRATE_mat <- base::table(data.frame(previous = df_NA$state[1:(n-1)],
                                          current = df_NA$state[2:n]))
      TRATE_mat <- TRATE_mat / sum(TRATE_mat)

      # ensure matrix contains all the states (b/c above filters may remove some)
      unique_states <- unique(seq_def_tidy$state) %>% as.vector()
      TRATE_filled <- tidyr::crossing(previous = unique_states, current = unique_states) %>%
        dplyr::left_join(dplyr::as_tibble(TRATE_mat),
                         by = c('previous', 'current')) %>%
        tidyr::replace_na(list(n = 0)) %>%
        tidyr::pivot_wider(names_from = previous, values_from = n)

      # add cluster name
      TRATE_filled$cluster <- df$cluster[[1]]

      return(TRATE_filled)
    }) %>%
    dplyr::bind_rows()

  # tidy the data
  TRATE_tidy <- tidyr::pivot_longer(cluster_rates,
                                    cols = -c('current', 'cluster'),
                                    names_to = "previous",
                                    values_to = "n"
  )

  return(TRATE_tidy)
}
