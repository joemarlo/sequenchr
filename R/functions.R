
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


#' Convert output of TraMineR::seqdef to tidy dataframe
#'
#' Converts the output matrix of TraMineR::seqdef into a tidy tibble
#'
#' @param sequence_data
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
    tidyr::pivot_longer(cols = base::setdiff(colnames(.), "sequenchr_seq_id")) %>%
    dplyr::mutate(period = as.numeric(name)) %>%
    dplyr::select(-name)

  return(tidy_df)
}
