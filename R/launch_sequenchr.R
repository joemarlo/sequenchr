#' Launches the sequenchr Shiny app
#'
#' Launches the sequenchr app in the default web browser. Optional covariates_data adds a tab containing plots of the marginal distributions of these covariates.
#'
#' @param sequence_data an object created from TraMineR::seqdef
#' @param covariates_data optional. a dataframe with the same number of rows as sequence_data
#'
#' @return launches the sequenchr Shiny app
#' @export
#' @import shiny
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
#' #launch_sequenchr(mvad.seq)
launch_sequenchr <- function(sequence_data, covariates_data = NULL){

  # quality checks
  is_TraMineR_object <- TraMineR::is.stslist(sequence_data)
  if (isFALSE(is_TraMineR_object)) stop("sequence_data must be a stslist object created from TraMineR::seqdef()")
  if (!is.null(covariates_data)){
    if (isFALSE(inherits(covariates_data, 'data.frame'))) stop("covariates_data must be a dataframe")
    is_all_numeric <- all(apply(covariates_data, 2, function(col) inherits(col, 'numeric')))
    if (isFALSE(is_all_numeric)) stop("All columns of covariates_data must be numeric")
    is_right_size <- nrow(covariates_data) == nrow(sequence_data)
    if (isFALSE(is_right_size)) stop("covariates_data must be same number of rows as sequence_data")
  }

  # data limit suggestion
  # TODO: is there a bette heuristic?
  if (nrow(sequence_data) * ncol(sequence_data) > 1e6) message("Consider taking a random sample of the data to increase performance")

  # pass the data to Shiny and launch the app
  shiny::shinyOptions(sequence_data = sequence_data,
                      covariates_data = covariates_data)
  shiny::runApp(system.file("sequenchr", package = "sequenchr"),
                launch.browser = TRUE)
}
