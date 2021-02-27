#' Launches the sequenchr Shiny app
#'
#' Launches the sequenchr app in the default web browser. Optional covariates_data produces plots of the marginal distributions.
#'
#' @param sequence_data an object created from TraMineR::seqdef
#' @param covariates_data optional. a dataframe with the same number of rows as sequence_data
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' launch_sequenchr()
#' }
launch_sequenchr <- function(MonteCarlo = 10){ #sequence_data, covariates_data = NULL){

  # quality checks
  # is_TraMineR_object <- TraMineR::is.stslist(sequence_data)
  # if (isFALSE(is_TraMineR_object)) stop("sequence_data must be a stslist object created from TraMineR::seqdef()")
  # if (!is.null(covariates_data)){
  #   if (isFALSE(methods::is(covariates_data, 'data.frame'))) stop("covariates_data must be a dataframe")
  #   is_all_numeric <- all(apply(covariates_data, 2, function(col) methods::is(col, 'numeric')))
  #   if (isFALSE(is_all_numeric)) stop("All columns of covariates_data must be numeric")
  #   is_right_size <- nrow(covariates_data) == nrow(sequence_data)
  #   if (isFALSE(is_right_size)) stop("covariates_data must be same number of rows as sequence_data")
  # }

  # getOption("shiny.launch.browser", TRUE)

  shiny::runApp(system.file("sequenchr", package = "sequenchr"),
                launch.browser = TRUE)

  # shiny::shinyOptions(sequence_data = sequence_data,
  #                     covariates_data = covariates_data)
  # source(file.path("sequenchr", "app.R"))$value
  # shiny::shinyOptions(MonteCarlo = MonteCarlo)
  # source(system.file("app.R", package = "sequenchr", local = TRUE, chdir = TRUE))$value
}
