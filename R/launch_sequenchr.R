launch_sequenchr <- function(sequence_data, covariates_data = NULL){
  #https://stackoverflow.com/questions/49470474/saving-r-shiny-app-as-a-function-with-arguments-passed-to-the-shiny-app
  
  # quality checks
  is_TraMineR_object <- TraMineR::is.stslist(sequence_data)
  if (isFALSE(is_TraMineR_object)) stop("sequence_data must be a stslist object created from TraMineR::seqdef()")
  if (!is.null(covariates_data)){
    if (isFALSE(is(covariates_data, 'data.frame'))) stop("covariates_data must be a dataframe")
    is_all_numeric <- all(apply(covariates_data, 2, function(col) is(col, 'numeric')))
    if (isFALSE(is_all_numeric)) stop("All columns of covariates_data must be numeric")
    is_right_size <- nrow(covariates_data) == nrow(sequence_data)
    if (isFALSE(is_right_size)) stop("covariates_data must be same number of rows as sequence_data")
  }
  
  shiny::shinyOptions(sequence_data = sequence_data,
                      covariates_data = covariates_data)
  # source(system.file("sequenchr/app.R", local = TRUE, chdir = TRUE))$value
  source(file.path("sequenchr", "app.R"))$value
}
