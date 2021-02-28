# set ggplot theme
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 16))

# retrieve pass-through variables
sequence_data <- shiny::getShinyOption("sequence_data")
covariates_data <- shiny::getShinyOption("covariates_data")


# pre-processing ----------------------------------------------------------

# establish color mapping
color_mapping <- viridis::viridis_pal()(length(alphabet(sequence_data)))
names(color_mapping) <- TraMineR::alphabet(sequence_data)

# tidy the data
tidy_data <- tidy_sequence_data(sequence_data)
if (!is.null(covariates_data)){
  tidy_cov_data <- covariates_data %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(sequenchr_seq_id = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = -sequenchr_seq_id)
} else tidy_cov_data <- NULL

# on exit, remove variables
# TODO this seems like a hack
onStop(function() {
  rm(list = c('sequence_data', 'tidy_data', 'color_mapping', 'tidy_cov_data'),
     envir = globalenv())
})
