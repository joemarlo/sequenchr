
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sequenchr

<!-- badges: start -->
<!-- badges: end -->

Sequence analysis tool for applied researchers. Designed for faster
analysis iterations or for whom just prefer a point-and-click interface.

## Installation

You can install the latest version of sequenchr via:

``` r
# install.packages("devtools")
devtools::install_github("joemarlo/sequenchr")
```

## Example

``` r
library(TraMineR)
library(sequenchr)

# load data and convert to a sequence object
data(mvad)
seqstatl(mvad[, 17:86])
mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school",
                   "training")
mvad.labels <- c("employment", "further education", "higher education",
                 "joblessness", "school", "training")
mvad.scodes <- c("EM", "FE", "HE", "JL", "SC", "TR")
mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet, states = mvad.scodes,
                   labels = mvad.labels, xtstep = 6)

# launch the sequenchr app
launch_sequenchr(mvad.seq)


# Or use the plotting functions directly ....

# tidy the data and compute color palette
seq_def_tidy <- tidy_sequence_data(mvad.seq)
color_mapping <- viridis::viridis_pal()(length(alphabet(mvad.seq)))
names(color_mapping) <- alphabet(mvad.seq)

# set the ggplot theme
ggplot2::theme_set(ggplot2::theme_minimal())

# plot the sequence index
plot_sequence_index(seq_def_tidy, color_mapping)

# cluster the data
dist_matrix <- TraMineR::seqdist(seqdata = mvad.seq, method = "DHD")
cluster_model <- fastcluster::hclust(d = as.dist(dist_matrix), method = 'ward.D2')
cluster_assignments <- stats::cutree(cluster_model, k = 5)

# plot the sequence index by cluster
plot_sequence_index(seq_def_tidy, color_mapping, cluster_assignments = cluster_assignments)
```
