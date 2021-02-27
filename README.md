
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sequenchr

<!-- badges: start -->
<!-- badges: end -->

Sequence analysis tool for applied researchers. Designed for faster
analysis iterations or for whom just prefer a point-and-click interface.

## Installation

You can install the lastest version sequenchr via:

``` r
# install.packages("devtools")
devtools::install_github("joemarlo/sequenchr")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(sequenchr)
library(TraMineR)

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
```

Or use the plotting functions directly ….

``` r
# TODO
```
