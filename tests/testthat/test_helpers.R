library(TraMineR)
data(mvad)
seqstatl(mvad[, 17:86])
mvad.alphabet <- c("employment", "FE", "HE", "joblessness", "school",
                   "training")
mvad.labels <- c("employment", "further education", "higher education",
                 "joblessness", "school", "training")
mvad.seq <- seqdef(mvad, 17:86, alphabet = mvad.alphabet, # states = mvad.scodes,
                   labels = mvad.labels, xtstep = 6)
dist_matrix <- TraMineR::seqdist(seqdata = mvad.seq, method = "DHD")
cluster_model <- hclust(d = as.dist(dist_matrix), method = 'ward.D2')

separation_metrics <- cluster_stats(
 dist_matrix = as.dist(dist_matrix),
 cluster_model = cluster_model,
 k_min = 2,
 k_max = 5
)
test_that("cluster_stats() output is rectangular", {
  expect_true(inherits(separation_metrics, 'data.frame'))
  expect_equal(nrow(separation_metrics), 4)
  expect_equal(ncol(separation_metrics), 5)
  expect_named(separation_metrics, c('k', 'ch', 'silhouette', 'ch_norm', 'silhouette_norm'))
})

seq_def_tidy <- tidy_sequence_data(mvad.seq)
test_that("tidy_sequence_data() output is rectangular", {
  expect_true(inherits(seq_def_tidy, 'data.frame'))
  expect_equal(nrow(seq_def_tidy), 49840)
  expect_equal(ncol(seq_def_tidy), 3)
  expect_named(seq_def_tidy, c('sequenchr_seq_id', 'period', 'state'))
})

entropy <- shannon_entropy(c('A', 'A', 'B', 'C', 'A', 'C'))
test_that("shannon_entropy() output is correct", {
  expect_true(inherits(entropy, "numeric"))
  expect_equal(entropy, 1.459148, tolerance = 1e-5)
})

cluster_labels <- label_clusters(cluster_model, k = 5)
test_that("label_clusters() output is correct", {
  expect_true(inherits(cluster_labels, 'factor'))
  expect_equal(length(cluster_labels), 712)
})

trans_tidy <- transition_matrix(seq_def_tidy, cluster_labels = cluster_labels)
test_that("transition_matrix() output is rectangular", {
  expect_true(inherits(trans_tidy, 'data.frame'))
  expect_equal(nrow(trans_tidy), 180)
  expect_equal(ncol(trans_tidy), 4)
  expect_named(trans_tidy, c('current', 'cluster', 'previous', 'n'))
})
trans_tidy_sums <- trans_tidy %>% dplyr::group_by(cluster) %>% dplyr::summarize(sum = sum(n)) %>% dplyr::pull(sum)
test_that("transition_matrix() output sums to one", {
  expect_equal(trans_tidy_sums, rep(1, 5))
})
