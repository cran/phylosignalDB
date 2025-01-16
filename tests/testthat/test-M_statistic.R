test_that("M_stat() works", {
  trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
  trait_dist <- gower_dist(x = trait_df)
  res <- M_stat(trait_dist, turtles$phylo)
  expect_equal(round(res, 7), c(M_stat = 0.6589147))
})

test_that("M_rand_perm() works", {
  trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
  trait_dist <- gower_dist(x = trait_df)
  set.seed(1314)
  res <- M_rand_perm(trait_dist, turtles$phylo, reps = 4)
  expect_equal(round(res$M_observed, 7), 0.6589147)
  expect_equal(round(res$M_permuted, 7),
               c(0.5348837, 0.6201550, 0.5930233, 0.5930233))
})

test_that("phylosignal_M() works", {
  trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
  trait_dist <- gower_dist(x = trait_df)
  set.seed(1314)
  res <- phylosignal_M(trait_dist, turtles$phylo, reps = 99)
  expect_equal(round(res$stat, 7), 0.6589147)
  expect_true(res$pvalue < 0.05)
})


test_that("Parallel computing works", {
  skip_on_cran()
  trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
  trait_dist <- gower_dist(x = trait_df)
  set.seed(1314)
  ptime_core1 <- system.time({res_core1 <- phylosignal_M(trait_dist, turtles$phylo, reps = 999, cores = 1)})
  expect_true(ptime_core1["elapsed"] < 15)
  #Parallel computing shows significant speed benefits when a phylogenetic tree has over 200 tips and more than 999 random permutations.
  # set.seed(1314)
  # ptime_parallel <- system.time({res_parallel <- phylosignal_M(trait_dist, turtles$phylo, reps = 999, cores = 0)})
  # expect_equal(res_core1$stat, res_parallel$stat)
  # expect_equal(res_core1$pvalue < 0.05, res_parallel$pvalue < 0.05)
  # expect_equal(ptime_core1["elapsed"] > ptime_parallel["elapsed"], c(elapsed = TRUE))
})

