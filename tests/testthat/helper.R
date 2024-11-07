
n_feat    <- 20
n_samples <- 100
withr::with_seed(101, {
  x <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
  y <- sample(1:2, n_samples, replace = TRUE)
})

colnames(x) <- paste0("feat", "_", head(letters, n_feat))

withr::with_options(list(signal.quiet = TRUE), {
  ss <- stability_selection(x, y, kernel = "l1-logistic", r_seed = 101)
  ss_perm10 <- stability_selection(x, y, kernel = "l1-logistic",
                                  num_perms = 10, r_seed = 101)
})
