
create_local_objects <- function(n = 100, p = 20, env = parent.frame()) {
  withr::local_options(list(signal.quiet = TRUE))
  withr::with_seed(101, {
    nms <- paste0("feat", "_", head(letters, p))
    x <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, nms))
    y <- sample(1:2L, n, replace = TRUE)
  })

  ss <- stability_selection(x, y, kernel = "l1-logistic", r_seed = 101)
  ss_perm10 <- stability_selection(x, y, kernel = "l1-logistic",
                                  num_perms = 10, r_seed = 101)
  assign("x", x, envir = env)
  assign("y", y, envir = env)
  assign("n_feat", p, envir = env)
  assign("n_samples", n, envir = env)
  assign("ss", ss, envir = env)
  assign("ss_perm10", ss_perm10, envir = env)
}

create_local_objects()
