
create_local_objects <- function(n = 100, p = 20, env = teardown_env()) {
  withr::local_options(list(signal.quiet = TRUE))
  nms       <- paste0("feat", "_", head(letters, p))
  n_feat    <- p
  n_samples <- n
  withr::with_seed(101, {
    x <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, nms))
    y <- sample(1:2L, n, replace = TRUE)
  })
  ss <- stability_selection(x, y, alpha = 0.8, Pw = 0.5, r_seed = 101)
  ss_perm <- stability_selection(x, y, alpha = 0.8, Pw = 0.5,
                                 n_perm = 10, r_seed = 101)
  get_env <- function() parent.frame()
  get("attach")(get_env(), name = "ss_testing", pos = 2L,
                warn.conflicts = FALSE)
  withr::defer(
    get("detach")("ss_testing", character.only = TRUE, force = TRUE),
    envir = env
  )
}

create_local_objects()
