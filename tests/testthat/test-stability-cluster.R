
clust_df <- stability_cluster(progeny_data, k = 3L,
                              n_iter = 250L, r_seed = 1)

test_that("`stability_cluster()` generates the correct object and dimensions", {
  expect_s3_class(clust_df, "tbl_df")
  expect_true(all(dplyr::select(clust_df, -prob_k) <= 1))
  expect_true(all(dplyr::select(clust_df, -prob_k) >= 0))
  expect_true(all(clust_df$prob_k %in% 1:3L))
})

test_that("`stability_cluster()` generates the expected output", {
  withr::local_options(list(pillar.sigfig = 6L))
  expect_snapshot( print(clust_df, n = Inf) )
})
