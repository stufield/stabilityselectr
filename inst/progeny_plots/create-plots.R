# ---------------------------------------------------
# saving cluster data plots from progeny clustering
# from the Progeny Clustering paper data set:
# Hu, C.W., Kornblau, S.M., Slater, J.H. and A.A. Qutub (2015).
# Progeny Clustering: A Method to Identify Biological Phenotypes. Scientific Reports,
# 5:12894. http://www.nature.com/articles/srep12894
# ---------------------------------------------------
suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))

# Progeny Clustering ----
# This can take a while!
pclust <- withr::with_seed(101,
  progeny_cluster(progeny_data, clust_iter = 2:9L,
                  reps = 50L, iter = 100L, size = 8, verbose = TRUE)
)
plot(pclust, file = "inst/progeny_plots/progeny_data_output.pdf")

# Stability Clustering ----
stab_clust <- withr::with_seed(101, stability_cluster(progeny_data,
                                                      k = 3L, iter = 1000L))
stab_clust <- dplyr::mutate(
  stab_clust,
  sample = dplyr::row_number(),
  pch    = rep(16:18, each = 50),
  pch    = dplyr::case_when(
    sample <= 50 & ProbK != 1L ~ 13L,                # cluster 1
    sample > 50 & sample <= 100 & ProbK != 2L ~ 13L, # cluster 2
    sample > 100 & ProbK != 3L ~ 13L,                # cluster 3
    TRUE ~ pch
  )
)

cols <- unlist(col_palette[c("purple", "lightgreen", "magenta")])

withr::with_par(list(mgp = c(2.00, 0.75, 0.00), mar = c(3, 4, 3, 1), mfrow = 1:2L), {
  figure("inst/progeny_plots/progeny_data_stability_scatter.pdf",
         width = 12, height = 6)
  plot(progeny_data, col = rep(cols, each = 50), pch = rep(16:18, each = 50),
       cex = 1.75, main = "Simulated 3 Cluster Data")
  plot(progeny_data, col = rep(cols, each = 50), pch = stab_clust$pch, cex = 1.75,
       main = "Stability Clustering")
close_figure("inst/progeny_plots/progeny_data_stability_scatter.pdf")
})
