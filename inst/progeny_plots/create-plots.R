# ---------------------------------------------------
# saving cluster data plots from progeny clustering
# from the Progeny Clustering paper data set:
# Hu, C.W., Kornblau, S.M., Slater, J.H. and A.A. Qutub (2015).
#   Progeny Clustering: A Method to Identify Biological Phenotypes.
#   Scientific Reports, **5**:12894. http://www.nature.com/articles/srep12894
# ---------------------------------------------------
suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))


# Progeny Clustering ----

# This can take a while!
pclust <- withr::with_seed(101,
  progeny_cluster(progeny_data, clust_iter = 2:9L,
                  reps = 25L, iter = 100L, size = 8)
)

p <- plot(pclust)

ggplot2::ggsave(
  filename = "inst/progeny_plots/progeny_data_output.pdf",
  plot = p, height = 5, width = 10
)


# Stability Clustering ----
stab_clust <- withr::with_seed(101, stability_cluster(progeny_data,
                                                      k = 3L, iter = 1000L))
stab_clust$true_cluster <- rep(1:3L, each = 50L)

cols <- c("#24135F", "#00A499", "#840B55")
file <- "inst/progeny_plots/progeny_data_stability_scatter.pdf"

withr::with_pdf(file, title = basename(file), width = 12, height = 6, {
  withr::with_par(list(mgp = c(2.00, 0.75, 0.00), mar = c(3, 4, 3, 1),
                       mfrow = 1:2L), {
   plot(progeny_data,
        col = cols[stab_clust$true_cluster],
        bg  = cols[stab_clust$true_cluster],
        pch = stab_clust$true_cluster + 20,
        lwd = 1, cex = 1.75, main = "Simulated 3 Cluster Data")
   plot(progeny_data,
        col = cols[stab_clust$ProbK],
        bg  = cols[stab_clust$true_cluster],
        pch = stab_clust$true_cluster + 20,
        lwd = 2.5, cex = 1.5, main = "Stability Clustering")
  })
})
