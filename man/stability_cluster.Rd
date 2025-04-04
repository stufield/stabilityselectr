% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability-cluster.R
\name{stability_cluster}
\alias{stability_cluster}
\title{Perform Stability Clustering}
\usage{
stability_cluster(data, k, iter = 100)
}
\arguments{
\item{data}{A (\verb{n x p}) data matrix containing \emph{n} samples and \emph{p} features.
Can also be a data frame where each row corresponds to a sample or
observation, whereas each column corresponds to a feature or variable.}

\item{k}{\code{integer(1)}. The number of clusters.}

\item{iter}{\code{integer(1)}. The number of random subset iterations to perform.}
}
\value{
A n x (k + 1) dimensional \code{tibble} of clustering probabilities for
each \code{k}, plus a final column named \code{"ProbK"}, which indicates
the "most probable" cluster membership for that sample.
}
\description{
Partitioning Around Medoids (PAM) is used both because
is uses a more robust measurement of the cluster centers (medoids) and
because this implementation keeps the cluster labels consistent across
runs, a key feature in calculating the across run stability. This does
not occur using \code{\link[=kmeans]{kmeans()}} where the initial cluster
labels are arbitrarily assigned.
}
\note{
How do we make sure clusters are indexed the same as what
comes out of k-means? Worried about index errors (but seems ok for now).
}
\examples{
stab_clust <- withr::with_seed(999, stability_cluster(progeny_data, k = 3, iter = 500))
table(actual = rep(1:3L, each = 50L), predicted = stab_clust$ProbK)

stab_clust$true_cluster <- rep(1:3L, each = 50L)

# View the stable clusters
stab_clust

# View the incorrectly clustered samples (n = 6)
filter(stab_clust, ProbK != true_cluster)

# Plot Stability Clusters
cols <- c("#24135F", "#00A499", "#840B55")
withr::with_par(list(mgp = c(2.00, 0.75, 0.0), mar = c(3, 4, 3, 1), mfrow = 1:2L), {
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
}
\references{
Hastie, et al. 2009.
}
\seealso{
\code{\link[=pam]{pam()}}

Other cluster: 
\code{\link{progeny_cluster}()}
}
\author{
Stu Field
}
\concept{cluster}
