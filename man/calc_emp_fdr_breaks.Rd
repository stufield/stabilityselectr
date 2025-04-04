% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc-emp-fdr-breaks.R
\name{calc_emp_fdr_breaks}
\alias{calc_emp_fdr_breaks}
\title{Calculate Empirical FDR Break Points}
\usage{
calc_emp_fdr_breaks(
  x,
  thresh_seq = seq(1, 0.1, by = -0.01),
  fdr_breaks = c(0.5, 1, 2, 3, 5)
)
}
\arguments{
\item{x}{A \code{stab_sel} class object generated
via \code{\link[=stability_selection]{stability_selection()}}.}

\item{thresh_seq}{\code{numeric(n)}. A sequence in
\verb{[0, 1]} specifying the thresholds to evaluate.}

\item{fdr_breaks}{\code{numeric(n)}. A vector specifying the desired mean number
of empirical false positives at which to calculate various thresholds.}
}
\value{
A list consisting of:
\item{n_selected}{A vector of the number of features selected at each
empirical stability selection threshold}
\item{meanFPs}{A vector of the mean number of false positive selected
features at each empirical stability selection threshold}
\item{breaks}{A \code{tibble} of containing empirical false positive
summary statistics at each FDR specified break point}
}
\description{
Calculates the stability selection threshold, the mean number of
false positive selected features (empirical), and the number of
selected features for specified FDR break points.
Relies on \code{\link[=calc_emp_fdr]{calc_emp_fdr()}} to calculate the mean false discovery
based on permutations during the stability selection algorithm.
}
\examples{
# l1-logistic
withr::with_seed(101, {
  n_feat      <- 20
  n_samples   <- 100
  x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
  colnames(x) <- paste0("feat", "_", head(letters, n_feat))
  y        <- sample(1:2, n_samples, replace = TRUE)
})
stab_sel <- stability_selection(x, y, "l1-logistic", num_iter = 25,
                                num_perms = 25,
                                r_seed = 101, parallel = TRUE)
calc_emp_fdr_breaks(stab_sel)
}
\seealso{
\code{\link[=get_stable_features]{get_stable_features()}}, \code{\link[=stability_selection]{stability_selection()}}

Other empirical FDR: 
\code{\link{calc_emp_fdr}()},
\code{\link{plot_emp_fdr}()}
}
\author{
Stu Field
}
\concept{empirical FDR}
