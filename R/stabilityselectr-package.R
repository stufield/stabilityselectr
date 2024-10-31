
#' @keywords internal package
"_PACKAGE"

#' @name stabilityselectr-package
#'
#' @details Stability selection is performed using [stabilitySelection()],
#'   which returns a `stab_sel` class object. The stability path can be
#'   plotted using [graphics::plot()]. A `tibble` with the highest
#'   selection probability can be created using [getStableFeatures()].
#'   See vignettes.
#'
#' @seealso [glmnet::glmnet()], [kmeans()]
#' @import helpr splyr
#' @references **Meinshausen, N. and Buhlmann, P.** (*2010*), Stability selection.
#'   Journal of the Royal Statistical Society: Series B (Statistical
#'   Methodology), 72: 417-473. doi: 10.1111/j.1467-9868.2010.00740.x
#' @references **Hu, C.W., Kornblau, S.M., Slater, J.H. and A.A. Qutub** (*2015*).
#'   Progeny Clustering: A Method to Identify Biological Phenotypes. Scientific Reports,
#'   5:12894. \url{http://www.nature.com/articles/srep12894}
NULL
