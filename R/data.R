#' Package Objects
#'
#' @name objects
NULL

#' @describeIn objects The original clustering example data
#'   from the Progeny Clustering paper. A 2 feature data set
#'   example (see references) in the example illustration in
#'   Figure 1, page 3.
#' @format `clust_data`: A 2 column data matrix containing the 2 features in the
#'   example, named "F1" and "F2", and containing 20 observations.
#' @references `clust_data`: Hu, C.W., Kornblau, S.M., Slater, J.H. and
#'   A.A. Qutub (*2015*). Progeny Clustering: A Method to Identify
#'   Biological Phenotypes. Scientific Reports, 5:12894.
#'   pg. 3. \url{http://www.nature.com/articles/srep12894}
#' @source `clust_data`: Hu, et. al.
#' @examples
#' head(clust_data)
#'
"clust_data"

#' @describeIn objects A simulated clustering data set generated
#'   to contain **3 true clusters** in bivariate space. There are 2 features,
#'   named "F1" and "F2" and 150 observations.
#' @format `progeny_data`: A 2 column data matrix containing 2 features,
#'   "F1" and "F2", 150 observations, and 3 clusters in bivariate space.
#' @seealso \code{\link[MASS]{mvrnorm}}
#' @source `progeny_data`: Stu Field
#' @examples
#' plot(progeny_data, col = rep(2:4, each = 50),
#'      pch = rep(16:18, each = 50), cex = 1.75,
#'      main = "Simulated 3 Cluster Data")
"progeny_data"
