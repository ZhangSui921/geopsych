#' Geographically Weighted Cronbach's Alpha
#'
#' This function computes geographically weighted cronbach's alpha (GWalpha)
#' for responses from georeferenced multi-item survey scales.
#'
#' @param x a character vector or numeric vector with the column names or indices
#'   of items in \code{data}
#' @param data a SpatialPointsDataFrame as defined in package \code{sp}, or a sf object
#'   defined in package \code{sf}
#' @param kernel type of kernel function used to weight responses. Available
#'   options: \code{"bisquare"} (default), \code{"gaussian"}, \code{"exponential"},
#'   or \code{"boxcar"}
#' @param adaptive logical; if \code{TRUE} (default), the bandwidth \code{bw}
#'   corresponds to the number of nearest neighbours; if \code{FALSE},
#'   \code{bw} is a fixed distance
#' @param bw bandwidth for weighting function, can be specified or obtained using
#'   [bw_gwalpha()]
#' @param ci logical; if \code{TRUE}, bootstrapped confidence intervals are
#'   computed
#' @param p the percentile for the upper confidence interval if \code{ci = TRUE}
#' @param nsims number of bootstrap iterations if \code{ci = TRUE}
#'
#' @return
#' a data frame with:
#' \item{\code{gwalpha}}{local estimates of reliability}
#' \item{\code{gwalpha_u}}{upper confidence interval if \code{ci = TRUE}}
#' \item{\code{coords}}{coordinates matrix for each responses}
#'
#' @references Zhang, S., and Z. Li. 2025. “ Geographically Weighted Cronbach's
#'  Alpha (GWalpha): An Exploratory Local Measure of Reliability for Scale Construction.”
#'   Geographical Analysis 57, no. 4: 758–772.
#'
#' @seealso [bw_gwalpha()]
#'
#' @examples
#' data(bes2011)
#' alpha100 <- gwalpha(
#'   x = c("willHelp", "closeKnit", "trust", "solComProb", "relGroups"),
#'   data = bes2011,
#'   bw = 100
#'   )
#'
#' @importFrom stats cov.wt dist quantile
#'
#' @export
gwalpha <- function(x, data, kernel = 'bisquare', adaptive = TRUE, bw,
                    ci = FALSE, p = .95, nsims = 1000) {
  coords <- sp::coordinates(data)
  data <- as.data.frame(data)
  n <- nrow(coords)
  x_v <- as.matrix(data[, x])
  n_v <- ncol(x_v)
  mv <- numeric(n)
  mcov <- numeric(n)
  gwalpha <- numeric(n)
  if (ci) {
    b_mv <- matrix(nrow = n, ncol = nsims)
    b_mcov <- matrix(nrow = n, ncol = nsims)
    b_gwalpha <- matrix(nrow = n, ncol = nsims)
  }

  d <- as.matrix(dist(coords))
  for (i in 1:n) {
    d_i <- d[, i]
    if (adaptive) {
      h <- sort(d_i)[bw]
    } else {
      h <- bw
    }
    if (kernel == 'bisquare') {
      w_i <- ifelse(d_i > h, 0, (1 - (d_i/h)^2)^2)
    } else if (kernel == 'gaussian') {
      w_i <- exp(-0.5 * (d_i/h)^2)
    } else if (kernel == 'exponential') {
      w_i <- exp(-d_i/h)
    } else {
      w_i <- ifelse(d_i > h, 0, 1)
    }

    wi <- w_i / sum(w_i)
    wcov <- cov.wt(x_v, wt = wi)$cov
    mv[i] <- sum(diag(wcov)) / n_v
    mcov[i] <- (sum(wcov) - sum(diag(wcov))) / (n_v * (n_v - 1))
    gwalpha[i] <- n_v * mcov[i] / (mv[i] + (n_v - 1) * mcov[i])

    if (ci) {
      b_wi <- wi[order(d_i)[1:bw]]
      for (j in 1:nsims) {
        b <- sample(order(d_i)[1:bw], replace = T)
        b_v <- x_v[b, , drop = FALSE]
        b_wcov <- cov.wt(b_v, wt = b_wi)$cov
        b_mv[i, j] <- sum(diag(b_wcov)) / n_v
        b_mcov[i, j] <- (sum(b_wcov) - sum(diag(b_wcov))) / (n_v * (n_v - 1))
        b_gwalpha[i, j] <- n_v * b_mcov[i, j] / (b_mv[i, j] + (n_v - 1) * b_mcov[i, j])
      }
    }
  }

  if (ci) {
    gwalpha_u <- apply(b_gwalpha, 1, function(x) quantile(x, p))
    res <- data.frame(gwalpha, gwalpha_u, coords)
  } else {
    res <- data.frame(gwalpha, coords)
  }

  return(res)
}
