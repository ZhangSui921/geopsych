#' Bandwidth Selection for Geographically Weighted Cronbach's Alpha
#'
#' A function for automatic bandwidth selection to calibrate GWalpha.
#'
#' @param crit predetermined criterion for reliability level
#' @param minmax a numeric vector of length 2 with the lower and upper bounds
#'   of the bandwidth search interval
#' @param adaptive logical; if \code{TRUE} (default), the bandwidth \code{bw}
#'   corresponds to the number of nearest neighbours; if \code{FALSE},
#'   \code{bw} is a fixed distance
#' @param tol convergence tolerance
#' @param max_iter maximum number of iterations
#' @param ... additional arguments passed on to [gwalpha()].
#'
#' @return
#' a adaptive or fixed distance bandwidth
#'
#' @seealso [gwalpha()]
#'
#' @export
bw_gwalpha <- function(crit, minmax, adaptive = TRUE, tol = 3, max_iter = 100, ...) {
  if (adaptive) tol <- 3 else tol <- 1e-03

  phi <- (sqrt(5) - 1) / 2
  a <- minmax[1]
  b <- minmax[2]
  x1 <-  b - phi * (b - a)
  x2 <- a + phi * (b - a)
  f1 <- sum_gwalpha(x1, crit, ...)
  f2 <- sum_gwalpha(x2, crit, ...)

  iter <- 0
  while ((b - a) > tol && iter < max_iter) {
    if (f1 > f2) {
      b  <- x2
      x2 <- x1
      f2 <- f1
      x1 <- b - phi * (b - a)
      f1 <- sum_gwalpha(x1, crit, ...)
    } else {
      a  <- x1
      x1 <- x2
      f1 <- f2
      x2 <- a + phi * (b - a)
      f2 <- sum_gwalpha(x2, crit, ...)
    }
    iter <- iter + 1
  }
  if (f1 > f2) opt <- x1 else opt <- x2
  if (adaptive) opt <- trunc(opt)

  return(opt)
}
