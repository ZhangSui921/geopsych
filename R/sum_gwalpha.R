#' Counting the Number of GWalpha below Criterion
#'
#' an internal function for counting how many responses have an upper confidence
#' bound below a user-specified criterion for reliability level.
#'
#' @param bw bandwidth for weighting function
#' @param crit predetermined criterion for reliability level
#' @param ... additional arguments passed on to [gwalpha()]
#'
#' @return
#' the number of responses where the upper confidence bound of GWalpha is below
#' \code{crit}.
#'
sum_gwalpha <- function(bw, crit, ...) {
  sum(gwalpha(bw = bw, ci = TRUE, ...)$gwalpha_u < crit)
}
