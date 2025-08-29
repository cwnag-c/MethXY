#' m2b
#' @param m Mval
#'
#' @return beta matrix
#' @export
#'
m2b <- function(m) {
  b <- 2^m / (2^m + 1)
  return(b)
}


#' b2m
#'
#' @param b beta matrix
#'
#' @return M matrix
#' @export
#'
b2m <- function(b) {
  b1 <- pmin(pmax(b, 1e-6), 1 - 1e-6)
  m <- log2(b1 / (1 - b1))
  return(m)
}
