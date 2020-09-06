#' Title
#'
#' @param x
#' @param mu
#' @param S
#'
#' @return
#' @export
#'
#' @note
#'
#'
#'
#' @examples
mahalanobisDistance <- function(x, mu, S){
  # Ensure that S is not ill conditioned beforehand
  t <- pracma::dot(((x - mu) %/% S), (x - mu))
  return(t)
}
