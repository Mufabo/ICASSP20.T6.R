#' Computes eta(t) of the Gaussian distribution
#'
#' @param t vector of length N. Squared Mahalanobis distances
#' @return vector of length N. eta of a Gaussian distribution
#'
#' @examples
#'
#' eta_gaus(c(1,2,3))
#'
#' @export
eta_gaus <- function(t) return(numeric(length(t)))

#' Computes eta(t) of the Huber distribution
#'
#' @param t Vector. Squared Mahalanobis distances
#' @param r Scalar. Dimension
#' @param lst List. Either empty, Contains field qH or the fields cH and bH. Default is set to empty
#'
#' @return eta of the Huber distribution
#'
#' @examples
#' eta_huber(rnorm(4), 2, list(.1))
#' eta_huber(rnorm(4), 2)
#'
#' @export
eta_huber <- function(t, r, lst=list()){
  if(length(lst)==0){
    qH <- 0.8
    cH <- sqrt(stats::qchisq(qH, r))
    bH <- stats::pchisq(cH^2, r+2) + cH^2 / r * (1 - stats::pchisq(cH^2, r))
  } else {
    if(length(lst)==1){
      qH = lst[[1]]
      cH <- sqrt(stats::qchisq(qH, r))
      bH <- stats::pchisq(cH^2, r+2) + cH^2 / r * (1 - stats::pchisq(cH^2, r))
    } else {
      if(length(lst)==2){
        cH = lst[[1]]
        bH = lst[[2]]
      } else {
        stop("Argument lst invalid. It should either be empty, contain field qH or the fields cH and bH")
      }
    }

  }

  eta <- numeric(length(t))
  if(all((2 * bH * t(t > cH^2)^2) == 0)){
    eta[t > cH^2] <- 0
  } else {
    eta[t > cH^2] <- -cH^2 / (2 * bH * t[t > cH^2]^2)}

  return(eta)
}

#'Computes eta(t) of the t distribution
#'
#' @param t Vector of length N. Squared Mahalanobis distances
#' @param r Scalar. Dimension
#' @param nu Scalar. Degree of Freedom
#'
#' @return eta(t) of the t distribution
#'
#' @examples
#' eta_t(c(-0.29022179, -0.58196155, -0.51125062, -0.09460791), .1, .3)
#'
#'
#' @export
eta_t <- function(t, r, nu){
  return(-0.5 * (nu + r) / (nu + t)^2)
}

#' Computes eta(t) of the Tukey loss function
#'
#' @param t Vector. Squared Mahalanobis distances
#' @param cT Scalar tuning parameter
#'
#' @return eta of Tukey
#'
#' @examples
#'
#' eta_tukey(c(-1.4456554,  0.9074437,  0.2274346, -0.8734823), .1)
#'
#' @export
eta_tukey <- function(t, cT){
  eta <- numeric(length(t))
  eta[t <= cT^2] <- t[t <= cT^2] / cT^4 - 1 / cT^2
  return(eta)
}
