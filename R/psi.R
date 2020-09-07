#### psi_gaus ####

#' Title
#'
#' @param t
#'
#' @return
#' @export
#'
#' @examples
psi_gaus <- function(t){
  return(rep(.5, length(t)))
}

#### psi_huber ####

#' Title
#'
#' @param t
#' @param r
#' @param lst
#'
#' @return
#' @export
#'
#' @examples
psi_huber <- function(t, r, lst){
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

  psi <- numeric(length(t))
  psi[t <= cH^2] <- 1/(2*bH)
  psi[t > cH^2] <- cH^2 / (2*bH*t[t>cH^2])
  return(psi)
}

#### psi_t ####
#' Title
#'
#' @param t
#' @param r
#' @param nu
#'
#' @return
#' @export
#'
#' @examples
psi_t <- function(t, r, nu){
  return(.5 * (nu + r) / (nu + t))
}

#### psi_tukey ####

#' Computes psi(t) of the Tukey loss function
#'
#' @param t vector[N] squared Mahalanobis distances
#' @param cT scalar, tuning parameter
#'
#' @return psi psi of the Tukey loss
#' @export
#'
#' @examples
#' a <- c(-0.4121757 ,-2.2634588 , 0.2893038,  0.1831577,  0.4861016)
#' psi_tukey(a, .3)
psi_tukey <- function(t, cT){
  psi <- numeric(length(t))
  psi[t <= cT^2] <- 0.5 * (1 - t[t<=cT^2]/cT^2)^2
  return(psi)
}
