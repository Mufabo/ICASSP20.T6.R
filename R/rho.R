#### rho_gaus ####

#' rho(t) of the Gaussian distributiomn
#'
#' @param t vector[N] squared Mahalanobis distances
#' @param r scaler, dimension
#'
#' @return rho rho of the Gaussian
#' @export
#'
#' @note
#'
#' "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
#' Christian A. Schroth and Michael Muma, Signal Processing Group, Technische Universität Darmstadt
#' submitted to IEEE Transactions on Signal Processing
#'
#' @examples
#' rho_gaus(rnorm(5), 1)
rho_gaus <- function(t, r){
  return(r/2 * log(2*pi) + t/2)
}

#### rho_huber ####

#' Computes rho(t) of the Huber distribution
#'
#' @param t vector[N] squared Mahalanobis distances
#' @param r scalar, dimension
#' @param lst List of scalar tuning parameters. There are 3 valid options of this argument
#'        \enumerate{
#'        \item Empty list. Default case
#'        \item List containing qH
#'        \item List containing cH, bH, aH
#'        }
#'
#' @return vector[N] rho(t) of Huber distribution
#' @export
#'
#' @examples
#' rho_huber(rnorm(5), 2)
rho_huber <- function(t, r, lst=list()){
  if(length(lst)==0){
    qH <- 0.8
    cH <- sqrt(stats::qchisq(qH, r))
    bH <- stats::pchisq(cH^2, r+2) + cH^2 / r * (1 - stats::pchisq(cH^2, r))
    aH <- aH <- gamma(r/2)/pi^(r/2) / ( (2*bH)^(r/2) * (gamma(r/2) - pracma::incgam(r/2, cH^2/(2*bH))) + (2*bH*cH^r*exp(-cH^2/(2*bH)))/(cH^2-bH*r))
  } else {
    if(length(lst)==1){
      qH = lst[[1]]
      cH <- sqrt(stats::qchisq(qH, r))
      bH <- stats::pchisq(cH^2, r+2) + cH^2 / r * (1 - stats::pchisq(cH^2, r))
      aH <- gamma(r/2)/pi^(r/2) / ( (2*bH)^(r/2) * (gamma(r/2) - pracma::incgam(r/2, cH^2/(2*bH))) + (2*bH*cH^r*exp(-cH^2/(2*bH)))/(cH^2-bH*r))
    } else {
      if(length(lst)==3){
        cH <- lst[[1]]
        bH <- lst[[2]]
        aH <- lst[[3]]
      } else {
        stop("Argument lst invalid. It should either be empty, contain field qH or the fields cH and bH")
      }
    }
  }

  rho <- numeric(length(t))
  rho[t <= cH^2] <- -log(aH) + t[t<=cH^2]/(2*bH)
  rho[t >  cH^2] <- -log(aH) + cH^2/(2*bH) * log(t[t>cH^2]) - cH^2 / bH * log(cH) + cH^2 / (2 * bH)

  return(rho)
}


#### rho_t ####

#' Computes rho(t) of the t distribution
#'
#' @param t vector[N] squared Mahalanobis distance
#' @param r scalar, dimension
#' @param nu scalar, degree of freedom
#'
#' @return vector[N] rho(t) of the t distribution
#' @export
#'
#' @examples
#' rho_t(rnorm(5),2,5)
rho_t <- function(t, r, nu){
  return(-log(gamma((nu + r)/2)) + log(gamma(nu/2)) + 0.5*r*log(pi*nu) + (nu+r)/2 *log(1+t/nu))

}


#### rho_tukey ####

#' Computes rho(t) of the Tukey loss function
#'
#' @param t vector[N] squared Mahalanobis distance
#' @param r scalar, dimension
#' @param cT scalar, tuning parameter
#'
#' @return rho vector[N] rho of the Tukey loss function
#'
#' @note
#' "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
#' Christian A. Schroth and Michael Muma, Signal Processing Group, Technische UniversitÃ¤t Darmstadt
#' submitted to IEEE Transactions on Signal Processing
#' @export
#'
#' @examples
#'
#' rho_tukey(rnorm(5),2,5)
rho_tukey <- function(t, r, cT){
  rho <- numeric(length(t))
  rho[t<=cT^2] <- r/2 * log(2*pi) + t[t<=cT^2]^3 /(6*cT^4) - t[t<=cT^2]^2 /(2*cT^2) + t[t <= cT^2]/2
  return(rho)
  }
