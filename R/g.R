#### g_gaus ####

#' Title
#'
#' @param t
#' @param r
#' @examples
#'
#' @return
#' @export
g_gaus <- function(t, r){
  clip <- 400 # clipping to avoid 0
  g <- (2 * pi)^(-r/2) * exp(-t/2)
  g[clip <= t] <- (2 * pi)^(-r/2) * exp(-clip / 2)
  return(g)
}

#### g_huber ####

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
g_huber <- function(t, r, lst){
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
  g <- numeric(length(t))
  g[t <= cH^2] <- aH * exp(-t[t<=cH^2]/(2*bH))
  g[t >  cH^2] <- aH * (exp(1) * t[t>cH^2]/cH^2)^(-cH^2/(2*bH))
  return(g)
}

#### g_t ####

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
g_t <- function(t, r, nu){
  return(gamma((nu + r)/2) / (gamma(nu/2)*(pi*nu)^(r/2)) * (1 + t/nu)^(-(nu+r)/2))

}
