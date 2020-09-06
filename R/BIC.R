#### BIC A ####

#' Computes the BIC of a RES distribution with the asymptotic penalty term
#'
#' @param S_est tensorA::tensor[r, r, ll]. Estimated scatter matrix of cluster ll
#' @param t Matrix[N, ll] Squared Mahalanobis distances of data points to cluster centers
#' @param mem Matrix[N, ll] Cluster memberships
#' @param rho Vector rho of density generator
#' @param psi Vector. psi of density generator
#' @param eta Vector. eta of density generator
#'
#' @return list
#' \enumerate{
#' \item bic
#' \item pen : Penalty term
#' \item like : Likelihood term
#'}
#'
#' @examples
#'
#'
#' @export
BIC_A <- function(S_est, t, mem, rho, psi, eta){
  N_m <- colSums(mem)
  r <- dim(S_est)[1]
  ll <- dim(S_est)[3]
  q <- 0.5 * r * (r + 3)

  temp_rho <- numeric(ll)
  temp_psi <- numeric(ll)
  temp_eta <- numeric(ll)
  logdetS <- numeric(ll)
  epsilon <- numeric(ll)

  for (m in 1:ll) {
    temp_rho[m] <- colSums(rho(t[mem[, m], m]))
    temp_psi[m] <- colSums(psi(t[mem[, m], m]))
    temp_eta[m] <- colSums(eta(t[mem[, m], m]))

    logdetS[m] <- log(det(S[,,m]))
    epsilon[m] <- max(abs(temp_psi[m]), abs(temp_eta[m]), N_m[m])
  }

  like <- - sum(temp_rho[temp_rho > 0]) + sum(N_m[N_m > 0] * log(N_m[N_m > 0])) - sum(N_m * logdetS)/2
  pen <- -0.5 * q * sum(epsilon[epsilon > 0])

  bic <- like + pen

  return(list('bic' = bic, 'pen' = pen, 'like' = like))
}

#### BIC F ####

#' computes the BIC of a RES distribution based on a finite sample penalty term
#'
#' @param data Matrix[N, r] data samples
#' @param S_est tensorA::Tensor[r, r, ll] Scatter matrices of the clusters
#' @param mu_est Matrix[r, ll] estimated mean vectors of all clusters
#' @param t Matrix[N, ll] Squared Mahalanobis distances to cluster centers
#' @param mem Matrix[N, ll] cluster memberships
#' @param rho Vector rho of density generator
#' @param psi Vector. psi of density generator
#' @param eta Vector. eta of density generator
#'
#' @return list
#' \enumerate{
#' \item bic
#' \item pen : Penalty term
#' \item like : Likelihood term
#'}
#' @note
#'
#'
#' "Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
#' Christian A. Schroth and Michael Muma, Signal Processing Group, Technische UniversitÃ¤t Darmstadt
#' submitted to IEEE Transactions on Signal Processing
#'
#' @export
BIC_F <- function(data, S_est, mu_est, t, mem, rho, psi, eta){
  N_m <- colSums(mem)

  r <- tensorA::dim(S_est)[1]
  ll<- tensorA::dim(S_est)[3]
  D <- ICASSP20.T6.R::duplicationMatrix(r)
  q <- 1/2*r*(r+3)

  temp_rho <- numeric(ll)
  logdetS <- numeric(ll)
  detJ <- numeric(ll)

  for (m in 1:ll) {
    x_hat_m <- t(data[mem[,m], 2:dim(data)[2]]) - mu_est[, m]
    t_m <- t[mem[,m], m]
    J <- FIM_RES(x_hat_m, t_m, S_est[,,m], psi, eta, D);
    detJ[m] <- det(J)
    temp_rho[m] = sum(rho(t(mem[,m], m)))
    logdetS[m] = log(det(S_est[,,m]))

    if(detJ[m] < 0){
      warning("negative determinant, J still not positive definite")
      detJ[m] <- detJ[m] + 10^-10
      if(detJ[m] < 0) detJ[m] <- abs(detJ[m])
    } else if(detJ[m] == 0 && N_m[m] == 0){
      warning("cluster without data point, zero determinant")
      detJ[m] <- 1
    } else if(detJ[m] == 0){
      warning("zero determinant")
      detJ[m] <- detJ[m] + 10^-10
      if(detJ[m] < 0) detJ[m] <- abs(detJ[m])
    }
  }

  like <- -sum(temp_rho) + sum(N_m[N_m > 0] * log(N_m[N_m > 0])) - sum(N_m[N_m > 0] * logdetS[N_m > 0]) / 2
  pen <- - 1/2 * sum(log(detJ)) + ll * q/2 * log(2*pi) - ll * log(ll)

  bic <- like + pen

  return(list(bic=bic, like=like, pen=pen))
}

#### BIC S ####

#' Title
#'
#' @param S_est
#' @param t
#' @param mem
#' @param rho
#'
#' @return
#' @export
#'
#' @examples
BIC_S <- function(S_est, t, mem, rho){
  N_m <- colSums(mem)
  r <- dim(S_est)[1]
  ll <- dim(S_est)[3]

  q <- .5 * r * (r+3)

  N <- dim(t)[1]
  N <- if(N == 0) 1 else N

  temp_rho <- numeric(ll)
  logdetS <- numeric(ll)

  for(m in 1:ll){
    temp_rho[m] <- sum(rho(t[mem[,m], m]))
    logdetS[m] <- log(det(S_est[,,m]))
  }

  like <- -sum(temp_rho) + sum(N_m[N_m > 0]) * log(N_m[N_m > 0]) - sum(N_m[N_m > 0] * logdetS[N_m > 0])/2
  pen <- - q * ll/2 * log(N)
  bic <- like + pen

  return(list(bic=bic, pen=pen, like=like))
  }


