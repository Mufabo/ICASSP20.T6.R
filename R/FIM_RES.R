#### FIM_RES ####


#' computes FIM of one cluster for a given RES distribution
#'
#' @param x_hat Matrix[r, N_m] data matrix
#' @param t_m vector[N_m] squared Mahalanobis distances
#' @param S_est Matrix[r, r] Scatter matrix
#' @param psi scalar, psi of density generator
#' @param eta scalar, eta of density generator
#' @param D Matrix[r^2, 1/2*r*(r+1)] duplication matrix
#'
#' @return J Matrix[q, q] FIM
#' @export
#'
#' @examples
#'
FIM_RES <- function(x_hat, t_m, S_est, psi, eta, D){
  r <- dim(S_est)[1]
  N_m <- length(t_m)

  #F_mumu
  temp_eta <- tensorA::to.tensor(0, c(r, r, N_m))
  for(n in 1:N_m){
    temp_eta[,,n] <- eta(t_m[n]) * x_hat[, n] %o% x_hat[, n]
  }

  F_mumu <- -4 * S_est %\% colSums(temp_eta[,,3]) %/% S_est - S_est %\% diag(1, r, r) * sum(psi(t_m)) * 2

  #F_muS
  temp_eta <- tensorA::to.tensor(0, c(r, r^2, N_m))

  for(n in 1:N_m){
    temp_eta[,, n] <- eta(t_m[n]) * kronecker((S_est%\%x_hat[,n] %o% x_hat[,n]) %/% S_est, t(x_hat[,n]) %/% S_est)
  }
  F_muS <- -2 * colSums(temp_eta[,,3]) %*% D

  F_Smu <- t(F_muS)

  #F_SS
  temp_eta <- tensorA::to.tensor(0, c(r^2, r^2, N_m))
  for(n in 1:N_m){
    temp_eta[,,n] <- eta(t_m[n]) * kronecker((S_est %\% x_hat[,n] %o% x_hat[,n]) %/% S_est, (S_est %\% x_hat[,n] %o% x_hat[,n]) %/% S_est)
  }

  F_SS <- - t(D) %*% colSums(temp_eta[,,3]) %*% D - N_m / 2 %*% t(D) %*% (kronecker(S_est, S_est) %\% diag(1, r^2, r^2)) %*% D

  J <- cbind(rbind(-F_mumu, -F_muS), rbind(-F_Smu, -F_SS))

  # Ensure that J is PSD
  J <- if(det(J) < 0) nearestSPD(J) else J

  return(J)
}
