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
  #browser()
  r <- dim(S_est)[1]
  N_m <- length(t_m)

  #F_mumu
  temp_eta <- tensorA::to.tensor(0, c(r, r, N_m))
  for(n in 1:N_m){
    # For whatever reason eta returns a matrix
    tmp <- if(is.matrix(eta(t_m[n]))) eta(t_m[n])[,]  else eta(t_m[n])
    temp_eta[,,n] <- tmp * x_hat[, n] %*% t(x_hat[, n])
  }

  F_mumu <- -4 * S_est %\% matSum(temp_eta) %/% S_est - S_est %\% diag(1, r, r) * sum(psi(t_m)) * 2

  #F_muS
  temp_eta <- tensorA::to.tensor(0, c(r, r^2, N_m))

  for(n in 1:N_m){
    tmp <- if(is.matrix(eta(t_m[n]))) eta(t_m[n])[,]  else eta(t_m[n])
    temp_eta[,, n] <- tmp * kronecker((S_est%\%x_hat[,n] %*% t(x_hat[,n])) %/% S_est, t(x_hat[,n]) %/% S_est)
  }

  F_muS <- -2 * matSum(temp_eta) %*% D

  F_Smu <- t(F_muS)

  #F_SS
  temp_eta <- tensorA::to.tensor(0, c(r^2, r^2, N_m))
  for(n in 1:N_m){
    tmp <- if(is.matrix(eta(t_m[n]))) eta(t_m[n])[,]  else eta(t_m[n])
    temp_eta[,,n] <- tmp * kronecker((S_est %\% x_hat[,n] %*% t(x_hat[,n])) %/% S_est, (S_est %\% x_hat[,n] %*% t(x_hat[,n])) %/% S_est)
  }

  F_SS <- - t(D) %*% matSum(temp_eta) %*% D - N_m / 2 * t(D) %*% (kronecker(S_est, S_est) %\% diag(1, r^2, r^2)) %*% D

  J <- rbind(cbind(-F_mumu, -F_muS), cbind(-F_Smu, -F_SS))

  # Ensure that J is PSD
  J <- if(det(J) < 0) nearestSPD(J) else J

  return(J)
}
