#' Title
#'
#' @param data Matrix[N, r] data without labels
#' @param ll scalar, number of clusters
#' @param g function(t), gauss
#' @param psi
#' @param limit
#' @param em_max_iter
#' @param reg_value
#'
#' @return
#' @export
#'
#' @examples
EM_RES <- function(data, ll, g, psi, limit = 1e-6, em_max_iter = 200, reg_value = 1e-6, test_args = NULL){
  # Variable initializations
  r <- dim(data)[2]
  N <- dim(data)[1]

  v <- matrix(0, N, ll)
  v_diff <- matrix(0, N, ll)
  tau <- numeric(ll)
  S_hat <- tensorA::to.tensor(0,c(r, r, ll))
  t <- matrix(0, N, ll)
  log_likelihood <- numeric(em_max_iter)

  ## Initialization using K-means++
  if (is.null(test_args)){
    tmp <- ClusterR::KMeans_rcpp(data, ll, num_init = 10)
    clu_memb_kmeans <- tmp$clusters
    mu_Kmeans <- tmp$centroids
    mu_hat <- t(mu_Kmeans) #stores centroids in columns
  } else{
    mu_hat <- test_args$mu.hat
    clu_memb_kmeans <- test_args$clu.memb.kmeans
  }

  for(m in 1:ll){
    tau[m] <- sum(clu_memb_kmeans == m)/N
    }

  for(m in 1:ll){
    x_hat <- matrix(data = data[clu_memb_kmeans == m,][,1], nrow = length(data[clu_memb_kmeans == m,][1,]), ncol = length(data[clu_memb_kmeans == m,][,1]), byrow = TRUE) - mu_hat[, m]
    N_m <- sum(clu_memb_kmeans == m)
    S_hat[,,m] <- (x_hat %*% t(x_hat))/N_m

    # Check if the sample covariance matrix is positive definite
    if( all( eigen(S_hat[,,m], only.values = TRUE)$values < 0) || pracma::cond(S_hat[,,m]) > 30){
      S_hat[,,m] <- 1 / (r*N_m) * sum(diag(x_hat %*% t(x_hat)))*diag(1,r,r)
      if(all( eigen(S_hat[,,m], only.values = TRUE)$values < 0)){
        S_hat[,,m] <- diag(1, r, r)
      }
    }
    t[, m] <-stats::mahalanobis(data, mu_hat[,m], tensorA::to.matrix.tensor(S_hat[,,m],1))
  }

  # EM Algorithm
  for(ii in 1:em_max_iter){
    # E-step
    v_lower <- matrix(0, N, ll)
    for(j in 1:ll){
      v_lower[,j] <- tau[j] * det(S_hat[,,j])^(-.5) * g(t[,j])
    }

    for(m in 1:ll){
      v[,m] <- tau[m] * det(S_hat[,,m])^-.5 * g(t[,m]) / rowSums(v_lower)
      v_diff[,m] <- v[,m] * psi(t[,m])
    }

    # M-step
    for(m in 1:ll){

      mu_hat[,m] <- colSums(v_diff[,m] * data) / sum(v_diff[,m])
      S_hat[,,m] <- 2 * (matrix(v_diff[,m], nrow=ncol(data), ncol=nrow(data),byrow = TRUE) * sweep(t(data), 1, mu_hat[,m])) %*% t(sweep(t(data), 1,mu_hat[,m])) / sum(v[,m]) + reg_value * diag(1, r, r)
      tau[m] <- sum(v[,m])/N
      t[,m] <- stats::mahalanobis(data, mu_hat[,m], tensorA::to.matrix.tensor(S_hat[,,m],1))
    }

    # Convergence check
    v_conv <- matrix(0, N, ll)
    for(m in 1:ll){
      v_conv[,m] <- tau[m] * det(S_hat[,,m])^-.5 * g(t[,m])
    }

    log_likelihood[ii] <- sum(log(rowSums(v_conv)))

    if(ii > 1){
      if(abs(log_likelihood[ii] - log_likelihood[ii-1]) < limit){
        break
      }
    }
  }

  # Calculate posterior probabilities
  R <- v_conv / rowSums(v_conv)

  # diagonal loading
  # If the estimated "Matrix is close to singular or badly scaled" it cannot be inverted.
  # https://math.stackexchange.com/questions/261295/to-invert-a-matrix-condition-number-should-be-less-than-what
  # If S_hat has a large condition number, a small number in comparison to the
  # matrix entries is added. This step should be subject for further tweaking.
  for(m in 1:ll){
    cond_S <- pracma::cond(S_hat[,,m])
    if(cond_S > 30){
      warning("S with large condition number")
      S_hat[,,m] <- S_hat[,,mn] + 0.01 * 0^floor(log10(trace(S_hat[,,m]))) * log10(cond_S) * diag(1,r,r)
    }
  }
  return(list(mu_hat = mu_hat, S_hat = S_hat, t=t, R=R))
}
