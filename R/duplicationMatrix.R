#### duplicationMatrix ####

#' The duplication matrix
#'
#' @param num_features int. Number of features
#'
#' @return dupmat Matrix[num_features^2, .5*num_features*(num_features+1)] The duplication matrix
#'
#' @note
#' [1] F. K. Teklehaymanot, M. Muma, and A. M. Zoubir, "Bayesian Cluster Enumeration Criterion for Unsupervised Learning",
#' IEEE Trans. Signal Process. (accepted),
#' [Online-Edition: https://arxiv.org/abs/1710.07954v2], 2018.
#'
#' [2] F. K. Teklehaymanot, M. Muma, and A. M. Zoubir, "Novel Bayesian Cluster
#'     Enumeration Criterion for Cluster Analysis With Finite Sample Penalty Term",
#' in Proc. 43rd IEEE Int. conf. on Acoustics, Speech and Signal Process. (ICASSP), pp. 4274-4278, 2018,
#' [Online-Edition: https://www.researchgate.net/publication/322918028]
#'
#'
#' Copyright (c) 2018 Freweyni K. Teklehaymanot. All rights reserved.
#' @export
duplicationMatrix <- function(num_features){

  dupmatTranspose <- matrix(data = 0, nrow = .5*num_features*(num_features+1), ncol = num_features^2)

  for (j in 1:num_features) {
    for (i in 1:num_features) {
      u <- numeric(.5*num_features*(num_features+1))
      idx <- (j-1)*num_features+i-.5*j*(j-1)
      u[idx] <- 1

      Y <- matrix(0, nrow = num_features, ncol = num_features)
      Y[i, j] <- 1
      Y[j, i] <- 1
      d <- u %o% c(Y)

      if(i >= j) dupmatTranspose <- dupmatTranspose + d

    }
  }

  return(t(dupmatTranspose))
}
