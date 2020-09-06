#### nearestSPD ####

#' Finds closest symmetric positive definite matrix
#'
#' @param A matrix
#'
#' @return Ahat closest SPD matrix to A
#' @export
#'
#' @examples
#'
#' m <- matrix(c(1,2,3,4,5,6,7,8,9),3,3,byrow = TRUE)
#' nearestSPD(m)
nearestSPD <- function(A){
  if(dim(A)[1] != dim(A)[2]){
    stop("A is not a square matrix!")
  }

  # if A is scalar, return
  if(is.null(dim(A)) && A <= 0) return(.Machine$double.eps)

  # symmetrize A into B
  B <- (A + t(A)) * .5

  # Compute the symmetric polar factor of B. Call it H.
  # Clearly H is itself SPD.
  duv <- svd(B)
  H <- duv$v %*% diag(duv$d) %*% t(duv$v)
  Ahat <- (B + H) / 2
  Ahat <- (Ahat + t(Ahat)) / 2 # Ensure symmetry
  # Ensure that Ahat is PD
  p <- FALSE
  k <- 0
  while(!p){
    vals <- eigen(Ahat)$values
    p <- all(vals > 0)
    k <- k + 1
    if(!p){
      mineig <- min(vals)
      Ahat <- Ahat + (-mineig * k^2 + .Machine$double.eps + 10^(-10)) * diag(1 ,dim(A))
    }
  }


  return(Ahat)
}
