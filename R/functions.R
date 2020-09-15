#### helper ####

'%\\%' <- function(A, B) pracma::mldivide(A, B)

'%/%' <- function(A, B) pracma::mrdivide(A, B)

#' sum of matrices of a rank 3 tensor
#'
#' @param tensor
#'
#' @return
#' @export
#'
#' @examples
matSum <- function(tensor){
  return(Reduce('+', comprehenr::to_list(for (i in 1:dim(tensor)[3]) tensor[,,i])))
}

#' Title
#'
#' @param data
#' @param ll
#' @param iter
#' @param rep
#'
#' @return
#' @export
#'
#' @examples
kmeanspp <- function(data, ll, iter = 25, rep = 5){
  best <- NULL
  for(i in 1:rep){
    init_centers <- ClusterR::KMeans_rcpp(data, ll)$centroids
    # Remove NaN's from init_centers
    init_centers[is.nan(init_centers)] <- 0
    km <- amap::Kmeans(data, init_centers, iter.max = iter, method="manhattan")
    if(is.null(best) || sum(km$withinss) < sum(best$withinss)) best <- km
  }

  return(list(centroids = best$centers, clusters = best$cluster))
}
