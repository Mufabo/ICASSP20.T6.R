#### helper ####


# See pracma::mldivide
'%\\%' <- function(A, B) pracma::mldivide(A, B)


# See pracma::mrdivide
'%/%' <- function(A, B) pracma::mrdivide(A, B)

#' sum of matrices of a rank 3 tensor
#'
#' @param tensor 3darray
#'
#' @return  Sum of matrices in tensor
#' @export
#'
matSum <- function(tensor){
  return(Reduce('+', comprehenr::to_list(for (i in 1:dim(tensor)[3]) tensor[,,i])))
}

#' Kmeans with kmeans++ initialization and cityblock distance
#'
#' @param data matrix. Data
#' @param ll int. Number of clusters
#' @param iter int. Number of iterations
#' @param rep int. Number of replications
#'
#' @return list
#' \enumerate{
#' \item centroids matrix. Cluster centers
#' \item clusters vector. Cluster membership of data samples
#' }
#' @export
#'
#' @examples
kmeanspp <- function(data, ll, iter = 25, rep = 50){
  best <- NULL
  for(i in 1:rep){
    set.seed(123)
    init_centers <- LICORS::kmeanspp(data, k= ll)
    km <- amap::Kmeans(data, init_centers$centers, method="manhattan")
    # To avoid empty cluster, only accept solutions where each cluster has
    # been allocated at least one member
    if(is.null(best) || sum(km$withinss) < sum(best$withinss) &
       !(any(is.na(km$centers)) | any(is.nan(km$centers)))) best <- km
  }

  return(list(centroids = best$centers, clusters = best$cluster))
}
