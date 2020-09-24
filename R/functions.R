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
    #set.seed(123)
    init_centers <- LICORS::kmeanspp(data, k= ll)
    km <- amap::Kmeans(data, init_centers$centers, method="manhattan")
    # To avoid empty cluster, only accept solutions where each cluster has
    # been allocated at least one member
    if(is.null(best) || sum(km$withinss) < sum(best$withinss) &
       !(any(is.na(km$centers)) | any(is.nan(km$centers)))) best <- km
  }

  return(list(centroids = best$centers, clusters = best$cluster))
}

# Function to plot color bar from
# https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)

  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
