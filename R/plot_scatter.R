#' Title
#'
#' @param data
#' @param K_true
#' @param r
#'
#' @return
#' @export
#'
#' @examples
plot_scatter <- function(data, K_true, r){
  for(k in 1:(K_true+1)){
    if(r == 1){
      if(k == (K_true+1)){
        graphics::plot(data[data[,1]==k, 2], numeric(sum(data[,1]==k)))
      } else {
        graphics::plot(data[data[,1]==k, 2], numeric(sum(data[,1]==k)), pch = 16)
      }
    } else if(r==2){
      if(k == (K_true + 1)){
        graphics::plot(data[data[,1]==k, 2], data[data[,1]==k, 3])
      } else graphics::plot(data[data[,1]==k, 2], data[data[,1]==k, 3], pch=16)
    } else if(r==3){
      if(k == (K_true + 1)){
        scatterplot3d::scatterplot3d(data[data[,1]==k, 3], data[data[,1]==k, 4])
      } else scatterplot3d::scatterplot3d(data[data[,1]==k, 3], data[data[,1]==k, 4], pch = 16)
    }

  }
}
