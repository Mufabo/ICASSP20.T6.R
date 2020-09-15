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
  if(r == 1){
    # data is Nx1
    graphics::plot.new()
    graphics::plot.window(c(-20, 20), c(-20, 15))
    graphics::axis(1)
    graphics::axis(2)
    for(k in 1:(K_true+1)){
      if(k == (K_true+1)){
        graphics::points(data[data[,1]==k, 2], numeric(sum(data[,1]==k)))
      } else {
        graphics::points(data[data[,1]==k, 2], numeric(sum(data[,1]==k)), pch = 16, col = k+1)
      }
    }
  }
  if(r == 2){
    # data is Nx1
    graphics::plot.new()
    graphics::plot.window(c(-20, 20), c(-20, 15))
    graphics::axis(1)
    graphics::axis(2)
    for(k in 1:(K_true+1)){
      if(k == (K_true+1)){
        graphics::points(data[data[,1]==k, 2], data[data[,1]==k, 3])
      } else {
        graphics::points(data[data[,1]==k, 2], data[data[,1]==k, 3], pch = 16, col = k+1)
      }
    }
  }
  if(r == 3){
    print("Implementation on demand ...")
    }
}
