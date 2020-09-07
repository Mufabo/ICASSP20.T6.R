#### helper ####

'%\\%' <- function(A, B) pracma::mldivide(A, B)

'%/%' <- function(A, B) pracma::mrdivide(A, B)

matSum <- function(tensor){
  return(Reduce('+', comprehenr::to_list(for (i in 1:dim(tensor)[3]) tensor[,,i])))
}
