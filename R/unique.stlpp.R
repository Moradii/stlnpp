#' @export
unique.stlpp <- function(X,...){
  if (!inherits(X, "stlpp")) stop("X should be from class stlpp")

  Y <- unique(as.data.frame(X$data),...)
  X <- as.stlpp(Y$x,Y$y,Y$t,L=domain(X))
  return(X)
}
