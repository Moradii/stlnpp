#' @export
unique.stlpp <- function(x,...){
  if (!inherits(x, "stlpp")) stop("X should be from class stlpp")

  Y <- unique(as.data.frame(x$data),...)
  X <- as.stlpp(Y$x,Y$y,Y$t,L=domain(x))
  return(X)
}
