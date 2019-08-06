#' @export
as.stlpp.lpp <- function(X){
  if(!any(class(X)=="stlpp")) stop("class(X) must be stlpp")
  Y <- as.lpp(x=X$data$x,y=X$data$y,L=X$domain,sparse = F)
  return(Y)
}
