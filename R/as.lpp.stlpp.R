#' @export
as.lpp.stlpp <- function(x,...){
  if(!any(class(x)=="stlpp")) stop("class(x) must be stlpp")
  Y <- as.lpp(x=x$data$x,y=x$data$y,L=x$domain,sparse = F,...)
  return(Y)
}
