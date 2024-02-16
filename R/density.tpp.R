#' Kernel estimation of intensity of one-dimensional point patterns
#'
#' Kernel estimation of intensity of one-dimensional point patterns.
#'
#' @usage \method{density}{tpp}(x,tbw,at=c("points","pixels"),...)
#'
#' @param x an object of class \code{\link{tpp}}
#' @param tbw time smoothing bandwidth
#' @param at string specifying whether to compute the intensity values at a grid of pixel locations (at="pixels") or only at the points of x (at="points"). default is to estimate the intensity at pixels
#' @param ... arguments passed to \link[stats]{density}
#' 
#' @seealso \code{\link{density}}, \code{\link{bw.nrd0}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com> and Ottmar Cronie
#' 
#' @returns 
#' If \code{at="points"}: a vector of intensity values at the data points of x.
#'
#' If \code{at="pixels"}: a vector of intensity values over a grid.
#'
#' @details 
#' A vector of intensity values.
#' 
#' 
#' @references Mateu, J., Moradi, M., & Cronie, O. (2019). Spatio-temporal point patterns on linear networks: Pseudo-separable intensity estimation. Spatial Statistics, 100400.
#' 
#' 
#' @examples  
#' X <- tpp(sample(c(1:24),200,replace = TRUE))
#' plot(density(X))
#' 
#' 
#' @export
density.tpp <- function(x,tbw,at=c("points","pixels"),...){
  
  if (!inherits(x, "tpp")) stop("x should an object of class tpp")
  
  n <- npoints(x) # Emerge number of points
  
  if (missing(tbw)) {
    d <- density(x$data$t,...)
    }
  else {
    d <- density(x,bw=tbw,...)
  }
  
  if(missing(at)) at <- "pixels"
  
  if(at=="points"){
    Tint <- d$y[findInterval(x$data$t, d$x)] * npoints(x)
  }
  if(at=="pixels"){
    Tint <- d$y * npoints(x)
  }
  
  out <- Tint
  
  attr(out,"tempden") <- d
  attr(out,"bw") <- d$bw
  attr(out,"time") <- x$data$t
  attr(out,"tpp") <- x
  attr(out,"tgrid") <- d$x
  if(at=="points"){
    class(out) <- c("numeric")  
  }else{
    class(out) <- c("tppint")
  }
  
  return(out)
}