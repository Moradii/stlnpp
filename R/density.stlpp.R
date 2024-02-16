#' Kernel estimation of intensity of spatio-temporal point patterns on a linear network
#'
#' Kernel density estimation of a spatio-temporal point pattern on a linear network.
#'
#' @usage \method{density}{stlpp}(x,lbw,tbw,at=c("points","pixels"),dimt=512,...)
#'
#' @param x an object of class \code{\link{stlpp}}
#' @param lbw network smoothing bandwidth
#' @param tbw time smoothing bandwidth
#' @param at string specifying whether to compute the intensity values at a grid of pixel locations and times (at="pixels") or only at the points of x (at="points"). default is to estimate the intensity at pixels
#' @param dimt the number of equally spaced points at which the temporal density is to be estimated. see \link[stats]{density}
#' @param ... arguments passed to \code{\link{density.lpp}}
#' 
#' @seealso \code{\link{density}}, \code{\link{density.lpp}}, \code{\link{bw.nrd0}}, \code{\link{bw.scott.iso}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com>
#' 
#' @returns 
#' If \code{at="points"}: a vector of intensity values at the data points of x.
#' If \code{at="pixels"}: a list of images on linear network. Each image represents an estimated spatio-temporal intensity at a fixed time.
#' Check the attributes for more accommodated outputs.
#' 
#' @details Kernel smoothing is applied to the spatio-temporal point pattern x using methods in Moradi et al (2019). The function computes estimated intensities assuming first-order separability. Estimated intensity values of the marginal spatial point pattern on the linear network will be obtained using the fast kernel smoothing technique of Rakshit et al. (2019) and function  \code{\link{densityQuick.lpp}}, whereas the estimated intensity values of the marginal temporal point pattern will be estimated using the function \code{\link{density}}.
#'
#' If lbw and tbw are not given, then they will be selected using \code{\link{bw.nrd0}} and \code{\link{bw.scott.iso}} respectively.
#' 
#' @references Moradi, M., & Mateu, J. (2020). First-and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics, 29(3), 432-443.
#' @examples  
#' X <- rpoistlpp(.2,a=0,b=5,L=easynet)
#' density(X)
#'
#' @import spatstat
#' @import spatstat.geom
#' @import spatstat.linnet
#' @import stats
#' @export
density.stlpp <- function(x,lbw,tbw,at=c("points","pixels"),dimt=512,...){
  
  if (!inherits(x, "stlpp")) stop("x should an object of class stlpp")
  
  if(missing(at)) at <- "pixels"
  ox <- x$data$x
  oy <- x$data$y
  ot <- x$data$t
  
  L <-  x$domain
  
  n <- npoints(x) # Emerge number of points
  
  stint <- 0 # define the vacant vector to save the densities values per points
  
  if (missing(tbw)) {
    d <- density(ot,n=dimt)
  }
  else{
    d <- density(ot,bw=tbw,n=dimt)
  }
  
  if(at=="points"){
    Tint <- d$y[findInterval(ot, d$x)] * n
  }
  else{
    Tint <- d$y * n
  }
  ############################################## space intensity
  
  pX <- as.lpp.stlpp(x)
  if (missing(lbw)) lbw <- bw.scott.iso(pX)
  
  if(at=="points"){
    ldens <- density.lpp(pX,sigma = lbw,distance="euclidean",...)
    Sint <- density.lpp(pX,sigma = lbw,at="points",distance="euclidean",...)
    stint <-  Sint*Tint/npoints(pX)
    out <- stint
  }
  
  if(at=="pixels"){
    
    ldens <- density.lpp(pX,sigma = lbw,distance="euclidean",...)
    out <- lapply(X=1:length(Tint), function(i){
      ldens*Tint[i]/npoints(pX)
    })
  }
  
  
  names(lbw) <- NULL
  if(at=="points") class(out) <- c("numeric")
  if(at=="pixels") class(out) <- c("list","stlppint")
  
  attr(out,"tempden") <- d
  attr(out,"netint") <- ldens
  attr(out,"time") <- ot
  attr(out,"bw") <- c("sigma_l"=lbw,"sigma_t"=d$bw)
  attr(out,"stlpp") <- x
  
  
  return(out)
  
}
