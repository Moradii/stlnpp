#' Intensity estimate of temporal point patterns using Voronoi-Dirichlet tessellation
#'
#' This function performs adaptive intensity estimation for temporal point patterns using Voronoi-Dirichlet tessellation.
#'
#' @usage \method{densityVoronoi}{tpp}(X, f = 1, nrep = 1, at=c("points","pixels"), dimt=128,...)
#'
#' @param X an object of class \code{\link{tpp}}
#' @param f fraction (between 0 and 1 inclusive) of the data points that will be used to build a tessellation for the intensity estimate
#' @param nrep number of independent repetitions of the randomised procedure
#' @param at string specifying whether to compute the intensity values at a grid of pixel locations and time (at="pixels") or only at the points of x (at="points"). default is to estimate the intensity at pixels
#' @param dimt the number of equally spaced points at which the temporal density is to be estimated. see \link[stats]{density}
#' @param ... arguments passed to \code{\link{densityVoronoi.lpp}}
#' 
#' @seealso \code{\link{densityVoronoi.lpp}}, \code{\link{density.stlpp}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com> and Ottmar Cronie
#' 
#' @returns 
#' If \code{at="points"}: a vector of intensity values at the data points of X.
#'
#' If \code{at="pixels"}: a vector of intensity values over a grid.
#'
#' @details 
#' This function computes intensity estimates for temporal point patterns using Voronoi-Dirichlet tessellation.
#' 
#' IF f<1, then nrep independent sub-samples of X are obtained using the function \code{\link{rthin.stlpp}}. Then for each of the obtained sub-samples, we calculate the Voronoi estimate. The final estimation is the sum of all obtained estimated intensities divided by (f*nrep). 
#' 
#' 
#' @references Mateu, J., Moradi, M., & Cronie, O. (2019). Spatio-temporal point patterns on linear networks: Pseudo-separable intensity estimation. Spatial Statistics, 100400.
#' 
#' 
#' @examples  
#' X <- rpoistlpp(0.2,a=0,b=5,L=easynet)
#' Y <- as.tpp.stlpp(X)
#' densityVoronoi(Y)
#' 
#' 
#' @export
densityVoronoi.tpp <- function(X, f = 1, nrep = 1,
                               at=c("points","pixels"),
                               dimt=128,...){
  
  if(!inherits(X, "tpp")) stop("X should an object of class tpp")
  
  if(missing(at)) at <- "pixels"
  
  n <- npoints(X)
  
  if(f<0 | f>1) stop("f should be between 0 and 1")
  
  Xt <- lpp(X=cbind(X$data$t,rep(0,n)), 
            L=linnet_interval(startp=X$time[1], endp=X$time[2]))
  
  out <- densityVoronoi.lpp(Xt,f=f,nrep=nrep,dimyx=dimt,...)
  
  if(at=="pixels"){
    out1 <- out$v[!is.na(out$v)]
    class(out1) <- c("tppint")
    attr(out1,"time") <- X$data$t
    
    attr(out1,"tgrid") <- out$xcol
    attr(out1,"tpp") <- X
    
    return(out1)
    }
  else{
    Tint <- out$v[!is.na(out$v)]
    ID <- findInterval(X$data$t, out$xcol)
    ID[which(ID==0)] <- 1
    Tintout <- Tint[ID]
    
    class(Tintout) <- c("numeric")
    attr(out,"time") <- X$data$t
    
    return(Tintout)
  }
}


