#' Intensity estimate of spatio-temporal point pattern using Voronoi-Dirichlet tessellation
#'
#' This function performs adaptive intensity estimation for spatio-temporal point patterns on linear networks using Voronoi-Dirichlet tessellation.
#'
#' @usage \method{densityVoronoi}{stlpp}(X, f = 1, nrep = 1, separable=FALSE, at=c("points","pixels"), dimt=128,...)
#'
#' @param X an object of class \code{\link{stlpp}}
#' @param f fraction (between 0 and 1 inclusive) of the data points that will be used to build a tessellation for the intensity estimate
#' @param nrep number of independent repetitions of the randomised procedure
#' @param separable logical. If FALSE, it then calculates a pseudo-separable estimate
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
#' If \code{at="pixels"}: a list of images on a linear network. Each image represents an estimated spatio-temporal intensity at a fixed time.
#' 
#' @details 
#' This function computes intensity estimates for spatio-temporal point patterns on linear networks using Voronoi-Dirichlet tessellation. Both first-order separability and pseudo-separability assumptions are accommodated in the function.
#'
#' If separable=TRUE, the estimated intensities will be a product of the estimated intensities on the network and those on time. Estimated intensity of the spatial component will be obtained using \code{\link{densityVoronoi.lpp}}, whereas estimated intensities of the temporal component will be obtained via \code{\link{densityVoronoi.tpp}}. If f=1, the function calculates the estimations based on the original Voronoi intensity estimator.
#'
#' If separable=FALSE, the estimated intensities will be calculated based on a sub-sampling technique explained in Mateu et al. (2019). nrep sub-samples will be obtained from X based on a given retention probability f, the function \code{\link{densityVoronoi.stlpp}}, considering separable=TRUE and f=1, will be applied to each obtained sub-sample, and finally, the estimated intensities will be the sum of all obtained estimated intensities from all sub-samples divided by the (f * nrep).
#' 
#' @references Mateu, J., Moradi, M., & Cronie, O. (2019). Spatio-temporal point patterns on linear networks: Pseudo-separable intensity estimation. Spatial Statistics, 100400.
#' 
#' 
#' @examples  
#' X <- rpoistlpp(.2,a=0,b=5,L=easynet)
#' densityVoronoi(X)
#' 
#' 
#' @import spatstat
#' @import spatstat.geom
#' @import spatstat.linnet
#' @importFrom spatstat.random rthin runifpointOnLines
#' @importFrom spatstat.explore densityVoronoi bw.scott.iso
#' @import stats
#' @export
densityVoronoi.stlpp <- function(X, f = 1, nrep = 1,
                                 separable=FALSE,at=c("points","pixels"),
                                 dimt=128,...){
  
  if(!inherits(X, "stlpp")) stop("X should an object of class stlpp")
  
  if(missing(at)) at <- "pixels"

  n <- npoints(X)
  
  if(f<0 | f>1) stop("f should be between 0 and 1")
  
  if(f==1) separable <- TRUE
  
  Xt <- lpp(X=cbind(X$data$t,rep(0,n)), 
            L=linnet_interval(startp=X$time[1], endp=X$time[2]))
  
  Xs <- as.lpp.stlpp(X)
  
if(separable){
    
    IntEstL <- densityVoronoi(X=Xs,f=f,nrep=nrep,...)
    IntEstT <- densityVoronoi(X=Xt,f=f,nrep=nrep,dimyx=dimt,...)
    
    tgrid <- IntEstT$xcol
    
    IntEstTv <- IntEstT$v[!is.na(IntEstT$v)]/n
    
    
    out <- lapply(X=(1:length(IntEstTv)), FUN=function(j){IntEstL*IntEstTv[j]}) 
    

}else{
    
Y <- rthin(X,P=f,nsim = nrep)
  
out.nonsep <- lapply(X=1:nrep, function(i){
                     densityVoronoi(Y[[i]],f=1,nrep = 1, separable = TRUE,dimt=dimt,...)
                   })
tgrid <- attr(out.nonsep[[1]],"tgrid")
  
out <- list()
k <- length(out.nonsep[[1]])
for (i in 1:k){
  
    out[[i]] <- out.nonsep[[1]][[1]]
    out[[i]]$v[!is.na(out[[i]]$v)] <- 0
   
for (j in 1:length(out.nonsep)){
      out[[i]] <- out[[i]] + out.nonsep[[j]][[i]]
                                     
                               }
    out[[i]] <- out[[i]]/(f*nrep)
                                     }
    }
  
if(at=="points"){
    t <- X$data$t
    id <- findInterval(t,tgrid)
    out1 <- c()
    for (i in 1:n){
      out1[i] <- out[[id[i]]][Xs[i]]
                  }
    out <- out1
                }

  if(at=="points") class(out) <- c("numeric","stlppint")
  if(at=="pixels") class(out) <- c("list","stlppint")
  
  if(separable){
    attr(out,"lint") <- IntEstL
    attr(out,"tint") <- IntEstT$v[!is.na(IntEstT$v)]
  }
  attr(out,"tgrid") <- tgrid
  attr(out,"time") <- X$data$t
  attr(out,"stlpp") <- X
return(out)
}




linnet_interval <- function(startp=0, endp=1,...){
  Wt <- c(startp, endp)
  vertices <- ppp(x=c(Wt[1],Wt[2]), y=c(0,0), 
                  window=owin(Wt,c(-0.5,0.5)))
  m <- matrix(data=c(FALSE,TRUE,TRUE,FALSE), nrow=2, ncol=2, 
              byrow=TRUE)
  out <- linnet(vertices=vertices, m=m,...)
  return(out)
}
