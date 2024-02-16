#' Inhomogeneous K-function for spatio-temporal  point processes on linear networks
#'
#' This function computes the inhomogeneous K-function for spatio-temporal  point patterns on linear networks.
#'
#' @usage STLKinhom(X,lambda=lambda,normalize=FALSE,r=NULL,t=NULL,nxy=10)
#'
#' @param X a spatio-temporal point pattern of class \code{\link{stlpp}}
#' @param lambda values of estimated intensity at data points
#' @param normalize normalization factor to be considered
#' @param r values of argument r where pair correlation function will be evaluated. optional
#' @param t values of argument t where pair correlation function will be evaluated. optional
#' @param nxy pixel array dimensions. optional
#' 
#' @seealso \code{\link{STLg}}, \code{\link{STLK}}, \code{\link{STLginhom}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com> 
#' 
#' @returns 
#' An object of class \code{sumstlpp}.
#'
#' @details 
#' This function calculates the inhomogeneous K-function for a spatio-temporal point patterns on a linear network.
#' 
#' @references Moradi, M., & Mateu, J. (2020). First-and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics, 29(3), 432-443.
#' 
#' 
#' @examples
#' X <- rpoistlpp(.2,a=0,b=5,L=easynet)
#' lambda <- density(X,at="points")
#' k <- STLKinhom(X,lambda=lambda,normalize=TRUE)
#' plot(k)
#' 
#' 
#' @export
STLKinhom <- function(X,lambda=lambda,normalize=FALSE,r=NULL,t=NULL,nxy=10){

  if (!inherits(X, "stlpp")) stop("X should be from class stlpp")


  Y <- as.lpp.stlpp(X)
  l <- domain(Y)
  tleng <- summary(l)$totlength
  n <- npoints(X)
  a <- X$time[1]
  b <- X$time[2]
  trange <- b-a
  timev <- X$data$t

  sdist <- pairdist.lpp(Y)
  tdist <- as.matrix(dist(timev))

  toler <- default.linnet.tolerance(l)
  ml <- matrix(1, n, n)
  for(j in 1:n) {
    ml[ -j, j] <- countends(l, Y[-j], sdist[-j,j], toler=toler)
  }

  mtplus <- matrix(timev,n,n,byrow = T)+tdist
  mtminus <- matrix(timev,n,n,byrow = T)-tdist
  mtedge <- (mtplus<=b) + (mtminus>=a)
  diag(mtedge) <- 1

  lamden <- outer(lambda,lambda,FUN = "*")
  diag(lamden) <- 1
  edgetl <- mtedge*ml*lamden

  maxs <- 0.7*max(sdist[!is.infinite(sdist)])
  maxt <- 0.7*(trange/2)

  if(is.null(r)) r <- seq((maxs/nxy),maxs,by=(maxs-(maxs/nxy))/(nxy-1))
  if(is.null(t)) t <- seq((maxt/nxy),maxt,by=(maxt-(maxt/nxy))/(nxy-1))

  K <- matrix(NA,nrow = nxy,ncol = nxy)

  for (i in 1:length(r)) {

    for (j in 1:length(t)) {
      out <- (sdist<=r[i])*(tdist<=t[j])
      diag(out) <- 0
      kout <- out/edgetl
      K[i,j] <- sum(kout[!is.na(kout) & !is.infinite(kout)])
    }
  }

  if(normalize){
    revrho <- outer(1/lambda,1/lambda,FUN = "*")
    appx <- (tleng*trange)/(sum(revrho[lower.tri(revrho, diag = FALSE)])*2)
    K <- K*appx
  }
  else{
    K <- K/(trange*tleng)
  }

  ##
  pixcor <- expand.grid(r,t)
  Kout <- list(Kinhom=K,Ktheo=matrix(pixcor[,1]*pixcor[,2],ncol = nxy),r=r,t=t)
  class(Kout) <- c("sumstlpp")
  attr(Kout,"nxy") <- nxy
  return(Kout)
}
