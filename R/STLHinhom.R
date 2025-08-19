#' Inhomogeneous nearest neighbour function for spatio-temporal point processes on linear networks
#'
#' This function computes the inhomogeneous nearest neighbour function for spatio-temporal point processes on linear networks.
#'
#' @usage STLHinhom(X,
#' lambda = NULL,
#' r = NULL,
#' t = NULL,
#' nxy = 10)
#'
#' @param X a spatio-temporal point pattern of class \code{\link{stlpp}}
#' @param lambda values of estimated intensity at data points
#' @param r values of argument r where pair correlation function will be evaluated. optional
#' @param t values of argument t where pair correlation function will be evaluated. optional
#' @param nxy pixel array dimensions. optional
#' 
#' @seealso \code{\link{STLFinhom}}, \code{\link{STLKinhom}}, \code{\link{STLginhom}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com> 
#' 
#' @returns 
#' An object of class \code{sumstlpp}.
#'
#' @details 
#' This function calculates the inhomogeneous nearest neighbour function for a spatio-temporal point patterns on a linear network.
#' 
#' @references Moradi, M., & Sharifi, A. (2024). Summary statistics for spatio-temporal point processes on linear networks. Spatial Statistics, 61, 100840.
#' 
#' 
#' @examples
#' \dontrun{
#' X <- rpoistlpp(0.2, a = 0, b = 5, L = easynet)
#' lambda <- density(X, at = "points")
#' h <- STLHinhom(X, lambda = lambda)
#' }
#' 
#' 
#' @export
STLHinhom <- function(X,
                      lambda = NULL,
                      r = NULL,
                      t = NULL,
                      nxy = 10){
  
  if (!inherits(X, "stlpp")) stop("X should be of class stlpp")
  
  Y <- as.lpp.stlpp(X)
  timev <- X$data$t
  
  n <- npoints(X)
  L <- as.linnet(Y)  
  
  a <- X$time[1]
  b <- X$time[2]
  trange <- b-a
  
  # r calculation
  if(is.null(r)) {
    rmaxdefault <- 0.5 * boundingradius(L)
    r <- seq(0,rmaxdefault,(rmaxdefault)/(nxy-1))
  }
  
  # t calculation 
  if(is.null(t)){
    tmax <- trange/2
    t <- seq(0, tmax, by = tmax/(nxy-1))
  }
  
  tdist <- as.matrix(dist(timev))
  
   
  
  tleng <- summary(L)$totlength
  
  dp <- pairdist.lpp(Y)
  
  # intensity at data/grid points
  lambdap <- lambda  # intensity at data points
  if(!(length(lambdap)==n) | anyNA(lambdap)) stop("some intensity values are missed or NA intensity value")
  
  lambdabar <- min(lambdap)
  
  Yborder <- terminalvertices(L)
  
  dcrossY <- crossdist.lpp(Y,Yborder)
  Ytoborder <- apply(dcrossY, 1, FUN=min)
  
  timetoborder <- data.frame(lowdist=timev-a,updist=b-timev)
  timetoborder <- apply(abs(timetoborder), 1, min)
  
  
  Hlpp <- lapply(X=1:length(r), function(i){
    out <- c()
    
    for (j in 1:length(t)) {
    
      sum <- 0
      OK <- (Ytoborder > r[i]) & (timetoborder>t[j]) # points after minus sampling
      H_seq_OK <- which(OK) # points after minus sampling
      
      if (sum(OK)==0) {
        out[j] <- 0
      }else{
        for (h in H_seq_OK) {
          row <- dp[,h]
          pr <- which(0<row & row<=r[i])
          
          rowt <- as.vector(tdist[,h])
          pt <- which(0<rowt & rowt<=t[j])
          
          p <- intersect(pr,pt)
          
          if(length(p)>=1){
            w <- countends.stlpp(X[h],r=dp[p,h],t=as.vector(tdist[p,h]))
            frac <- 1 - (lambdabar/(lambdap[p]*as.numeric(unlist(w))))
            frac[frac==Inf] <- 0
            sum <- sum + prod(frac)
          }else{
            sum <- sum+1
          }
        }
        out[j] <- (sum/sum(OK))
      }
    }
    return(out)
  })

  out.H <- 1 - do.call(rbind, Hlpp)
  
  Htheo <- 1 - exp(-lambdabar*outer(r,t))
  
  Hout <- list(Hinhom = out.H, Htheo = Htheo, r = r, t = t)
  
  class(Hout) <- c("sumstlpp")
  attr(Hout,"nxy") <- nxy
  return(Hout)
}
