#' Inhomogeneous empty space function for spatio-temporal  point processes on linear networks
#'
#' This function computes the inhomogeneous empty space function for spatio-temporal  point patterns on linear networks.
#'
#' @usage STLFinhom(X,
#' lambda = NULL,
#' r = NULL,
#' t = NULL,
#' nxy = 10,
#' dimt = 10,
#' dimyx = 12)
#'
#' @param X a spatio-temporal point pattern of class \code{\link{stlpp}}
#' @param lambda values of estimated intensity at data points
#' @param r values of argument r where pair correlation function will be evaluated. optional
#' @param t values of argument t where pair correlation function will be evaluated. optional
#' @param nxy pixel array dimensions. optional
#' @param dimt description
#' @param dimyx description
#' 
#' @seealso \code{\link{STLJinhom}}, \code{\link{STLHinhom}}, \code{\link{STLKinhom}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com> 
#' 
#' @returns 
#' An object of class \code{sumstlpp}.
#'
#' @details 
#' This function calculates the inhomogeneous empty space function for a spatio-temporal point patterns on a linear network.
#' 
#' @references Moradi, M., & Sharifi, A. (2024). Summary statistics for spatio-temporal point processes on linear networks. Spatial Statistics, 61, 100840.
#' 
#' 
#' @examples
#' \dontrun{
#' X <- rpoistlpp(.2,a=0,b=5,L=easynet)
#' lambda <- density(X,at="points")
#' k <- STLKinhom(X,lambda=lambda,normalize=TRUE)
#' plot(k)
#' }
#' 
#' @export
STLFinhom <- function(X,
                      lambda = NULL,
                      r = NULL,
                      t = NULL,
                      nxy = 10,
                      dimt = 10,
                      dimyx = 12
                      ){
  
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
  
  grid <- grid.stlpp(L, a, b, dimt = dimt, dimyx = dimyx)
  timegrid <- grid$data$t
  
  dc <- crossdist.lpp(Y, as.lpp.stlpp(grid))
  
  # intensity at data/grid points
  lambdap <- lambda  # intensity at data points
  if(!(length(lambdap)==n) | anyNA(lambdap)) stop("some intensity values are missed or NA intensity value")
  
  lambdabar <- min(lambdap)
  
  # border extraction on the network
  Yborder <- terminalvertices(L)
  
  dcrossgrid <- crossdist.lpp(as.lpp.stlpp(grid), Yborder) # distance of grid points to border
  gridtoborder <- apply(dcrossgrid, 1, FUN=min)
  
  timegridtoborder <- data.frame(lowdist=grid$data$t-a,updist=b-grid$data$t)
  timegridtoborder <- apply(abs(timegridtoborder), 1, min)
  
  tXG <- t(abs(sapply(timev,"-",grid$data$t)))
  
  Flpp <- lapply(X=1:length(r), function(i){
    out <- c()
    
    for (j in 1:length(t)) {
      #print(paste0("Flpp:",i,"-",j))
      
      sum <- 0
      
      OK <- (gridtoborder > r[i]) & (timegridtoborder > t[j]) # points after minus sampling
      F_seq_OK <- which(OK) # points after minus sampling
      
      if (sum(OK)==0) {
        out[j] <- 0
      }else{
        for (h in F_seq_OK) {
          row <- dc[,h]
          pr <- which(0<row & row<=r[i])
          
          rowt <- as.vector(tXG[,h])
          pt <- which(0<rowt & rowt<=t[j])
          
          p <- intersect(pr,pt)
          
          if(length(p)>=1){
            w <- countends.stlpp(grid[h],r=dc[p,h],t=as.vector(tXG[p,h]))
            frac <- 1 - (lambdabar/(lambdap[p]*as.numeric(unlist(w))))
            # frac <- 1 - (1/(as.numeric(unlist(w))))
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
  
  out.f <- 1 - do.call(rbind, Flpp)
  Ftheo <- 1 - exp(-lambdabar*outer(r,t))
  

  Fout <- list(Finhom = out.f, Ftheo = Ftheo, r = r, t = t)
  
  
  class(Fout) <- c("sumstlpp")
  attr(Fout,"nxy") <- nxy
  return(Fout)
}
