#' Simulating spatio-temporal Poisson point processes on a linear network
#'
#' This function simulates realisations of a spatio-temporal Poisson point process on a linear network.
#'
#' @usage rpoistlpp(lambda,a,b,L,check=FALSE,lmax=NULL,nsim=1)
#'
#' @param lambda intensity of the point process. it can be either a number, function of location and time, or an abject of class \code{stlppint}
#' @param a lower bound of time period
#' @param b upper bound of time period
#' @param L a linear network
#' @param check logical value indicating whether to check that all the (x,y) points lie inside the specified window. see \code{\link{ppp}}
#' @param lmax upper bound for the values of \code{labmda}. this is optional
#' @param nsim number of simulated patterns to generate
#' 
#' @seealso \code{\link{density.stlpp}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com> 
#' 
#' @returns 
#' an object of class \code{\link{stlpp}} if nsim=1, otherwise a list of objects of class \code{\link{stlpp}}.
#'
#' @details 
#' This function generates realisations of a spatio-temporal poisson point process on a linear network based on an intensity function lambda and lower/upper bounds a and b.
#' 
#' 
#' @references Moradi, M., & Mateu, J. (2020). First-and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics, 29(3), 432-443.
#' 
#' 
#' @examples  
#' X <- rpoistlpp(0.2,a=0,b=5,L=easynet)
#' X
#' 
#' @export
rpoistlpp <-  function(lambda,a,b,L,check=FALSE,lmax=NULL,nsim=1){
  
    
  if (!is.numeric(lambda) & !is.function(lambda) & !any(class(lambda)=="stlppint"))
    stop(" lambda should be a number, a function, or an object of class stlppint")
  
  if(nsim > 1) {
    out <- list()
    for (i in 1:nsim) {
      out[[i]] <- rpoistlpp(lambda,a=a,b=b,L=L,nsim=1,check=check,lmax=lmax)
    }
    return(out)
  }
  
  if (is.numeric(lambda)){
    
    if (!inherits(L,"linnet"))   stop("L should be a linear network")
    
    if (a >= b)  stop("lower bound must be smaller than upper bound")
    
    
    n <- rpois(1,lambda*volume(L)*(b-a))
    X <- runiflpp(n,L)
    t <- runif(npoints(X),a,b)
    stlpp <- data.frame(x=X$data$x,y=X$data$y,t)
  }
  else if(is.function(lambda)){
    
    if (!inherits(L,"linnet"))   stop("L should be a linear network")
    
    if (a >= b)  stop("lower bound must be smaller than upper bound")
    
    
    if(is.null(lmax)){
      Llines <- as.psp(L)
      linemask <- as.mask.psp(Llines)
      lineimage <- as.im(linemask)
      xx <- raster.x(linemask)
      yy <- raster.y(linemask)
      mm <- linemask$m
      xx <- as.vector(xx[mm])
      yy <- as.vector(yy[mm])
      pixelcentres <- ppp(xx, yy, window=as.rectangle(linemask), check=check)
      pixelcentres <- unique.ppp(pixelcentres)
      pixdf <- data.frame(xc=xx, yc=yy)
      p2s <- project2segment(pixelcentres, Llines)
      projloc <- as.data.frame(p2s$Xproj)
      projmap <- as.data.frame(p2s[c("mapXY", "tp")])
      projdata <- cbind(pixdf, projloc, projmap)
      gridx <- p2s$Xproj$x
      gridy <- p2s$Xproj$y
      df <- data.frame(gridx,gridy)
      df <- df[!duplicated(df), ]
      grid <- lpp(df,L)
      grid <- unique(grid)
      t0 <- runif(npoints(grid),a,b)
      
      lmax=max(lambda(grid$data$x,grid$data$y,t0))
    }
    
    
    mean <- lmax*volume(L)*(b-a)
    
    n <- rpois(1,mean)
    
    unipoint <- runiflpp(n,L)
    
    hlpp <- cbind(unipoint$data$x,unipoint$data$y,runif(n,a,b))
    
    prob <- lambda(hlpp[,1],hlpp[,2],hlpp[,3])/lmax
    
    if(check) {
      if(any(prob < 0))
        warning("Negative values of lambda obtained")
      if(any(prob > 1))
        warning("lmax is not an upper bound for lambda")
    }
    
    u <- runif(length(hlpp[,1]))
    
    retain <-  (u <= prob)
    
    stlpp <- hlpp[retain,]
    stlpp <- data.frame(x=stlpp[,1],y=stlpp[,2],t=stlpp[,3])
  }
  else if(any(class(lambda)=="stlppint")){
    
    Y <- attr(lambda,"stlpp")
    a <- Y$time[1]
    b <- Y$time[2]
    L <- Y$domain
    
    lmax <- max(unlist(lapply(lambda, max)))
    
    mean <- lmax*volume(L)*(b-a)
    
    n <- rpois(1,mean)
    
    unipoint <- runiflpp(n,L)
    
    hlpp <- stlpp(unipoint,T=runif(npoints(unipoint),a,b),L=L)
    
    prob <- lambda[hlpp]/lmax
    
    if(check) {
      if(any(prob < 0))
        warning("Negative values of lambda obtained")
      if(any(prob > 1))
        warning("lmax is not an upper bound for lambda")
    }
    
    u <- runif(npoints(hlpp))
    
    retain <-  (u <= prob)
    
    stlpp <- hlpp[retain]
  }
  
  out <- ppx(data=stlpp,domain = L,coord.type = c("s","s","t"))
  class(out) <- c("stlpp","ppx")
  out$time <- c(a,b)
  return(out)
}
