#' Simulating one-dimensional Poisson point patterns
#'
#' This function simulates realisations of an one-dimensional Poisson point process.
#'
#' @usage rpoistpp(lambda,a,b,check=FALSE,lmax=NULL,nsim=1)
#'
#' @param lambda intensity of the point process. it can be either a number, a function of location and time, or an object of class \code{tppint}
#' @param a lower bound of time period
#' @param b upper bound of time period
#' @param check logical value indicating whether to check that all the (x,y) points lie inside the specified time period.
#' @param lmax upper bound for the values of \code{labmda}. this is optional
#' @param nsim number of simulated patterns to generate
#' 
#' @seealso \code{\link{rpoistlpp}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com> 
#' 
#' @returns 
#' an object of class \code{\link{tpp}} if nsim=1, otherwise a list of objects of class \code{\link{tpp}}.
#'
#' @details 
#' This function generates realisations of a temporal poisson point process based on a given intensity function lambda and lower/upper bounds a and b.
#' 
#' 
#' @references Moradi, M., & Mateu, J. (2020). First-and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics, 29(3), 432-443.
#' 
#' 
#' @examples  
#' f <- function(t){0.1*exp(t)}
#' X <- rpoistpp(f,a=1,b=10)
#'
#' @export
rpoistpp <- function(lambda,a=NULL,b=NULL,check=FALSE,lmax=NULL,nsim=1){
  
  
  # if (!is.numeric(lambda) & !is.function(lambda) & !class(lambda)=="tppint")
    if (!inherits(lambda,"numeric") & !inherits(lambda,"function") & !inherits(lambda,"tppint"))
    stop(" lambda should be a number or a function")
  
  if(nsim > 1) {
    out <- list()
    for (i in 1:nsim) {
      out[[i]] <- rpoistpp(lambda,a=a,b=b,check=check,lmax=lmax,nsim=1)
    }
    return(out)
  }
  
  if (inherits(lambda,"numeric")){
    
    if (a >= b)  stop("lower bound must be smaller than upper bound")
    
    n <- rpois(1,lambda*(b-a))
    t <- runif(n,a,b)
    out <- tpp(t,a=a,b=b)
  }
  else if(inherits(lambda,"function")){
    
    if (a >= b)  stop("lower bound must be smaller than upper bound")
    
    if(is.null(lmax)){
      tgrid <- seq(a,b,length.out = 128)
      lmax <- max(lambda(tgrid))
    }
    
    mean <- lmax*(b-a)
    
    n <- rpois(1,mean)
    
    t <- runif(n,a,b)
    
    prob <- lambda(t)/lmax
    
    if(check) {
      if(any(prob < 0))
        warning("Negative values of lambda obtained")
      if(any(prob > 1))
        warning("lmax is not an upper bound for lambda")
    }
    
    u <- runif(length(t))
    
    retain <-  (u <= prob)
    
    t <- t[retain]
    out <- tpp(t,a=a,b=b)
    
  }else if(inherits(lambda,"tppint")){
    
    lmax <- max(lambda)
       
       a <- attr(lambda,"tpp")$time[1]
       
       b <- attr(lambda,"tpp")$time[2]
       
    mean <- lmax*(b-a)
       
       n <- rpois(1,mean)
       
       t <- runif(n,a,b)
       
    prob <- lambda[tpp(t)]/lmax
       
       if(check) {
         if(any(prob < 0))
           warning("Negative values of lambda obtained")
         if(any(prob > 1))
           warning("lmax is not an upper bound for lambda")
       }
       
       u <- runif(length(t))
       
       retain <-  (u <= prob)
       
       t <- t[retain]
       out <- tpp(t,a=a,b=b)
  }
  
  class(out) <- c("tpp","ppx")
  
  out$time <- c(a,b)
  
  return(out)
}