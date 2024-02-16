#' Random thinning
#'
#' This function applies independent random thinning to a spatio-temporal point pattern on a linear network.
#'
#' @usage \method{rthin}{stlpp}(X, P = P, nsim = 1)
#'
#' @param X a spatio-temporal point pattern of class \code{\link{stlpp}}
#' @param P retention probability
#' @param nsim number of simulated realisations to be generated
#' 
#' @seealso \code{\link{stlpp}}, \code{\link{rthin}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com> 
#' 
#' @returns 
#' An object of the same kind as X if nsim=1, or a list of such objects if nsim > 1.
#'
#' @details 
#' See \code{\link{rthin}}.
#' 
#' @references Moradi, M., & Mateu, J. (2020). First-and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics, 29(3), 432-443.
#' 
#' 
#' @examples  
#' data(Medellin)
#' rthin(Medellin,P=.5)
#'
#' @export
rthin.stlpp <- function(X,P=P,nsim=1){

  if(!any(class(X)=="stlpp")) stop("class(X) must be stlpp")
  if (P<0 | P>1) stop("P should be a numeric value between 0 and 1")

  if (nsim==1){
    u <- runif(npoints(X),0,1)
    OK <- (u<P)
    return(sub.stlpp(X,OK))

  }else {
    out <- list()
    for (i in 1:nsim) {
      u <- runif(npoints(X),0,1)
      OK <- (u<P)
      out[[i]] <- sub.stlpp(X,OK)
    }
    class(out) <- "list"
    return(out)
  }
}

sub.stlpp <- function(x,i){
  d <- as.data.frame(x$data[i,])
  as.stlpp(d$x,d$y,d$t,L=x$domain)
}