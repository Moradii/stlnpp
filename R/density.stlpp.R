#' @export
density.stlpp <- function(X,lbw,tbw,...){

  if (!inherits(X, "stlpp")) stop("X should an object of class stlpp")

    ox <- X$data$x
    oy <- X$data$y
    ot <- X$data$t

    L <-  X$domain

    n <- npoints(X) # Emerge number of points

    stint <- 0 # define the vacant vector to save the densities values per points

    if (missing(tbw)) {d <- density(ot)}
    else {
      d <- density(ot,bw=tbw)
    }

    Tint <- d$y[findInterval(ot, d$x)] * length(ot)
    ############################################## space intensity

    pX <- as.stlpp.lpp(X)
    if (missing(lbw)) lbw <- bw.scott.iso(pX)
    Sint <- density.lpp(pX,sigma = lbw,at="points",distance="euclidean",...)

    stint <-  Sint*Tint/npoints(pX)
    
  out <- stint
  ldens <- density.lpp(pX,sigma = lbw,distance="euclidean",...)

  names(lbw) <- NULL
  attr(out,"tempden") <- d
  attr(out,"netint") <- ldens
  attr(out,"time") <- ot
  attr(out,"bw") <- c("sigma_l"=lbw,"sigma_t"=d$bw)
  class(out) <- c("stlppint")
  return(out)

}
