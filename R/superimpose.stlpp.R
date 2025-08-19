#' Methods for spatio-temporal point patterns on a linear network
#'
#' Superimpose several spatio-temporal point patterns on a linear network.
#'
#' @usage \method{superimpose}{stlpp}(...)
#'
#' @param ... any number of arguments, each of which represents a spatio-temporal point pattern on a linear network
#' 
#' @seealso \code{\link{stlpp}},
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com>
#' 
#' @returns An object of class \code{\link{stlpp}}.
#' 
#' @details This function is used to superimpose several spatio-temporal point patterns on a linear network, each of which of class \code{\link{stlpp}}. They should be one same linear network and have the same time interval.
#' 
#' @examples  
#' data(easynet)
#' X <- rpoistlpp(2, L = easynet, a = 0, b = 2)
#' Y <- rpoistlpp(2, L = easynet, a = 0, b = 2)
#' superimpose(X,Y)
#' 
#' @export
superimpose.stlpp <- function(...) {
  objs <- list(...)
  
  # Check all inputs are stlpp
  if(!all(sapply(objs, function(x) inherits(x, "stlpp")))) {
    stop("All inputs must be stlpp objects")
  }
  
  # Combine the data frames
  combined_data <- as.data.frame(do.call(rbind, lapply(objs, function(x) x$data)))
  
  # Keep the domain of the first object (or you could decide how to handle differing domains)
  domain <- objs[[1]]$domain
  time <- objs[[1]]$time
  
  # Create a new stlpp object
  out <- as.stlpp(x=combined_data, L = domain, T = time)
  
  return(out)
}