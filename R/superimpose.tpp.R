#' Methods for temporal point patterns
#'
#' Superimpose several temporal point patterns.
#'
#' @usage \method{superimpose}{tpp}(...)
#'
#' @param ... any number of arguments, each of which represents a temporal point pattern
#' 
#' @seealso \code{\link{tpp}},
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com>
#' 
#' @returns An object of class \code{\link{tpp}}.
#' 
#' @details This function is used to superimpose several temporal point patterns, each of which of class \code{\link{tpp}}. They should have the same time interval.
#' 
#' @examples  
#' data(easynet)
#' X <- tpp(runif(10), a = 0, b = 1)
#' Y <- tpp(runif(10), a = 0, b = 1)
#' superimpose(X,Y)
#' 
#' @export
superimpose.tpp <- function(...) {
  objs <- list(...)
  
  # Check all inputs are stlpp
  if(!all(sapply(objs, function(x) inherits(x, "tpp")))) {
    stop("All inputs must be stlpp objects")
  }
  
  # Combine the data frames
  combined_data <- as.data.frame(do.call(rbind, lapply(objs, function(x) x$data)))$t
  
  # Keep the domain of the first object (or you could decide how to handle differing domains)
  time <- objs[[1]]$time
  
  # Create a new stlpp object
  out <- tpp(X = combined_data, a = time[1], b = time[2])
  
  return(out)
}