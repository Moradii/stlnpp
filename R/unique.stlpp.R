#' Extract unique points from a spatio-temporal point pattern on a linear network
#'
#' This function extracts unique points from a spatio-temporal point pattern on a linear network.
#'
#' @usage \method{unique}{stlpp}(x,...)
#'
#' @param x a spatio-temporal point pattern of class \code{\link{stlpp}}
#' @param ... arguments for \code{\link{unique}}
#' 
#' @seealso \code{\link{unique}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com> 
#' 
#' @returns 
#' A spatio-temporal point pattern on a linear network with no duplicated point.
#'
#' @details 
#' This function extracts unique points from a spatio-temporal point pattern on a linear network.
#' 
#' @references Moradi, M., & Mateu, J. (2020). First-and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics, 29(3), 432-443.
#' 
#' 
#' @examples  
#' X <-  rpoistlpp(0.1,0,5,L=easynet)
#' df <- as.data.frame(X)
#' df_dup <- df[sample(nrow(df), 20,replace = TRUE), ]
#' Y <- as.stlpp(df_dup,L=easynet)
#' npoints(Y)
#' npoints(unique(Y))
#'
#' 
#' 
#' @export
unique.stlpp <- function(x,...){
  if (!inherits(x, "stlpp")) stop("X should be from class stlpp")

  Y <- unique(as.data.frame(x$data),...)
  X <- as.stlpp(Y$x,Y$y,Y$t,L=domain(x))
  return(X)
}
