
#### Function "print.stlpp" prints an object of class stlpp and
#### "print.stlppint" prints an objetc of class stlppint (e.g. estimated intensity of an
### object of class stlpp).
#' @export
print.stlpp <- function(X)
{
  cat("Spatio-temporal point pattern on linear network \n");
 if(npoints(X)>1){cat(paste0(npoints(X)," ", "points"),"\n")}
  else{cat(paste0(npoints(X)," ", "point"),"\n")};
  print(X$domain);cat(paste0("Time period: [",X$time[1],", ", X$time[2],"]"),"\n")
  }

#' @export
print.stlppint <- function(X) {
  if(any(class(X)=="list")) print("List of images on linear network")
  if(any(class(X)=="numeric")) print(as.vector(X))
}

#' @export
print.sumstlpp <- function(X)
{
  cat("summary statistics for class stlpp \n");
  if(any(names(X)=="Kest")){cat(paste0("homogeneous K-function"),"\n")}
  else if (any(names(X)=="gest")){cat(paste0("homogeneous pair correlation function"),"\n")}
  else if (any(names(X)=="Kinhom")){cat(paste0("inhomogeneous K-function"),"\n")}
  else if (any(names(X)=="ginhom")){cat(paste0("inhomogeneous pair correlation function"),"\n")}
  else{cat(paste0("strange"),"\n")}
}
