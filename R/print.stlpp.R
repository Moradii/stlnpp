
#### Function "print.stlpp" prints an object of class stlpp and
#### "print.stlppint" prints an objetc of class stlppint (e.g. estimated intensity of an
### object of class stlpp).
#' @export
print.stlpp <- function(x,...)
{
  cat("Spatio-temporal point pattern on a linear network \n");
 if(npoints(x)>1){cat(paste0(npoints(x)," ", "points"),"\n")}
  else{cat(paste0(npoints(x)," ", "point"),"\n")};
  print(x$domain,...);cat(paste0("Time period: [",x$time[1],", ", x$time[2],"]"),"\n")
  }

#' @export
print.stlppint <- function(x,...) {
  if(any(class(x)=="list")) print("List of images on linear network")
  if(any(class(x)=="numeric")) print(as.vector(x),...)
}

#' @export
print.sumstlpp <- function(x,...)
{
  cat("summary statistics for class stlpp \n");
  if(any(names(x)=="Kest")){cat(paste0("homogeneous K-function"),"\n")}
  else if (any(names(x)=="gest")){cat(paste0("homogeneous pair correlation function"),"\n")}
  else if (any(names(x)=="Kinhom")){cat(paste0("inhomogeneous K-function"),"\n")}
  else if (any(names(x)=="ginhom")){cat(paste0("inhomogeneous pair correlation function"),"\n")}
  else{cat(paste0("strange"),"\n")}
}
