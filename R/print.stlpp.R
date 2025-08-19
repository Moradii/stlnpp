
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
print.stlppint <- function(x, ...) {
  if (!inherits(x, "stlppint")) stop("x must be of class 'stlppint'")
  
  stlpp_obj <- attr(x, "stlpp")  # object containing $time
  timex <- stlpp_obj$time
  cat("Object of class 'stlppint'\n")
  
  # number of images / time points
  n_time <- length(x)
  if(!is.null(attr(x, "tgrid"))){
    t <- attr(x, "tgrid")
  }else{
    t <- attr(x, "tempden")$x
  }
  # t <- attr(x, "time")
  cat("Number of time points/images:", n_time, "\n")
  cat("Spatial domain:", "\n")
  print(stlpp_obj$domain)
  cat("Tempral domain:", "\n")
  cat("[","",timex[1],",",timex[2], "]", "\n", sep="")
  
  # compute average intensity per spatial image
  avg_spatial <- sapply(x, function(im) mean(as.vector(as.im(im)), na.rm = TRUE))
  
  if (!is.null(avg_spatial)) {
    max_idx <- which.max(avg_spatial)
    min_idx <- which.min(avg_spatial)
    
    cat(sprintf("Highest average spatio-temporal intensity at time %.2f (%.3f)\n",
                t[max_idx], avg_spatial[max_idx]))
    cat(sprintf("Lowest average spatio-temporal intensity at time %.2f (%.3f)\n",
                t[min_idx], avg_spatial[min_idx]))
  }
  
  # temporal intensity if available
  if (!is.null(attr(x, "tempden"))) {
    
    # get the time of highest/lowest intensity using stlpp attribute
    
    stlpp_time <- stlpp_obj$time
    
    d <- attr(x, "tempden")
    tgrid <- d$x
    int <- npoints(stlpp_obj) * d$y
    
    
    
    max_idx <- which.max(int)
    min_idx <- which.min(int)
    
    cat(sprintf("Highest temporal intensity during time period %.2f (%.3f)\n",
                tgrid[max_idx], int[max_idx]))
    cat(sprintf("Lowest temporal intensity during time period %.2f (%.3f)\n",
                tgrid[min_idx], int[min_idx]))
    
  } else if (!is.null(attr(x, "tint"))) {
    d <- attr(x, "tint")
    tgrid <- attr(x, "tgrid")
    
    stlpp_obj <- attr(x, "stlpp")
    stlpp_time <- stlpp_obj$time
    
    max_idx <- which.max(d)
    min_idx <- which.min(d)
    
    cat(sprintf("Highest temporal intensity during time period %.2f (%.3f)\n",
                tgrid[max_idx], d[max_idx]))
    cat(sprintf("Lowest temporal intensity during time period %.2f (%.3f)\n",
                tgrid[min_idx], d[min_idx]))
  } else if(!is.null(attr(x, "tgrid"))){
    d <- as.numeric(as.tppint.stlppint(x))
    tgrid <- attr(x, "tgrid")
    stlpp_time <- stlpp_obj$time
    
    max_idx <- which.max(d)
    min_idx <- which.min(d)
    
    cat(sprintf("Highest temporal intensity during time period %.2f (%.3f)\n",
                tgrid[max_idx], d[max_idx]))
    cat(sprintf("Lowest temporal intensity during time period %.2f (%.3f)\n",
                tgrid[min_idx], d[min_idx]))
  }
  
  invisible(x)
} 

# print.stlppint <- function(x,...) {
#   if(any(class(x)=="list")) print("List of images on linear network")
#   if(any(class(x)=="numeric")) print(as.vector(x),...)
# }

#' @export
print.sumstlpp <- function(x,...)
{
  cat("summary statistics for class stlpp \n");
  if(any(names(x)=="Kest")){cat(paste0("homogeneous K-function"),"\n")}
  else if (any(names(x)=="gest")){cat(paste0("homogeneous pair correlation function"),"\n")}
  else if (any(names(x)=="Kinhom")){cat(paste0("inhomogeneous K-function"),"\n")}
  else if (any(names(x)=="ginhom")){cat(paste0("inhomogeneous pair correlation function"),"\n")}
  else if (any(names(x)=="Hinhom")){cat(paste0("inhomogeneous nearest neighbour function"),"\n")}
  else if (any(names(x)=="Finhom")){cat(paste0("inhomogeneous empty space function"),"\n")}
  else if (any(names(x)=="Jinhom")){cat(paste0("inhomogeneous J-function"),"\n")}
  else{cat(paste0("strange"),"\n")}
}
