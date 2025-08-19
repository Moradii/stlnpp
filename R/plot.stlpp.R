#' @import graphics
#' @export
plot.stlpp <- function(x, onlypoints = FALSE, onlytime = FALSE,
                       args.xy = list(), args.time = list(),...) {
  
  if (!inherits(x, "stlpp")) stop("x should be an object of stlpp")
  oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar))
  
  if (onlypoints) {
    if (is.null(args.xy$main)) args.xy$main <- "xy-locations on linear network"
    do.call(plot, c(list(as.lpp.stlpp(x)), args.xy))
    
  } else if (onlytime) {
    xx <- sort(as.data.frame(x$data[, 3])[, 1], index.return = TRUE)
    x1 <- as.data.frame(x$data[xx$ix, ])
    
    if (is.null(args.time$xlab)) args.time$xlab <- "time"
    if (is.null(args.time$ylab)) args.time$ylab <- ""
    if (is.null(args.time$main)) args.time$main <- "cumulative number"
    
    do.call(plot, c(list(x1[, 3], cumsum(x1[, 3]), type = "l",
                         las = 1, xlim = c(min(x1[, 3]), max(x1[, 3]))),
                    args.time))
    
  } else {
    par(mfrow = c(1, 2), pty = "s")
    
    if (is.null(args.xy$main)) args.xy$main <- "xy-locations on linear network"
    do.call(plot, c(list(as.lpp.stlpp(x)), args.xy))
    
    xx <- sort(as.data.frame(x$data[, 3])[, 1], index.return = TRUE)
    x1 <- as.data.frame(x$data[xx$ix, ])
    
    if (is.null(args.time$xlab)) args.time$xlab <- "time"
    if (is.null(args.time$ylab)) args.time$ylab <- ""
    if (is.null(args.time$main)) args.time$main <- "cumulative number"
    
    do.call(plot, c(list(x1[, 3], cumsum(x1[, 3]), type = "l",
                         las = 1, xlim = c(min(x1[, 3]), max(x1[, 3]))),
                    args.time))
  }
}

#' 
# plot.tpp=function(x, onlypoints=FALSE, onlytime=FALSE, xlab=xlab,...){
#   
#   if (inherits(x, "tpp") == FALSE) stop("x should be an object of tpp")
#   oldpar <- par(no.readonly = TRUE)
#   on.exit(par(oldpar)) 
#   
#   t <- x$data$t
#   
#   if(onlypoints){
#     y_values <- rep(0, npoints(x))
#     plot(t, y_values, yaxt = "n", bty = "n", ...)
#     
#   }else if(onlytime){
#     xx = sort(t, index.return = TRUE)
#     t = t[xx$ix]
#     plot(t, cumsum(t), type = "l", las = 1,
#          xlim=c(min(t),max(t)), ...)
#   }
#   else{
#     par(mfrow = c(1, 2), pty = "s")
#     y_values <- rep(0, npoints(x))
#     plot(t, y_values, yaxt = "n", bty = "n", ylab="", ...)
#     
#     xx = sort(t, index.return = TRUE)
#     x1 = t[xx$ix]
#     if (missing(xlab)) xlab <- "time"
#     plot(x1, cumsum(x1), type = "l", xlab=xlab,
#          ylab = "", main = "cumulative number", las = 1
#          ,xlim=c(min(x1),max(x1)),...)
#   }
# }


#' @export
plot.stlppint <- function(x,
               style = NULL,
               xlab = NULL, ylab = NULL, main = NULL,
               xlim = NULL, ylim = NULL,
               bar = TRUE,
               args.spatial = list(),
               args.temporal = list(),
               all.images = FALSE,
               per.page = 9,
               ...) {
  
  if (!inherits(x, "stlppint"))
    stop("x must be of class 'stlppint'")
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # ------------------------------------------------------------
  # CASE 1: plot all images over time (spatial only)
  # ------------------------------------------------------------
  if (all.images) {
    nImgs <- length(x)
    nPages <- ceiling(nImgs / per.page)
    
    # global zlim across all linims
    zrange <- range(unlist(lapply(x, range)), na.rm = TRUE)
    
    for (p in seq_len(nPages)) {
      subset <- ((p - 1) * per.page + 1):min(p * per.page, nImgs)
      n <- length(subset)
      nrow <- ceiling(sqrt(n))
      ncol <- ceiling(n / nrow)
      par(mfrow = c(nrow, ncol))
      
      for (i in subset) {
        args.spat <- utils::modifyList(list(
          x = x[[i]],
          main = if (is.null(main)) paste("time =", i) else main,
          zlim = zrange   # <--- fix: same color scale for all
        ), args.spatial)
        if (!is.null(style)) args.spat$style <- style
        do.call(plot, args.spat)
      }
      
      if (p < nPages) {
        readline("Press [Enter] to continue to next page...")
      }
    }
    return(invisible(NULL))
  }
  
  # ------------------------------------------------------------
  # CASE 2: default mode (spatial + temporal summary)
  # ------------------------------------------------------------
  t <- attr(x, "time")
  par(mfrow = c(1, 2))
  
  if (is.null(attr(x, "tempden")) & is.null(attr(x, "tint"))) {
    onL <- as.linim.stlppint(x)
    onT <- as.tppint.stlppint(x)
    tgrid <- attr(onT, "tgrid")
    
    if (is.null(xlim)) xlim <- range(tgrid)
    if (is.null(xlab)) xlab <- "time"
    if (is.null(ylim)) ylim <- c(0, max(onT, table(round(t))))
    if (is.null(ylab)) ylab <- expression(hat(lambda)[time])
    
    args.temp <- utils::modifyList(list(
      x = tgrid, y = onT,
      type = "l", ylim = ylim, xlab = xlab, ylab = "", main = "", xlim = xlim, ...
    ), args.temporal)
    do.call(plot, args.temp)
    if (bar) points(table(round(t)))
    title(ylab = ylab, line = 2, cex = 1.2)
    
    args.spat <- utils::modifyList(list(
      x = onL, main = "",
      zlim = range(onL, na.rm = TRUE)  # consistent scaling
    ), args.spatial)
    if (!is.null(style)) args.spat$style <- style
    do.call(plot, args.spat)
    
  } else if (!is.null(attr(x, "tempden"))) {
    d <- attr(x, "tempden")
    int <- length(t) * d$y
    
    if (is.null(xlim)) xlim <- range(d$x)
    if (is.null(ylim)) ylim <- c(0, max(int, table(round(t))))
    if (is.null(xlab)) xlab <- "time"
    if (is.null(ylab)) ylab <- expression(hat(lambda)[time])
    OK <- d$x >= range(xlim)[1] & d$x <= range(xlim)[2]
    
    args.temp <- utils::modifyList(list(
      x = d$x[OK], y = int[OK],
      type = "l", ylim = ylim, xlab = xlab, ylab = "", main = "", xlim = xlim, ...
    ), args.temporal)
    do.call(plot, args.temp)
    if (bar) points(table(round(t)))
    title(ylab = ylab, line = 2, cex = 1.2)
    
    args.spat <- utils::modifyList(list(
      x = attr(x, "netint"), main = "",
      zlim = range(attr(x, "netint"), na.rm = TRUE)
    ), args.spatial)
    if (!is.null(style)) args.spat$style <- style
    do.call(plot, args.spat)
    
  } else {
    d <- attr(x, "tint")
    tgrid <- attr(x, "tgrid")
    
    if (is.null(xlim)) xlim <- range(tgrid)
    if (is.null(xlab)) xlab <- "time"
    if (is.null(ylim)) ylim <- c(0, max(d, table(round(t))))
    if (is.null(ylab)) ylab <- expression(hat(lambda)[time])
    
    args.temp <- utils::modifyList(list(
      x = tgrid, y = d,
      type = "l", ylim = ylim, xlab = xlab, ylab = "", main = "", xlim = xlim, ...
    ), args.temporal)
    do.call(plot, args.temp)
    if (bar) points(table(round(t)))
    title(ylab = ylab, line = 2, cex = 1.2)
    
    args.spat <- utils::modifyList(list(
      x = attr(x, "lint"), main = "",
      zlim = range(attr(x, "lint"), na.rm = TRUE)
    ), args.spatial)
    if (!is.null(style)) args.spat$style <- style
    do.call(plot, args.spat)
  }
  
  invisible(NULL)
}

# plot.stlppint <- function(x,style=style,xlab=xlab,xlim=xlim,ylim=ylim,bar=TRUE,...){
# 
#   if (inherits(x, "stlppint") == FALSE) stop(" x must be from class stlppint")
# 
#   oldpar <- par(no.readonly = TRUE)
#   on.exit(par(oldpar))
# 
#   t <- attr(x,"time")
#   par(mfrow=c(1,2))
# 
#   if(is.null(attr(x,"tempden")) & is.null(attr(x,"tint"))){
#     onL <- as.linim.stlppint(x)
#     onT <- as.tppint.stlppint(x)
# 
#     tgrid <- attr(onT,"tgrid")
# 
#     if (missing(xlim)) xlim <- range(tgrid)
#     if (missing(xlab)) xlab <- "time"
#     if (missing(ylim)) ylim = c(0,max(onT,table(round(t))))
# 
#     plot(tgrid,onT,ylab="",main="",
#          type="l",ylim = ylim,
#          xlab=xlab,xlim = xlim,...)
#     if(bar) points(table(round(t)))
#     title(ylab=expression(hat(lambda)[time]), line=2,cex=3,...)
# 
#     if (missing(style)) {plot(onL,main="",...)}
#     else {plot(onL,main="",style=style,...)}
#   }
#   else if(!is.null(attr(x,"tempden"))){
#     d <- attr(x,"tempden")
#     int <- length(t)*d$y
# 
#     if (missing(xlim)) xlim <- range(d$x)
#     if (missing(ylim)) ylim = c(0,max(int,table(round(t))))
#     OK <- d$x>=range(xlim)[1] & d$x<=range(xlim)[2]
#     if (missing(xlab)) xlab <- "time"
#     plot(d$x[OK],int[OK],
#          ylab="",main="",type="l",ylim = ylim,xlab=xlab,xlim = xlim,...)
#     if(bar) points(table(round(t)))
#     title(ylab=expression(hat(lambda)[time]), line=2,cex=3,...)
# 
#     if (missing(style)) {plot(attr(x,"netint"),main="",...)}
#     else {plot(attr(x,"netint"),main="",style=style,...)}
#   }
#   else{
#     d <- attr(x,"tint")
#     tgrid <- attr(x,"tgrid")
# 
#     if (missing(xlim)) xlim <- range(tgrid)
#     if (missing(xlab)) xlab <- "time"
#     if (missing(ylim)) ylim = c(0,max(d,table(round(t))))
# 
#     plot(tgrid,d,ylab="",main="",
#          type="l",ylim = ylim,
#          xlab=xlab,xlim = xlim,...)
#     if(bar) points(table(round(t)))
#     title(ylab=expression(hat(lambda)[time]), line=2,cex=3,...)
# 
#     if (missing(style)) {plot(attr(x,"lint"),main="",...)}
#     else {plot(attr(x,"lint"),main="",style=style,...)}
#   }
# 
# }

#' @export
plot.sumstlpp <- function(x,
                          style = c( "perspective",  "contour","level"),
                          theta = 35, phi = 10, facets = FALSE,
                          ticktype = "detailed", resfac = 5,
                          xlab = NULL, ylab = NULL, zlab = NULL, main = NULL,
                          args.level = list(),
                          args.contour = list(),
                          args.perspective = list(),...) {
  
  style <- match.arg(style)
  
  if (!inherits(x, "sumstlpp")) stop("x is not an object of class sumstlpp")
  
  # Required packages
  if (!requireNamespace("lattice", quietly = TRUE)) stop("lattice required")
  if (!requireNamespace("graphics", quietly = TRUE)) stop("graphics required")
  if (!requireNamespace("plot3D", quietly = TRUE)) stop("plot3D required")
  
  # Identify which statistic is present
  stats <- c("Kest", "gest", "Kinhom", "ginhom", "Hinhom", "Finhom", "Jinhom")
  stat_name <- stats[stats %in% names(x)][1]
  if (is.na(stat_name)) stop("No recognized summary statistic found in object")
  
  # Extract values
  z <- x[[stat_name]]
  ztheo_name <- paste0(tolower(sub(".*", "theo", stat_name)))
  ztheo <- if (ztheo_name %in% names(x)) x[[ztheo_name]] else NULL
  
  # Default labels if not provided
  if (is.null(xlab)) xlab <- "r = distance"
  if (is.null(ylab)) ylab <- "t = time"
  if (is.null(zlab)) zlab <- ""
  if (is.null(main)) {
    main_labels <- list(
      Kest = expression(italic({hat(K)[L]^{ST}}(r,t))),
      gest = expression(italic({hat(g)[L]^{ST}}(r,t))),
      Kinhom = expression(italic({hat(K)[LI]^{ST}}(r,t))),
      ginhom = expression(italic({hat(g)[LI]^{ST}}(r,t))),
      Hinhom = expression(italic({hat(H)[LI]^{ST}}(r,t))),
      Finhom = expression(italic({hat(F)[LI]^{ST}}(r,t))),
      Jinhom = expression(italic({hat(J)[LI]^{ST}}(r,t)))
    )
    main <- main_labels[[stat_name]]
  }
  
  # Automatic zlim
  zlim <- if (!is.null(ztheo)) range(z, ztheo, na.rm = TRUE) else range(z, na.rm = TRUE)
  
  # Plot depending on style
  if (style == "perspective") {
    args <- utils::modifyList(list(x = x$r, y = x$t, z = z,
                                   xlab = xlab,
                                   ylab = ylab,
                                   zlab = zlab,
                                   main = main,
                                   zlim = zlim,
                                   theta = theta,
                                   phi = phi,
                                   facets = facets,
                                   ticktype = ticktype,
                                   resfac = resfac), args.perspective)
    do.call(plot3D::persp3D, args)
  }
  else if (style == "contour") {
    args <- utils::modifyList(list(x = x$r, y = x$t, z = z,
                                   xlab = xlab,
                                   ylab = ylab,
                                   main = main), args.contour)
    do.call(graphics::contour, args)
    
  }
  else if (style == "level") {
    args <- utils::modifyList(list(x = z,
                            row.values = x$r,
                            column.values = x$t,
                            xlab = xlab,
                            ylab = ylab,
                            main = main), args.level)
    do.call(lattice::levelplot, args)
    
  } 
}

