grid.stlpp <- function(L,
                       a = a,
                       b = b,
                       dimyx = 64,
                       dimt = 100,
                       ...){
  
  if (!inherits(L, "linnet")) stop("X should be from class linnet")
  
  Llines <- as.psp(L)
  
  # creating a grid in L
  linemask <- as.mask.psp(Llines,dimyx=dimyx)
  lineimage <- as.im(linemask)
  
  xx <- raster.x(linemask)
  yy <- raster.y(linemask)
  
  mm <- linemask$m
  xx <- as.vector(xx[mm])
  yy <- as.vector(yy[mm])
  
  pixelcentres <- ppp(xx, yy, window=as.rectangle(linemask), check=FALSE)
  pixdf <- data.frame(xc=xx, yc=yy)
  p2s <- project2segment(pixelcentres, Llines)
  
  
  gridx <- p2s$Xproj$x
  gridy <- p2s$Xproj$y
  grid  <- lpp(data.frame(gridx,gridy),L)
  grid  <- unique(grid,warn = F)
  
  xf <- rep(grid$data$x,each=dimt)
  yf <- rep(grid$data$y,each=dimt)
  tseq <- seq(a,b,length.out=dimt)
  tf <- rep(tseq,each=npoints(grid))
  
  return(as.stlpp(xf,yf,tf,L))
  
}