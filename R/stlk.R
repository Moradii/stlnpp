#' @export
STLK <- function(X,r=NULL,t=NULL,nxy=10){
  if (!inherits(X, "stlpp")) stop("X should be from class stlpp")

  Y <- as.lpp.stlpp(X)
  l <- domain(Y)
  tleng <- summary(l)$totlength
  n <- npoints(X)
  a <- X$time[1]
  b <- X$time[2]
  trange <- b-a
  timev <- X$data$t

  norm <- (tleng*trange)/(n*(n-1))

  sdist <- pairdist.lpp(Y)
  tdist <- as.matrix(dist(timev))

  ##
  toler <- default.linnet.tolerance(l)
  ml <- matrix(1, n, n)
  for(j in 1:n) {
    ml[ -j, j] <- countends(l, Y[-j], sdist[-j,j], toler=toler)
  }

  mtplus <- matrix(timev,n,n,byrow = T)+tdist
  mtminus <- matrix(timev,n,n,byrow = T)-tdist
  mtedge <- (mtplus<=b) + (mtminus>=a)
  diag(mtedge) <- 1

  edgetl <- mtedge*ml

  # maxs <- 0.98*boundingradius.linnet(l)

  maxs <- 0.7*max(sdist[!is.infinite(sdist)])
  maxt <- 0.7*(trange/2)

  if(is.null(r)) r <- seq((maxs/nxy),maxs,by=(maxs-(maxs/nxy))/(nxy-1))
  if(is.null(t)) t <- seq((maxt/nxy),maxt,by=(maxt-(maxt/nxy))/(nxy-1))

  K <- matrix(NA,nrow = nxy,ncol = nxy)
  no <- sdist == 0 & tdist == 0 | sdist==Inf

  for (i in 1:length(r)) {

    for (j in 1:length(t)) {
      out <- (sdist<=r[i])*(tdist<=t[j])
      kout <- out[!no]/edgetl[!no]

      # kout <- out/edgetl
      K[i,j] <- sum(kout[!is.na(kout) & !is.infinite(kout)])
    }
  }

  K <- K*norm

  ##
  pixcor <- expand.grid(r,t)
  Kout <- list(Kest=K,Ktheo=matrix(pixcor[,1]*pixcor[,2],ncol = nxy),r=r,t=t)
  class(Kout) <- c("sumstlpp")
  attr(Kout,"nxy") <- nxy
  return(Kout)
}
