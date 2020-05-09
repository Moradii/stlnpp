#' @export
STLginhom <- function(X,lambda,normalize=FALSE,r=NULL,t=NULL,nxy=10){


  if (!inherits(X, "stlpp")) stop("X should be from class stlpp")


  Y <- as.lpp.stlpp(X)
  l <- domain(Y)
  tleng <- summary(l)$totlength
  n <- npoints(Y)
  a <- X$time[1]
  b <- X$time[2]
  trange <- b-a
  timev <- X$data$t

  sdist <- pairdist(Y)
  tdist <- as.matrix(dist(timev))

  toler <- default.linnet.tolerance(l)
  ml <- matrix(1, n, n)
  for (j in 1:n) {
    ml[-j, j] <- countends(l, Y[-j], sdist[-j, j], toler = toler)
  }

  mtplus <- matrix(timev, n, n, byrow = T) + tdist
  mtminus <- matrix(timev, n, n, byrow = T) - tdist
  mtedge <- (mtplus <= b) + (mtminus >= a)
  diag(mtedge) <- 1

  lamden <- outer(lambda,lambda,FUN = "*")

  edgetl <- mtedge * ml * lamden

  maxs <- 0.7*max(sdist[!is.infinite(sdist)])
  maxt <- 0.7*(trange/2)

  if(is.null(r)) r <- seq((maxs/nxy),maxs,by=(maxs-(maxs/nxy))/(nxy-1))
  if(is.null(t)) t <- seq((maxt/nxy),maxt,by=(maxt-(maxt/nxy))/(nxy-1))

  g <- matrix(NA, nrow = nxy, ncol = nxy)
  no <- sdist == 0 & tdist == 0 | sdist==Inf | sdist>maxs | tdist>maxt
  bwl <- bw.nrd0(as.numeric(sdist[!no]))
  bwt <- bw.nrd0(as.numeric(tdist[!no]))

  for (i in 1:length(r)) {
    for (j in 1:length(t)) {

      outl <- dkernel(as.numeric(sdist[!no] - r[i]),
                      sd = bwl)
      outt <- dkernel(as.numeric(tdist[!no] - t[j]),
                      sd = bwt)
      g1 <- outl * outt/(edgetl[!no])
      g[i, j] <- sum(g1[!is.na(g1) & !is.infinite(g1)])
    }
  }

  if(normalize){
    revrho <- outer(1/lambda,1/lambda,FUN = "*")
    appx <- (tleng*trange)/(sum(revrho[lower.tri(revrho, diag = FALSE)])*2)
    gval <- g * appx
  }
  else {
    gval <- g/(tleng*trange)
  }

  gout <- list(ginhom = gval, gtheo = matrix(rep(1, length(t) *
                                                 length(r)), ncol = nxy), r = r, t = t)
  class(gout) <- c("sumstlpp")
  attr(gout,"nxy") <- nxy

  return(gout)
}
