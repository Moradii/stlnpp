\name{stlpp}
\alias{stlpp}
\title{
Create spatio-temporal point pattern on linear network
}
\description{
Create spatio-temporal point pattern on linear network.
}
\usage{
stlpp(X, L, T, ...)
}
\arguments{
  \item{X}{
Locations of the points. A matrix or data frame of coordinates, or a point pattern object (of class "ppp") or other data acceptable to \code{\link{as.ppp}} or  \code{\link{lpp}}}
  \item{L}{
linear network (object of class "linnet"). No need to be given if X is an object of \code{\link{lpp}}
}
  \item{T}{
time vector
}
  \item{\dots}{
ignored
}
}
\details{
This function creates an object of class stlpp. For details about X see \code{\link{lpp}}. T represents the time occurrence of data points.

}
\value{
an object of class stlpp.
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}
\seealso{
\code{\link{as.stlpp}}, \code{\link{lpp}}
}
\examples{
X <- rpoislpp(1,easynet)
t <- runif(npoints(X))
stlpp(X,T=t,L=easynet)
}
