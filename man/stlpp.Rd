\name{stlpp}
\alias{stlpp}
\title{
Create spatio-temporal point pattern on linear network
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stlpp(X, L, T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Locations of the points. A matrix or data frame of coordinates, or a point pattern object (of class "ppp") or other data acceptable to \code{\link{as.ppp}} or  \code{\link{lpp}}}
  \item{L}{
linear network (object of class "linnet")
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
}
