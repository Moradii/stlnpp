\name{plot.stlpp}
\alias{plot.stlpp}
\title{Methods for space-time point patterns on a linear network.}
\usage{
plot(X,...)
}
\description{
plot objecs of class \code{\link{stlpp}}.}

\arguments{
  \item{X}{an object of classes \code{\link{stlpp}}.}
  \item{...}{arguments passed to \code{\link{plot}}}
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{as.stlpp}}
}
\examples{
X = rpoistlpp(.2,a=0,b=5,L=easynet)
plot(X)
}
