\name{plot.stlppint}
\alias{plot.stlppint}
\title{Methods for space-time point patterns on a linear network.}
\usage{
plot.stlppint(X)
}

\description{
plot an estimated intensity function of a space-time point patterns on a linear network.
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{as.stlpp}}, \code{\link{density.stlpp}}
}
\examples{
X = rpoistlpp(.2,a=0,b=5,L=easynet)
d = density(X)
plot(d)
}
