\name{print.stlppint}
\alias{print.stlppint}
\title{Methods for space-time point patterns on a linear network.}
\usage{
print.stlppint(x)
}
\description{
Methods to print an object of class \code{\link{stlppint}}.
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{as.stlpp}}, \code{\link{density.stlpp}}
}
\examples{
X = rpoistlpp(.2,a=0,b=5,L=easynet)
density(X)
}
