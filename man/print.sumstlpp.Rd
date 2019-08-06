\name{print.sumstlpp}
\alias{print.sumstlpp}
\title{Methods for space-time point patterns on a linear network.}
\usage{
print.sumstlpp(x)
}
\description{
Methods to print an object of class \code{\link{stlpp}}.
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{STLg}}, \code{\link{STLK}}, \code{\link{STLginhom}}, \code{\link{STLKinhom}}
}
\examples{
X = rpoistlpp(.2,a=0,b=5,L=easynet)
STLK(X)
STLg(X)
}
