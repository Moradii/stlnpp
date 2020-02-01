\name{rthin.stlpp}
\alias{rthin.stlpp}
\title{
Random thinning
}
\description{
Applies independent random thinning to a spatio-temporal point pattern on a linear network.}
\usage{
\method{rthin}{stlpp}(X, P = P, nsim = 1)
}
\arguments{
  \item{X}{
a spatio-temporal point pattern of class \code{\link{stlpp}}}
  \item{P}{
retention probability
}
  \item{nsim}{
number of simulated realisations to be generated}
}
\details{
See \code{\link{rthin}}.}
\value{
An object of the same kind as X if nsim=1, or a list of such objects if nsim > 1.
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{stlpp}}, \code{\link{rthin}}
}
\examples{
data(Medellin)
rthin(Medellin,P=.5)
}
