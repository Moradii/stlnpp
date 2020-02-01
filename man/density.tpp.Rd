\name{density.tpp}
\alias{density.tpp}
\title{
Kernel estimation of intensity of temporal point patterns
}
\description{
Kernel estimation of intensity of temporal point patterns.
}
\usage{
\method{density}{tpp}(x,tbw,at=c("points","pixels"),...)
}
\arguments{
  \item{x}{
a temporal point pattern (of class \code{\link{tpp}})
}
\item{tbw}{time smoothing bandwidth}
\item{at}{string specifying whether to compute the intensity values at a grid of pixel locations (at="pixels") or only at the points of x (at="points"). default is to estimate the intensity at pixels}
\item{...}{arguments passed to \link[stats]{density}}
}
\value{
A vector of intensity values.
}

\references{
Moradi, M.M. and Mateu, J. (2019). First and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics. In press.
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{density}}, \code{\link{bw.nrd0}}
}

\examples{
X <- tpp(sample(c(1:24),200,replace = TRUE))
plot(density(X))
}