\name{unique.stlpp}
\alias{unique.stlpp}
\title{Extract unique points from a spatio-temporal point pattern on a linear network}
\description{
Extract unique points from a spatio-temporal point pattern on a linear network.
}
\usage{
\method{unique}{stlpp}(x,...)
}
\arguments{
  \item{x}{a realisation of a spatio-temporal point processes on a linear networks. }
  \item{...}{arguments for \code{\link{unique}}.}
  }
\details{
This function calculates the inhomogeneous pair correlation function for a spatio-temporal point processes on a linear network.
}

\value{
A spatio-temporal point pattern on a linear network with no duplicated point.
}
\references{
Moradi, M.M. and Mateu, J. (2019). First and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics. In press.
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{unique}}
}

\examples{
  X <-  rpoistlpp(0.1,0,5,L=easynet)
  df <- as.data.frame(X)
  df_dup <- df[sample(nrow(df), 20,replace = TRUE), ]
  Y <- as.stlpp(df_dup,L=easynet)
  npoints(Y)
  npoints(unique(Y))
}
