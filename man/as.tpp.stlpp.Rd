\name{as.tpp.stlpp}
\alias{as.tpp.stlpp}

\title{
Methods for space-time point patterns on a linear network.
}
\description{
This function converts an object of class \code{\link{stlpp}} to  class \code{\link{tpp}}.}
\usage{
\method{as.tpp}{stlpp}(X)
}
\arguments{
\item{X}{an object of class \code{\link{stlpp}}}
}
\details{
This function projects the space-time point pattern X on L times T into the time interval T.
}
\value{
An object of class tpp.
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}


\seealso{
\code{\link{as.stlpp}}, \code{\link{lpp}}, \code{\link{as.lpp}}
}

\examples{
X <- rpoistlpp(10,1,2,easynet)
as.tpp.stlpp(X)
}