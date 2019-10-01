\name{as.stlpp.tpp}
\alias{as.stlpp.tpp}

\title{
Methods for space-time point patterns on a linear network.
}
\description{
This function converts an object of class \code{\link{stlpp}} to  class \code{\link{tpp}}.}
\usage{
as.stlpp.tpp(X)
}
\arguments{
  \item{X}{an object of class \code{\link{stlpp}}}
}
\details{
This function projects the space-time point pattern X on L times T into the time interval T.
}
\value{
an object of class tpp.
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}


\seealso{
\code{\link{as.stlpp}}, \code{\link{lpp}}, \code{\link{as.lpp}}
}
