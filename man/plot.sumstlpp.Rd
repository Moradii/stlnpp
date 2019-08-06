\name{plot.sumstlpp}
\alias{plot.sumstlpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods for space-time point patterns on a linear network.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot(X,style=c("level","contour","perspective"), theta = 35, phi = 10,
facets = FALSE, ticktype = "detailed", resfac = 5,
xlab="r = distance",ylab="t = time",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
an object of class sumstlpp.
}
\item{style}{type of plot to be exhibited.}
\item{...}{arguments passed to \code{\link{levelplot}},\code{\link{contour}} or \code{\link{persp3D}} .}

}
\details{
This function plots K-function and pair correlation function for spatio-temporal point patterns on linear networks using different styles.
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{levelplot}}, \link[graphics]{contour}, \code{\link{persp3D}}
}
