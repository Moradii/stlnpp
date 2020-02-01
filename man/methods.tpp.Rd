\name{methods.tpp}
\alias{methods.tpp}
\alias{print.tpp}
\alias{print.tppint}
\alias{plot.tpp}
\alias{plot.tppint}
\title{
Methods for temporal point patterns
}
\description{
Methods for temporal point patterns.
}
\usage{
\method{plot}{tpp}(x,...)
\method{plot}{tppint}(x,...)
\method{print}{tpp}(x)
\method{print}{tppint}(x)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class tpp or tppint.
}
\item{...}{graphics parameters passed to plot function.}
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}
\examples{
X <- tpp(sample(c(1:24),200,replace = T))
plot(X)
plot(density(X))
}

