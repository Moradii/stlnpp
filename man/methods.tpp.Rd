\name{methods.tpp}
\alias{methods.tpp}
\alias{print.tpp}
\alias{plot.tpp}
\title{
Methods for temporal point patterns.
}
\description{
Methods for temporal point patterns.
}
\usage{
\method{plot}{tpp}(x)
\method{print}{tpp}(x)
\method{plot}{tppint}(x)
\method{print}{tppint}(x)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class tpp or tppint.}
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}
\examples{
X <- tpp(sample(c(1:24),200,replace = T))
plot(X)
plot(density(X))
}

