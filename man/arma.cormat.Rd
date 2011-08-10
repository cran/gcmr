\name{arma.cormat}
\alias{arma.cormat}
\title{
Gaussian Copula Regression with ARMA(p,q) Correlation Matrix
}
\description{
Sets ARMA(p,q) correlation in Gaussian copula regression models.
}
\usage{
arma.cormat(p, q)
}
\arguments{
  \item{p}{
order of the autoregressive component.
}
  \item{q}{
order of the moving average component.
}
}
\value{
An object of class \code{\link{cormat.gcmr}} representing a correlation matrix with ARMA(p,q) structure.
}
\references{
Masarotto, G. and Varin, C. (2011). Gaussian copula marginal regression. Preprint.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}
}
\keyword{regression}
\keyword{nonlinear}
