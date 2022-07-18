\name{arma.cormat}
\alias{arma.cormat}
\title{
ARMA(p,q) Correlation 
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
Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. 

Masarotto, G. and Varin C. (2017). Gaussian Copula Regression in R. \emph{Journal of Statistical Software}, \bold{77}(8), 1--26.  
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}.
}
\keyword{regression}
\keyword{nonlinear}
