\name{matern.cormat}
\alias{matern.cormat}
\title{
Matern Spatial Correlation
}
\description{
Sets a Matern spatial correlation matrix in Gaussian copula marginal regression models.
}
\usage{
matern.cormat(D, alpha = 0.5)
}
\arguments{
  \item{D}{
matrix with values of the distances between pairs of data locations.
}
  \item{alpha}{
value of the shape parameter of the Matern correlation class. The default {alpha = 0.5} corresponds to an exponential correlation model. 
}
}
\value{
An object of class \code{\link{cormat.gcmr}} representing a Matern correlation matrix.
}
\details{
The Mat\'ern correlation function is inherited from the \code{geoR} package (Diggle and Ribeiro, 2007).
}
\references{
Diggle, P. and Ribeiro, P.J. (2007). \emph{Model-based Geostatistics}. Springer.

Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}, \code{\link{matern}}.
}
\keyword{regression}
\keyword{nonlinear}