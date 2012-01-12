\name{matern.cormat}
\alias{matern.cormat}
\title{
Gaussian Copula Regression with Matern Spatial Correlation Matrix
}
\description{
Sets a Matern spatial correlation matrix in Gaussian copula marginal regression models.
}
\usage{
matern.cormat(D, alpha=0.5)
}
\arguments{
  \item{D}{
matrix with values of the distances between pairs of data locations.
}
  \item{alpha}{
value of the smoothness parameter of the Matern correlation class.
}
}
\value{
An object of class \code{\link{cormat.gcmr}} representing a Matern correlation matrix.
}
\details{
The Matern correlation function is inherited from the \code{geoR} package (Diggle and Ribeiro, 2007).
}
\references{
Diggle, P. and Ribeiro, P.J. (2007). \emph{Model-based Geostatistics}. Springer.

Masarotto, G. and Varin, C. (2011). Gaussian copula marginal regression. Preprint. 
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}, \code{\link{matern}}
}
\keyword{regression}
\keyword{nonlinear}