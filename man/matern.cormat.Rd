\name{matern.cormat}
\alias{matern.cormat}
\title{
Gaussian Copula Regression with Matern Spatial Correlation Matrix
}
\description{
Sets a Matern spatial correlation matrix in Gaussian copula marginal regression models.
}
\usage{
matern.cormat(D, k=0.5)
}
\arguments{
  \item{D}{
matrix with values of the distances between pairs of data locations.
}
  \item{k}{
value of the smoothness parameter of the Matern correlation class.
}
}
\value{
An object of class \code{\link{cormat.gcmr}} representing a Matern correlation matrix.
}
\details{
The correlation matrix is inherited from function \code{\link{matern}} of the geoR package.
}
\references{
Masarotto, G. and Varin, C. (2011). Gaussian copula marginal regression. Preprint.

Ribeiro Jr., P.J. and Diggle, P.J. (2001). geoR: A package for geostatistical analysis. \emph{R-NEWS} \bold{1}, 15--18. 
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{matern}}, \code{\link{gcmr}}
}
\keyword{regression}
\keyword{nonlinear}