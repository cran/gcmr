\name{cluster.cormat}
\alias{cluster.cormat}
\title{
Gaussian Copula Regression with Longitudinal/Clustered Data Correlation Matrix
}
\description{
Sets longitudinal/clustered data correlation in Gaussian copula regression models.
}
\usage{
cluster.cormat(id, type = c("ar1", "ma1", "exch", "unstr"))
}
\arguments{
  \item{id}{
subject id. This is a vector of the same lenght of the number of observations. Please note that data must be
sorted in way that observations from the same cluster are contiguous.
}
  \item{type}{
a character string specifying the correlation structure. At the moment, the following are implemented: 
\tabular{ll}{
\code{ar1} \tab autoregressive of order 1. \cr
\code{ma1} \tab moving average of order 1. \cr
\code{exch} \tab exchangeable. \cr
\code{unstr} \tab unstructured.
}
}
}
\value{
An object of class \code{\link{cormat.gcmr}} representing a correlation matrix for longitudinal or clustered data.
}
\details{
The correlation matrices are inherited from the \code{\link{nlme}} package (Pinheiro and Bates, 2000).
}
\references{
Masarotto, G. and Varin, C. (2011). Gaussian copula regression models. Preprint.

Pinheiro, J.C. and Bates, D.M. (2000). \emph{Mixed-Effects Models in S and S-PLUS}. Springer.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}, \code{\link{nlme}}
}
\keyword{regression}
\keyword{nonlinear}