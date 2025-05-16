\name{cluster.cormat}
\alias{cluster.cormat}
\title{
Longitudinal/Clustered Data Correlation
}
\description{
Sets longitudinal/clustered data correlation in Gaussian copula regression models.
}
\usage{
cluster.cormat(id, type = c("independence", "ar1", "ma1", 
		   	    "exchangeable", "unstructured"))
}
\arguments{
  \item{id}{
subject id. This is a vector of the same lenght of the number of observations. Please note that data must be
sorted in way that observations from the same cluster are contiguous.
}
  \item{type}{
a character string specifying the correlation structure. At the moment, the following are implemented: 
\tabular{ll}{
\code{independence} \tab working independence. \cr
\code{ar1} \tab autoregressive of order 1. \cr
\code{ma1} \tab moving average of order 1. \cr
\code{exchangeable} \tab exchangeable. \cr
\code{unstructured} \tab unstructured.
}
}
}
\value{
An object of class \code{\link{cormat.gcmr}} representing a correlation matrix for longitudinal or clustered data.
}
\details{
The correlation matrices are inherited from the \code{\link[nlme]{nlme}} package (Pinheiro and Bates, 2000).
}
\references{
Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549.

Masarotto, G. and Varin C. (2017). Gaussian Copula Regression in R. \emph{Journal of Statistical Software}, \bold{77}(8), 1--26.

Pinheiro, J.C. and Bates, D.M. (2000). \emph{Mixed-Effects Models in S and S-PLUS}. Springer.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}, \code{\link[nlme]{nlme}}.
}
\keyword{regression}
\keyword{nonlinear}