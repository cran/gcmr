\name{gcmr-package}
\alias{gcmr-package}
\docType{package}
\title{
Gaussian Copula Marginal Regression
}
\description{
 Fits Gaussian copula marginal regression models described by Song (2000) and Masarotto and Varin (2012). 
}
\details{
Gaussian copula marginal regression models (Song, 2000; Masarotto and Varin, 2012) provide a flexible general framework for modelling dependent responses of any type. Gaussian copulas combine the simplicity of interpretation in marginal modelling with the flexibility in the specification of the dependence structure in multivariate normal models. 

This package contains \R functions related to the paper Masarotto and Varin (2012). The main function is \code{\link{gcmr}} that fits Gaussian copula marginal regression models. Inference is performed through a likelihood approach. Computation of the exact likelihood is possible only for continuous responses, otherwise the likelihood function is approximated by importance sampling. See Masarotto and Varin (2012) for details.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\references{
Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics}, \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.

Song, P. X.-K. (2000). Multivariate dispersion models generated from Gaussian copula. \emph{Scandinavian Journal of Statistics} \bold{27}, 305--320.
}
\keyword{package}
