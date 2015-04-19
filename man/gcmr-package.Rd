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
Gaussian copula regression models (Song, 2000; Pitt, 2006; Masarotto and Varin, 2012) provide a mathematically convenient framework to account for dependence in regression analysis. Gaussian copulas allow to combine the simplicity of interpretation of marginal regression with the flexibility in the specification of the dependence structure of multivariate normal models. 

This package contains \R functions related to the papers Masarotto and Varin (2012) and Guolo and Varin (2014). The main function is \code{\link{gcmr}}, which fits Gaussian copula marginal regression models. Inference is performed through a likelihood approach. Computation of the exact likelihood is possible only for continuous responses, otherwise the likelihood function is approximated by importance sampling. See Masarotto and Varin (2012) for details.

}
\author{
Guido Masarotto and Cristiano Varin.
}
\references{
Guolo, A. and Varin, C. (2014). Beta regression for time series analysis of bounded data, with application to Canada Google Flu Trends. \emph{The Annals of Applied Statistics} \bold{8}, 74--88.

Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.

Pitt, M., Chan D. and Kohn, R. (2006). Efficient Bayesian inference for Gaussian copula regression models. \emph{Biometrika} \bold{93}, 537--554.

Song, P. X.-K. (2000). Multivariate dispersion models generated from Gaussian copula. \emph{Scandinavian Journal of Statistics} \bold{27}, 305--320.
}
\keyword{package}
