\name{gcmr-package}
\alias{gcmr-package}
\docType{package}
\title{
Gaussian Copula Marginal Regression
}
\description{
 Fits Gaussian copula marginal regression models described in Song (2000) and Masarotto and Varin (2012; 2017). 
}
\details{
Gaussian copula models are frequently used to extend univariate regression models to the multivariate case. The principal merit of the approach is that the specification of the regression model is conveniently separated from the dependence structure described in the familiar form of the correlation matrix of a multivariate Gaussian distribution (Song 2000). This form of flexibility has been successfully employed in several complex applications including longitudinal data analysis, spatial statistics, genetics and time series. Some useful references can be found in Masarotto and Varin (2012; 2017), Song et al. (2013) and Nikoloulopoulos (2015). 

This package contains \R functions that implement the methodology discussed in Masarotto and Varin (2012) and Guolo and Varin (2014). The main function is \code{\link{gcmr}}, which fits Gaussian copula marginal regression models. Inference is performed through a likelihood approach. Computation of the exact likelihood is possible only for continuous responses, otherwise the likelihood function is approximated by importance sampling. See Masarotto and Varin (2017) for details.

}
\author{
Guido Masarotto and Cristiano Varin.
}
\references{
Guolo, A. and Varin, C. (2014). Beta regression for time series analysis of bounded data, with application to Canada Google Flu Trends. \emph{The Annals of Applied Statistics} \bold{8}, 74--88.

Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.

Masarotto, G. and Varin C. (2017). Gaussian Copula Regression in R. \emph{Journal of Statistical Software}, \bold{77}(8), 1--26. \doi{10.18637/jss.v077.i08}.

Nikoloulopoulos A (2015). Efficient Estimation of High-dimensional Multivariate Normal Copula Models with Discrete Spatial Responses. \emph{Stochastic Environmental Research and Risk Assessment}, to appear.

Song, P. X.-K. (2000). Multivariate dispersion models generated from Gaussian copula. \emph{Scandinavian Journal of Statistics} \bold{27}, 305--320.

Song, P. X.-K., Li, M. and Zhang, P. (2013). Copulae in Mathematical and Quantitative Finance. In \emph{Vector Generalized Linear Models: A Gaussian Copula Approach}, 251--276. Springer Berlin Heidelberg.
}
\keyword{package}
