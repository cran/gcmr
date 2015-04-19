\name{summary.gcmr}

\alias{summary.gcmr}
\alias{print.summary.gcmr}

\title{
Summarizing Gaussian Copula Marginal Regression Models
}
\description{
These functions are all \code{\link{methods}} for class \code{\link{gcmr}} or \code{summary.gcmr} objects. 
}
\usage{
\method{summary}{gcmr}(object, type = "hessian", \dots)
\method{print}{summary.gcmr}(x, digits = max(3, getOption("digits") - 3), \dots) 
}
\arguments{
\item{object}{a fitted marginal regression model of class \code{\link{gcmr}}.}
\item{x}{an object of class "summary.glm", usually, a result of a call to \code{summary.glm}.}
\item{type}{the method used to compute standard errors, see \code{\link{vcov.gcmr}}.}
\item{digits}{the number of significant digits to use when printing.}
\item{\dots}{further arguments passed to or from other methods.}
}
\value{
The function \code{summary.gcmr} returns an object of class "summary.glm", a list with some components of the \code{gcmr} object, plus 
\item{se.type}{the method used to compute standard errors, see \code{\link{vcov.gcmr}}.}
\item{coefficients}{a list with components \code{marginal} and \code{copula} containing the maximum likelihood estimates of the marginal and Gaussian copula parameters, respectively.}
\item{aic}{Akaike Information Criterion.}
}
\seealso{\code{\link{gcmr}}}
\examples{
data(epilepsy)
fit <- gcmr(counts ~ offset(log(time)) + visit + trt + visit:trt, data = epilepsy, 
subset = (id != 49), marginal = negbin.marg, cormat = cluster.cormat(id, "ar1"), 
options=gcmr.options(seed=123, nrep=c(25,100) ))
summary(fit)
}