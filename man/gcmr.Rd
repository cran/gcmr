\name{gcmr}
\alias{gcmr}
\alias{gcmr.fit}
\title{
Fitting Gaussian Copula Marginal Regression Models by Maximum (Simulated) Likelihood.
}
\description{
Fits Gaussian copula marginal regression models by maximum (simulated) likelihood.   
}
\usage{

gcmr(formula, data, subset, offset, contrasts=NULL, 
     marginal, cormat, start, fixed, options=gcmr.options())

gcmr.fit(x=rep(1,NROW(y)), y, offset=rep(0,NROW(y)), 
	 marginal, cormat, start, fixed, options=gcmr.options())

}

\arguments{
\item{formula}{
an object of class \code{"\link{formula}"} (or one that
can be coerced to that class): a symbolic description of the model to be fitted.
}

\item{data}{
an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model.  If not found in \code{data}, the variables are taken from \code{environment(formula)}.
}

\item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}

\item{offset}{
this can be used to specify an \emph{a priori} known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more \code{\link{offset}} terms can be included in the formula instead or as well, and if more than one is specified their sum is used.  See \code{\link{model.offset}}.
}

\item{contrasts}{
an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}.
}

\item{x}{a design matrix of dimension \code{n * p}.
}

\item{y}{vector of observations of length \code{n}.
}

\item{marginal}{
an object of class \code{\link{marginal.gcmr}} specifying the marginal part of the model. 
}

\item{cormat}{
an object of class \code{\link{cormat.gcmr}} representing the correlation matrix of the errors. 
}

\item{start}{
optional numeric vector with starting values for the model parameters.
}

\item{fixed}{
optional numeric vector of the same length as the total number of parameters. If supplied, only \code{NA} entries in fixed will be varied. 
}

\item{options}{
list of options passed to function \code{\link{gcmr.options}}.
}
}

\details{
Gaussian copula marginal regression models (Song, 2000; Masarotto and Varin, 2011) provide a flexible general framework for modelling dependent responses of any type. Gaussian copulas combine the simplicity of interpretation in marginal modelling with the flexibility in the specification of the dependence structure in multivariate normal models. 

This package contains \R functions related to the paper Masarotto and Varin (2011). The main function is \code{\link{gcmr}} that fits Gaussian copula marginal regression models. Inference is performed through a likelihood approach. Computation of the exact likelihood is possible only for continuous responses, otherwise the likelihood function is approximated by importance sampling. See Masarotto and Varin (2011) for details.

For binomial marginals specified by \code{bn.marg} the response is specified as a factor when the first level denotes failure and all others success or as a two-column matrix with the columns giving the numbers of successes and failures. 

\code{gcmr.fit} is the workhorse function: it is not normally called directly but can be more efficient where the response vector and design matrix have already been calculated.
}
\value{
An object of class \code{"gcmr"} with the following components:
\item{estimate}{the vector of parameter estimates.}
\item{maximum}{the maximum (simulated) likelihood.}
\item{hessian}{(minus) the Hessian at the maximum likelihood estimate.}
\item{jac}{the Jacobian at the maximum likelihood estimate.}
\item{y}{the y vector used.}
\item{x}{the model matrix used.}
\item{offset}{the offset used.}
\item{n}{the number of observations.}
\item{call}{the matched call.}
\item{not.na}{the vector of binary indicators of missing observations.}
\item{marginal}{the marginal model used.}
\item{cormat}{the correlation matrix used.}
\item{fixed}{the numeric vector indicating which parameters are constants.}
\item{ibeta}{the indices of marginal parameters.}
\item{igamma}{the indices of dependence parameters.}
\item{options}{the fitting options used, see \code{\link{gcmr.options}}.}

Functions \code{\link{coefficients}}, \code{\link{logLik}}, \code{\link{vcov.gcmr}}, \code{\link{se}} and \code{\link{residuals.gcmr}} can be used to extract various useful features of the value returned by \code{\link{gcmr}}.
}
\references{
Masarotto, G. and Varin, C. (2011). Gaussian copula marginal regression. Preprint. 

Song, P. X.-K. (2000). Multivariate dispersion models generated from Gaussian copula. \emph{Scandinavian Journal of Statistics} \bold{27}, 305--320.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{cormat.gcmr}}, \code{\link{marginal.gcmr}}, \code{\link{gcmr.options}}
}
\examples{
## Warning: Likelihood approximated using only a limited number 
## of Monte Carlo replications.
## Polio data. Marginal negative binomial model with ARMA(2,1) correlation matrix.
data(polio)
gcmr(y~., data=polio, marginal=negbin.marg, cormat=arma.cormat(2,1), 
     options=list(seed=71271, nrep=100))
## Scotland lip cancer data. 
## Marginal negative binomial model with Matern correlation matrix.
data(scotland)
D.scotland <- spDists(cbind(scotland$longitude, scotland$latitude), longlat=TRUE)
gcmr(observed~offset(log(expected))+AFF+I(latitude/100), data=scotland, 
     marginal=negbin.marg, cormat=matern.cormat(D.scotland), 
     options=list(seed=71271, nrep=100))
}
\keyword{regression}
\keyword{nonlinear}
