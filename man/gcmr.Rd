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
gcmr(formula, data, subset, offset, marginal,
     cormat, start, fixed, options=gcmr.options(...), model=TRUE,...)

gcmr.fit(x=rep(1,NROW(y)), y, z=NULL, offset=NULL,
	 marginal, cormat, start, fixed, options=gcmr.options())		   
}

\arguments{
\item{formula}{
a symbolic description of the model to be fitted of type \code{y ~ x} or 
\code{y ~ x | z}, for details see below.
}

\item{data}{
an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model.  If not found in \code{data}, the variables are taken from \code{environment(formula)}.
}

\item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}

\item{offset}{
optional numeric vector with an a priori known component to be included in the linear predictor for the mean. When appropriate, offset may also be a list of two offsets for the mean and precision equation, respectively.
}

\item{x}{design matrix.
}

\item{y}{vector of observations.
}

\item{z}{optional design matrix for the dispersion/shape.
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
\item{model}{
logical. If \code{TRUE}, then the model frame is returned. 
}

\item{...}{
arguments passed to \code{\link{gcmr.options}}.
}
}

\details{
Function \code{\link{gcmr}} computes maximum likelihood estimation in Gaussian copula marginal regression models. Computation of the exact likelihood is possible only for continuous responses, otherwise the likelihood function is approximated by importance sampling. See Masarotto and Varin (2012) for details.

Standard formula \code{y ~ x1 + x2} indicates that the mean response is modelled as a function of covariates \code{x1} and \code{x2} through an appropriate link function. Extended formula  \code{y ~ x1 + x2 | z1 + z2} indicates that the dispersion (or the shape) parameter of the marginal distribution is modelled as a function of covariates \code{z1} and \code{z2}. Dispersion (or shape) parameters are always modelled on logarithm scale. The model specification is inspired by beta regression as implemented in \code{\link{betareg}} (Cribari-Neto and Zeileis, 2010) through extended \code{\link{Formula}} objects (Zeileis and Croissant, 2010). 

For binomial marginals specified by \code{\link{binomial.marg}} the response is specified as a factor when the first level denotes failure and all others success or as a two-column matrix with the columns giving the numbers of successes and failures. 

\code{gcmr.fit} is the workhorse function: it is not normally called directly but can be more efficient where the response vector and design matrix have already been calculated.
}
\value{
An object of class \code{"gcmr"} with the following components:
\item{estimate}{the maximum likelihood estimate.}
\item{maximum}{the maximum likelihood value.}
\item{hessian}{(minus) the Hessian at the maximum likelihood estimate.}
\item{jac}{the Jacobian at the maximum likelihood estimate.}
\item{fitted.values}{the fitted values.}
\item{marginal}{the marginal model used.}
\item{cormat}{the correlation matrix used.}
\item{fixed}{the numeric vector indicating which parameters are constants.}
\item{ibeta}{the indices of marginal parameters.}
\item{igamma}{the indices of dependence parameters.}
\item{nbeta}{the number of marginal parameters.}
\item{ngamma}{the number of dependence parameters.}
\item{options}{the fitting options used, see \code{\link{gcmr.options}}.}
\item{call}{the matched call.}
\item{formula}{the model formula.}
\item{terms}{the terms objects for the fitted model.}
\item{levels}{the levels of the categorical regressors.}
\item{model}{the model frame, returned only if \code{model=TRUE}.}
\item{contrasts}{the contrasts corresponding to \code{levels}.}
\item{y}{the y vector used.}
\item{x}{the model matrix used for the mean response.}
\item{z}{the (optional) model matrix used for the dispersion/shape.}
\item{offset}{the offset used.}
\item{n}{the number of observations.}
\item{not.na}{the vector of binary indicators of the available observations (not missing).}


Functions \code{\link{coefficients}}, \code{\link{logLik}}, \code{\link{fitted}}, \code{\link{vcov.gcmr}} and \code{\link{residuals.gcmr}} can be used to extract various useful features of the value returned by \code{\link{gcmr}}. Function \code{\link{plot.gcmr}} produces various diagnostic plots for fitted \code{gcmr} objects.
}
\references{
Cribari-Neto, F. and Zeileis, A. (2010). Beta regression in R. \emph{Journal of Statistical Software} \bold{34}, 1--24. \url{http://www.jstatsoft.org/v34/i02/}.

Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.

Rocha, A.V. and Cribari-Neto, F. (2009). Beta autoregressive moving average models. \emph{Test} \bold{18}, 529--545.

Zeileis, A. and Croissant, Y. (2010). Extended model formulas in R: Multiple parts and multiple responses. \emph{Journal of Statistical Software} \bold{34}, 1--13. \url{http://www.jstatsoft.org/v34/i01/}.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{cormat.gcmr}}, \code{\link{marginal.gcmr}}, \code{\link{gcmr.options}}, \code{\link{Formula}}, \code{\link{betareg}}.
}
\examples{
## negative binomial model for longitudinal data
data(epilepsy)
gcmr(counts ~ offset(log(time)) + visit + trt + visit:trt, data = epilepsy, 
subset = (id != 49), marginal = negbin.marg, cormat = cluster.cormat(id, "ar1"), 
options=gcmr.options(seed=123, nrep=100 ))
## Hidden Unemployment Rate (HUR) data (Rocha and Cribari-Neto, 2009)
## beta regression with ARMA(1,3) errors
data(HUR)
trend <- scale(time(HUR))
gcmr(HUR ~ trend | trend, marginal = beta.marg, cormat = arma.cormat(1, 3))
}
\keyword{regression}
\keyword{nonlinear}
