\name{profile.gcmr}
\alias{profile.gcmr}
\title{
Profile Log-Likelihood for Gaussian Copula Marginal Regression Models
}
\description{
Computes the profile log-likelihood for mean response parameters of a Gaussian copula marginal regression model.
}
\usage{
\method{profile}{gcmr}(fitted, which, low, up, npoints = 10, 
        display = TRUE, alpha = 0.05, progress.bar = TRUE, ...)
}
\arguments{
  \item{fitted}{
a fitted Gaussian copula marginal regression model of class \code{\link{gcmr}}.
}
  \item{which}{
the index of the regression parameter which should be profiled. 
}
  \item{low}{
the lower limit used in computation of the profile log-likelihood. If this is \code{missing}, then the lower limit is set equal to the estimate minus three times its standard error.
}
  \item{up}{
the upper limit used in computation of the profile log-likelihood. If this is \code{missing}, then the upper limit is set equal to the estimate plus three times its standard error.
}
  \item{npoints}{
number of points used in computation of the profile log-likelihood. Default is \code{10}.
}
  \item{display}{
should the profile log-likelihood be displayed or not? default is \code{TRUE}. 
}
  \item{alpha}{
the significance level, default is \code{0.05}. 
}
\item{progress.bar}{
logical. If TRUE, a text progress bar is displayed. 
}
\item{...}{
further arguments passed to \code{plot}.
}
}
\details{
If the display is requested, then the profile log-likelihood is smoothed by cubic spline interpolation. 
}
\value{
A list with the following components:
\item{points}{points at which the profile log-likelihood is evaluated.}
\item{profile}{values of the profile log-likelihood.}
}
\references{
Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. 

Masarotto, G. and Varin C. (2017). Gaussian Copula Regression in R. \emph{Journal of Statistical Software}, \bold{77}(8), 1--26. 
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}
}
\examples{
## spatial binomial data
\dontrun{
data(malaria)
D <- sp::spDists(cbind(malaria$x, malaria$y))/1000
m <- gcmr(cbind(cases, size-cases) ~ netuse+I(green/100)+phc, data=malaria,
marginal=binomial.marg, cormat=matern.cormat(D), options=gcmr.options(seed=987))
prof <- profile(m, which = 2)
prof
}
}
\keyword{regression}
\keyword{nonlinear}
