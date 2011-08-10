\name{profile.gcmr}
\alias{profile.gcmr}
\title{
Profile likelihood for Gaussian Copula Marginal Regression Models
}
\description{
Computes the profile likelihood for a fitted Gaussian copula marginal regression model.
}
\usage{
\method{profile}{gcmr}(fitted, which, low, up, npoints = 10, 
        display = TRUE, alpha = 0.05, ...)
}
\arguments{
  \item{fitted}{
a fitted Gaussian copula marginal regression model of class \code{\link{gcmr}}.
}
  \item{which}{
the index of the parameter which should be profiled. 
}
  \item{low}{
the lower limit used in computation of the profile likelihood.
}
  \item{up}{
the upper limit used in computation of the profile likelihood.
}
  \item{npoints}{
number of points used in computation of the profile likelihood. Default is \code{10}.
}
  \item{display}{
should the profile likelihood be displayed or not? default is \code{TRUE}. 
}
  \item{alpha}{
the significance level, default is \code{0.05}. 
}
\item{...}{
further argument. 
}
}
\details{
If display requested, then the profile likelihood is smoothed by cubic spline interpolation. 
}
\value{
A list with the following components:
\item{points}{points at which the profile likelihood is evaluated.}
\item{profile}{values of the profile likelihood.}
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}
}
\examples{
## marginal negative binomial model with Matern correlation
data(scotland)
D.scotland <- spDists(cbind(scotland$longitude, scotland$latitude), longlat=TRUE)
m <- gcmr(observed~offset(log(expected))+AFF+I(latitude/100), data=scotland, 
          marginal=nb.marg(), cormat=matern.cormat(D.scotland), 
          options=list(seed=71271, nrep=100))
prof <- profile(m, which=2, low=0, up=8)
prof
}
\keyword{regression}
\keyword{nonlinear}
