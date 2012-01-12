\name{gcmr.options}
\alias{gcmr.options}
\title{
Setting Options for Fitting Gaussian Copula Marginal Regression Models
}
\description{
Sets options which affect fitting of Gaussian copula marginal regression models.
}
\usage{
gcmr.options(seed = round(runif(1, 1, 1e+05)), nrep = c(100, 1000), 
             no.se = FALSE, opt = gcmr.opt)
}
\arguments{
  \item{seed}{
seed of the pseudorandom generator used in the importance sampling algorithm  for likelihood approximation in case of noncontinuous responses.
}
  \item{nrep}{
Monte Carlo size of the importance sampling algorithm  for likelihood approximation in case of noncontinuous responses. This could be a vector, in this case different Monte Carlo sizes are used in sequence with starting values for optimization of the likelihood taken from previous fitting. A good strategy could be to consider first a small Monte Carlo size to get sensible starting values and then move to a larger Monte Carlo size.
}
  \item{no.se}{
logical. Should standard errors be computed and returned or not?
}
  \item{opt}{
    the function used to optimize the likelihood. The default uses
    \code{nlminb}. See \code{gcmr:::gcmr.opt}.
}
}
\value{
A list containing the options.
}
\references{
Masarotto, G. and Varin. C. (2011). Gaussian copula marginal regression. Preprint.
}
\author{
Guido Masarotto and Cristiano Varin. 
}
\seealso{
\code{\link{gcmr}}
}
\keyword{regression}
\keyword{nonlinear}
