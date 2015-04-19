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
    no.se = FALSE, method = c("BFGS", "Nelder-Mead", "CG")) 
}
\arguments{
  \item{seed}{
seed of the pseudorandom generator used in the importance sampling algorithm  for likelihood approximation in case of discrete responses.
}
  \item{nrep}{
Monte Carlo size of the importance sampling algorithm  for likelihood approximation in case of discrete responses. \code{nrep} can be a vector so that the model is fitted with a sequence of different Monte Carlo sizes. In this case, the starting values for optimization of the likelihood are taken from the previous fitting. A reasonable strategy is to fit the model with a small Monte Carlo size to obtain sensible starting values and then refit with a larger Monte Carlo size. The default value is \code{100} for the first optimization and \code{1000} for the second and definitive optimization. 
}
  \item{no.se}{
logical. Should standard errors be computed and returned or not?
}
  \item{method}{
   a character string specifying the method argument passed to \code{optim}. The default optimization routine is the quasi-Newton algorithm \code{BFGS}. See \code{\link{optim}} for details. 
}
}
\value{
A list containing the options.
}
\references{
Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.
}
\author{
Guido Masarotto and Cristiano Varin. 
}
\seealso{
\code{\link{gcmr}}
}
\keyword{regression}
\keyword{nonlinear}
