\name{hausman}
\alias{hausman}
\title{
 Hausman-type specification test
}
\description{
Computes the Hausman-type specification test for validation of the multivariate normal assumption in Gaussian copula marginal regression models. The p-value of the test is estimated via parametric Bootstrap. 
}
\usage{
hausman(m, method = c("one-step", "full"), 
        nboot = if (method == "one-step") 1000 else 200, 
        nrep = 200)
}
\arguments{
  \item{m}{
a fitted marginal regression model of class \code{\link{gcmr}}.
}
  \item{method}{
\code{"one-step"} for a first-order approximation of the solution or \code{"full"} for a complete optimization for each Bootstrap sample. 
}
  \item{nboot}{
number of Bootstrap replicates. 
}
  \item{nrep}{
sets the value of parameter \code{nrep} in function \code{\link{gcmr.options}}.
}
}
\details{
The option \code{method="full"} gives more precise answers but at considerable higher computational cost. 
}
\value{
A list with the following components:
\tabular{ll}{
\code{overall} \tab the overall Hausman-type specification test. \cr
\code{parameters} \tab the Hausman-type specification test for each single marginal parameter. \cr
\code{options} \tab the options used (method, nboot, nrep).
}
}
\references{
Hausman, J.A. (1978). Specification tests in econometrics. \emph{Econometrica} \bold{46}, 1251--1271.

Masarotto, G. and Varin, C. (2011). Gaussian copula marginal regression. Preprint.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}
}
\keyword{regression}
\keyword{nonlinear}
