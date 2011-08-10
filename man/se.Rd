\name{se}
\alias{se}
\title{
Standard Errors of Parameters Estimators in Gaussian Copula Marginal Regression Models
}
\description{
Returns the standard errors vector of a fitted Gaussian copula marginal regression model.
}
\usage{
se(x, ...)
}
\arguments{
  \item{x}{
a fitted Gaussian copula marginal regression model of class \code{\link{gcmr}}.
}
  \item{\dots}{
further arguments passed to function \code{\link{vcov.gcmr}}.
}
}
\value{
A vector containing the estimated standard errors of the parameter of a fitted \code{\link{gcmr}} model.
}
\references{
Masarotto, G. and Varin. C. (2011). Gaussian copula marginal regression. Preprint. 
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}, \code{\link{vcov.gcmr}}
}
\examples{
data(polio)
## ARMA(2,1)      
m <- gcmr(y~., data=polio, marginal=nb.marg(), cormat=arma.cormat(2,1), 
          options=list(seed=71271, nrep=100))
m
## various type of standard errors
res <- cbind(se(m,"he"),se(m,"v"),se(m,"s"),se(m,"hac"))
colnames(res) <- c("hessian", "vscore", "sandwich", "HAC")
res
}
\keyword{regression}
\keyword{nonliner}
