## Marginals

## Marginal must have class marginal.mr and the following elements:
## - start(y,x): compute initial estimates ignoring correlations
## - dp(y,x,lambda): evaluate [d,p]
## - q(u,y,x,lambda): evaluate quantiles
## - npar(x): number of parameters
## - type: response type (integer or numeric)

# Gaussian
gaussian.marg <- function(link = "identity" ) {
  fm <- gaussian( substitute( link ) )
  ans <- list()
  ans$start <- function(y, x, offset) {
    m <- glm.fit( x , y, offset=offset, family=fm )
    ans <- c( coef( m ) , log(sqrt( mean( residuals( m )^2 ) ) ))
    names(ans) <- c(dimnames(as.matrix(x))[[2L]],"log.sigma")
    ans
  }
  ans$npar <- function(x) NCOL(x)+1
  ans$dp <- function(y, x, offset, lambda) {
    nb <- length(lambda)
    beta <- lambda[-nb]
    mu <- fm$linkinv( x %*% beta + offset )
    sd <- exp( lambda[nb] )
    cbind( dnorm( y , mu , sd ) , pnorm( y , mu , sd ) )
  }
  ans$sim <- function(u, x, offset, lambda) {
    nb <- length(lambda)
    beta <- lambda[-nb]
    mu <- fm$linkinv( x %*% beta + offset )
    sd <- exp( lambda[nb] )
    qnorm( u , mu , sd )
  }
  ans$type <- "numeric"
  class(ans) <- c( "marginal.gcmr")
  ans
}

# Binomial
# y is [#success,#failure]
binomial.marg <- function(link = "logit") {
  fm <- binomial( substitute( link ) )
  ans <- list()
  sizes <- 1
  ans$init <- function(y, x=NULL){
    if(NCOL(y)==1){
      y <- as.factor(y)
      y <- y!=levels(y)[1L]
      y <- cbind(y, 1-y)
    } 
    sizes <<- y[,1]+y[,2]
  }
  ans$start <- function(y, x, offset) {
    res <- coef(glm.fit( x, y, offset=offset, family=fm ))
    names(res) <- dimnames(as.matrix(x))[[2L]]
    res
  }
  ans$npar <- function(x) NCOL(x)
  ans$dp <- function(y, x, offset, lambda) {
    if(NCOL(y)==1){
      y <- as.factor(y)
      y <- y!=levels(y)[1L]
      y <- cbind(y, 1-y)
    }
    beta <- lambda
    mu <- fm$linkinv( x %*% beta + offset )
    cbind(dbinom( y[,1], sizes, mu ) ,
          pbinom( y[,1], sizes, mu ) )
  }
  ans$sim <- function(u, x, offset, lambda) {
    if(NCOL(y)==1){
      y <- as.factor(y)
      y <- y!=levels(y)[1L]
      y <- cbind(y, 1-y)
    }
    beta <- lambda
    mu <- fm$linkinv( x %*% beta + offset )
    sim <- qbinom( u, sizes, mu )
    cbind( sim, sizes-sim )
  }
  ans$type <- "integer"
  class(ans) <- c( "marginal.gcmr")
  ans
}

# Poisson
poisson.marg <- function(link = "log") {
  fm <- poisson( substitute( link ) )
  ans <- list()
  ans$start <- function(y, x, offset) {
    res <- coef(glm.fit( x , y, offset=offset, family=fm ))
    names(res) <- dimnames(as.matrix(x))[[2L]]
    res
  }
  ans$npar <- function(x) NCOL(x)
  ans$dp <- function(y, x, offset, lambda) {
    beta <- lambda
    mu <- fm$linkinv( x %*% beta + offset )
    cbind( dpois( y , mu ) , ppois( y , mu ) )
  }
  ans$sim <- function(u, x, offset, lambda) {
    beta <- lambda
    mu <- fm$linkinv( x %*% beta + offset )
    qpois( u , mu )
  }
  ans$type <- "integer"
  class(ans) <- c( "marginal.gcmr")
  ans
}


# Negative binomial
# var(y) = E(y) + k*E(y)^2 (k>0)
negbin.marg <- function(link = "log" ) {
  fm <- poisson( substitute( link ) )
  eps <- .Machine$double.eps
  ans <- list()
  ans$start <- function(y, x, offset) {
    m <- glm.fit( x , y, offset=offset, family=fm )
    mu <- fitted(m)
    res <- c( coef( m ) , max(0.01, mean(((y-mu)^2-mu)/mu^2)))
    names(res) <- c(dimnames(as.matrix(x))[[2L]],"dispersion")
    res
  }
  ans$npar <- function(x) NCOL(x)+1
  ans$dp <- function(y, x, offset, lambda) {
    nb <- length(lambda)
    beta <- lambda[-nb]
    if (lambda[nb]<=0) rep(NA,NROW(y))
    mu <- fm$linkinv( x %*% beta + offset )
    size <- pmax(eps, 1 / lambda[nb])
    cbind( dnbinom( y, mu=mu, size=size) , pnbinom( y, mu=mu, size=size) )
  }
  ans$sim <- function(u, x, offset, lambda) {
    nb <- length(lambda)
    beta <- lambda[-nb]
    if (lambda[nb]<=0) rep(NA,NROW(x))
    mu <- fm$linkinv( x %*% beta + offset )
    size <- 1 / beta[nb]
    qnbinom( u, mu=mu, size=size)
  }
  ans$type <- "integer"
  class(ans) <- c( "marginal.gcmr")
  ans
}

# next lines neeeded for back-compatibility with gcmr version 0.3
gs.marg <- gaussian.marg
bn.marg <- binomial.marg
ps.marg <- poisson.marg
nb.marg <- negbin.marg

# Weibull
weibull.marg <- function(link = "log"){
  fm <- Gamma( substitute( link ) ) # ;-)
  ans <- list()
  ans$start <- function(y, x, offset){
    nb <- length(beta)
    # starting values based on gamma glm
    m <- glm.fit(x , y, offset=offset, family=fm)
    disp <- sum( residuals(m, "pearson")^2 )/( NROW(y)-NCOL(x) )
    shape <- log( log(2) )/log( median( y/fitted(m) ) )
    ans <- c( coef(m), log(shape) )
    names(ans) <- c(dimnames(as.matrix(x))[[2L]], "log.shape")
    ans
  }
  ans$npar <- function(x) NCOL(x)+1
  ans$dp <- function(y, x, offset, lambda){
    nb <- length(lambda)
    beta <- lambda[-nb]
    scale <- fm$linkinv( x %*% beta + offset )
    shape <- exp( lambda[nb] )
    cbind( dweibull(y, shape=shape, scale=scale) , pweibull(y, shape=shape, scale=scale) )
  }
  ans$sim <- function(u, x, offset, lambda){
    nb <- length(lambda)
    beta <- lambda[-nb]
    scale <- fm$linkinv( x %*% beta + offset )
    shape <- exp( lambda[nb] )
    qweibull(u, shape=shape, scale=scale)
  }
  ans$type <- "numeric"
  class(ans) <- c( "marginal.gcmr")
  ans
} 
 
# Gamma
Gamma.marg <- function(link = "inverse"){
  fm <- Gamma( substitute( link ) )
  ans <- list()
  ans$start <- function(y, x, offset) {
    m <- glm.fit(x , y, offset=offset, family=fm)
    disp <- sum( residuals(m, "pearson")^2 )/( NROW(y)-NCOL(x) )
    shape <- 1/disp
    ans <- c( coef(m), log(shape) )
    names(ans) <- c(dimnames(as.matrix(x))[[2L]], "log.shape")
    ans
  }
  ans$npar <- function(x) NCOL(x)+1
  ans$dp <- function(y, x, offset, lambda){
    nb <- length(lambda)
    beta <- lambda[-nb]
    mu <- fm$linkinv( x %*% beta + offset )
    shape <- exp( lambda[nb] )
    cbind( dgamma(y, shape=shape, rate=shape/mu) , pgamma(y, shape=shape, rate=shape/mu) )
  }
  ans$sim <- function(u, x, offset, lambda){
    nb <- length(lambda)
    beta <- lambda[-nb]
    mu <- fm$linkinv( x %*% beta + offset )
    shape <- exp( lambda[nb] )
    qgamma(u, shape=shape, rate=shape/mu)
  }
  ans$type <- "numeric"
  class(ans) <- c( "marginal.gcmr")
  ans
}

 



 
