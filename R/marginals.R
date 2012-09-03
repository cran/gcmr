## Marginals

## Marginal must have class marginal.mr and the following elements:
## - start(y,x,offset): compute initial estimates ignoring correlations
##      optional attributes lower and upper can be used
##      to specify box constrained parameters
## - dp(y,x,offset,lambda): evaluate [d,p]
## - q(p,x,offset,lambda): evaluate quantiles
## - npar(x): number of parameters
## - type: response type (integer or numeric)

# Gaussian
gaussian.marg <- function(link = "identity" ) {
    fm <- gaussian( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x, offset) {
        eps <- sqrt(.Machine$double.eps)
        m <- glm.fit( x , y, offset=offset, family=fm )
        lambda <- c( coef( m ) , max(10*eps,sd(residuals(m))))
        names(lambda) <- c(dimnames(as.matrix(x))[[2L]],"sigma")
        attr(lambda,"lower") <- c(rep(-Inf,NCOL(x)), eps)
        lambda
    }
    ans$npar <- function(x) NCOL(x)+1
    ans$dp <- function(y, x, offset, lambda) {
        nb <- length(lambda)
        beta <- lambda[-nb]
        sd <- lambda[nb]
        mu <- fm$linkinv( x %*% beta + offset )
        cbind( dnorm( y , mu , sd ) , pnorm( y , mu , sd ) )
    }
    ans$q <- function(p, x, offset, lambda) {
        nb <- length(lambda)
        beta <- lambda[-nb]
        sd <- lambda[nb]
        mu <- fm$linkinv( x %*% beta + offset )
        qnorm( p , mu , sd )
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
    sizes  <- 1
    ans$start <- function(y, x, offset) {
        if(NCOL(y)==1) {
            y <- as.factor(y)
            y <- y!=levels(y)[1L]
            y <- cbind(y, 1-y)      
        } 
        sizes <<- y[,1]+y[,2]
        lambda <- coef(glm.fit( x, y, offset=offset, family=fm ))
        names(lambda) <- dimnames(as.matrix(x))[[2L]]
        lambda
    }
    ans$npar <- function(x) NCOL(x)
    ans$dp <- function(y, x, offset, lambda) {
        if(NCOL(y)==1){
            y <- as.factor(y)
            y <- y!=levels(y)[1L]
            y <- cbind(y, 1-y)
        }
        mu <- fm$linkinv( x %*% lambda + offset )
        cbind(dbinom( y[,1], sizes, mu ) ,
              pbinom( y[,1], sizes, mu ) )
    }
    ans$q <- function(p, x, offset, lambda) {
        mu <- fm$linkinv( x %*% lambda + offset )
        q <- qbinom( p, sizes, mu )
        cbind( q, sizes-q )
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
        lambda <- coef(glm.fit( x , y, offset=offset, family=fm ))
        names(lambda) <- dimnames(as.matrix(x))[[2L]]
        lambda
    }
    ans$npar <- function(x) NCOL(x)
    ans$dp <- function(y, x, offset, lambda) {
        mu <- fm$linkinv( x %*% lambda + offset )
        cbind( dpois( y , mu ) , ppois( y , mu ) )
    }
    ans$q <- function(p, x, offset, lambda) {
        mu <- fm$linkinv( x %*% lambda + offset )
        qpois( p , mu )
    }
    ans$type <- "integer"
    class(ans) <- c( "marginal.gcmr")
    ans
}


# Negative binomial
# var(y) = E(y) + k*E(y)^2 (k>0)
negbin.marg <- function(link = "log" ) {
    fm <- poisson( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x, offset) {
        eps <- sqrt(.Machine$double.eps)
        m <- glm.fit( x , y, offset=offset, family=fm )
        mu <- fitted(m)
        lambda <- c( coef( m ) , max(10*eps , mean(((y-mu)^2-mu)/mu^2)))
        names(lambda) <- c(dimnames(as.matrix(x))[[2L]],"dispersion")
        attr(lambda,"lower") <- c(rep(-Inf,NCOL(x)),eps)
        lambda
    }
    ans$npar <- function(x) NCOL(x)+1
    ans$dp <- function(y, x, offset, lambda) {
        nb <- length(lambda)
        beta <- lambda[-nb]
        size <- 1/lambda[nb]
        mu <- fm$linkinv( x %*% beta + offset )
        cbind( dnbinom( y, mu=mu, size=size) , pnbinom( y, mu=mu, size=size) )
    }
    ans$q <- function(p, x, offset, lambda) {
        nb <- length(lambda)
        beta <- lambda[-nb]
        size <- 1 / lambda[nb]
        mu <- fm$linkinv( x %*% beta + offset )
        qnbinom( p, mu=mu, size=size)
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
    ans$start <- function(y, x, offset) {
        eps <- sqrt(.Machine$double.eps)
        m <- glm.fit(x , y, offset=offset, family=fm)
        shape <- max(10*eps, 1.2/sqrt(mean(log(y/fitted(m))^2)))
        lambda <- c( coef(m), shape )
        names(lambda) <- c(dimnames(as.matrix(x))[[2L]], "shape")
        attr(lambda,"lower") <- c(rep(-Inf,NCOL(x)),eps)
        lambda
    }
    ans$npar <- function(x) NCOL(x)+1
    ans$dp <- function(y, x, offset, lambda){
        nb <- length(lambda)
        beta <- lambda[-nb]
        shape <- lambda[nb]
        scale <- fm$linkinv( x %*% beta + offset )
        cbind(dweibull(y, shape=shape, scale=scale) ,
              pweibull(y, shape=shape, scale=scale) )
    }
    ans$q <- function(p, x, offset, lambda){
        nb <- length(lambda)
        beta <- lambda[-nb]
        shape <- lambda[nb]
        scale <- fm$linkinv( x %*% beta + offset )
        qweibull(p, shape=shape, scale=scale)
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
        eps <- sqrt(.Machine$double.eps)
        m <- glm.fit(x , y, offset=offset, family=fm)
        disp <- sum( residuals(m, "pearson")^2 )/( NROW(y)-NCOL(x) )
        lambda <- c( coef(m), max(10*eps,1/disp) )
        names(lambda) <- c(dimnames(as.matrix(x))[[2L]], "shape")
        attr(lambda,"lower") <- c(rep(-Inf,NCOL(x)),eps)
        lambda
    }
    ans$npar <- function(x) NCOL(x)+1
    ans$dp <- function(y, x, offset, lambda){
        nb <- length(lambda)
        beta <- lambda[-nb]
        shape <- lambda[nb]
        mu <- fm$linkinv( x %*% beta + offset )
        cbind(dgamma(y, shape=shape, rate=shape/mu) ,
              pgamma(y, shape=shape, rate=shape/mu) )
    }
    ans$q <- function(p, x, offset, lambda){
        nb <- length(lambda)
        beta <- lambda[-nb]
        shape <- lambda[nb]
        mu <- fm$linkinv( x %*% beta + offset )
        qgamma(p, shape=shape, rate=shape/mu)
    }
    ans$type <- "numeric"
    class(ans) <- c( "marginal.gcmr")
    ans
}

 



 
