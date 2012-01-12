## Marginals

## Marginal must have class marginal.mr and the following elements:
## - start(y,x): compute initial estimates ignoring correlations
## - dp(y,x,beta): evaluate [d,p]
## - q(u,y,x,beta): evaluate quantiles
## - npar(x): number of parameters
## - type: response type (integer or numeric)

# Gaussian
gs.marg <- function(link = "identity" ) {
    fm <- gaussian( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x, offset) {
        m <- glm.fit( x , y, offset=offset, family=fm )
        ans <- c( coef( m ) , log(sqrt( mean( residuals( m )^2 ) ) ))
        names(ans) <- c(dimnames(as.matrix(x))[[2L]],"log.sigma")
        ans
    }
    ans$npar <- function(x) NCOL(x)+1
    ans$dp <- function(y, x, offset, beta) {
        nb <- length(beta)
        mu <- fm$linkinv( x %*% beta[-nb] + offset )
        sd <- exp( beta[nb] )
        cbind( dnorm( y , mu , sd ) , pnorm( y , mu , sd ) )
    }
    ans$sim <- function( u, y, x, offset, beta) {
        nb <- length(beta)
        mu <- fm$linkinv( x %*% beta[-nb] + offset )
        sd <- exp( beta[nb] )
        qnorm( u , mu , sd )
    }
    ans$type <- "numeric"
    class(ans) <- c( "marginal.gcmr")
    ans
}

# Binomial
# y is [#success,#failure]
bn.marg <- function(link = "logit") {
    fm <- binomial( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x, offset) {
        res <- coef(glm.fit( x , y, offset=offset, family=fm ))
        names(res) <- dimnames(as.matrix(x))[[2L]]
        res
    }
    ans$npar <- function(x) NCOL(x)
    ans$dp <- function(y, x, offset, beta) {
        if(NCOL(y)==1){
          y <- as.factor(y)
          y <- y!=levels(y)[1L]
          y <- cbind(y, 1-y)
        } 
        mu <- fm$linkinv( x %*% beta + offset )
        cbind(dbinom( y[,1] , y[,1]+y[,2], mu ) ,
              pbinom( y[,1] , y[,1]+y[,2] , mu ) )
    }
    ans$sim <- function(u, y, x, offset, beta) {
        if(NCOL(y)==1){
          y <- as.factor(y)
          y <- y!=levels(y)[1L]
          y <- cbind(y, 1-y)
        } 
        mu <- fm$linkinv( x %*% beta + offset )
        sim <- qbinom( u, y[,1]+y[,2] , mu )
        cbind( sim, y[,1]+y[,2]-sim )
    }
    ans$type <- "integer"
    class(ans) <- c( "marginal.gcmr")
    ans
}

# Poisson
ps.marg <- function(link = "log") {
    fm <- poisson( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x, offset) {
        res <- coef(glm.fit( x , y, offset=offset, family=fm ))
        names(res) <- dimnames(as.matrix(x))[[2L]]
        res
    }
    ans$npar <- function(x) NCOL(x)
    ans$dp <- function(y, x, offset, beta) {
        mu <- fm$linkinv( x %*% beta + offset )
        cbind( dpois( y , mu ) , ppois( y , mu ) )
    }
    ans$sim <- function(u, y, x, offset, beta) {
        mu <- fm$linkinv( x %*% beta + offset )
        qpois( u , mu )
    }
    ans$type <- "integer"
    class(ans) <- c( "marginal.gcmr")
    ans
}


# Negative binomial
# var(y) = E(y) + k*E(y)^2 (k>0)
nb.marg <- function(link = "log" ) {
    fm <- poisson( substitute( link ) )
    h <- 2
    h2 <- 2 - h
    eps <- .Machine$double.eps
    ans <- list()
    ans$start <- function(y, x, offset) {
        m <- glm.fit( x , y, offset=offset, family=fm )
        mu <- fitted(m)
        res <- c( coef( m ) , max(0.01, mean(((y-mu)^2-mu)/mu^h)))
        names(res) <- c(dimnames(as.matrix(x))[[2L]],"k")
        res
    }
    ans$npar <- function(x) NCOL(x)+1
    ans$dp <- function(y, x, offset, beta) {
        nb <- length(beta)
        if (beta[nb]<=0) rep(NA,NROW(y))
        mu <- fm$linkinv( x %*% beta[-nb] + offset )
        size <- pmax(eps, mu^h2 / beta[nb])
        cbind( dnbinom( y , mu=mu , size=size) , pnbinom( y , mu=mu , size=size) )
    }
    ans$sim <- function( u, y, x, offset, beta) {
        nb <- length(beta)
        if (beta[nb]<=0) rep(NA,NROW(y))
        mu <- fm$linkinv( x %*% beta[-nb] + offset )
        size <- mu^h2 / beta[nb]
        qnbinom( u , mu=mu, size=size)
    }
    ans$type <- "integer"
    class(ans) <- c( "marginal.gcmr")
    ans
}

 


 
