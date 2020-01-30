## Copyright (C) P. J. Ribeiro and P. J. Diggle
## Distributed under the GPL-2 or GPL-3 licence
## Originally included in package geoR
.matern <- function (u, phi, kappa) 
{
    if (is.vector(u)) 
        names(u) <- NULL
    if (is.matrix(u)) 
        dimnames(u) <- list(NULL, NULL)
    uphi <- u/phi
    uphi <- ifelse(u > 0, (((2^(-(kappa - 1)))/ifelse(0, Inf, 
        gamma(kappa))) * (uphi^kappa) * besselK(x = uphi, nu = kappa)), 
        1)
    uphi[u > 600 * phi] <- 0
    return(uphi)
}
