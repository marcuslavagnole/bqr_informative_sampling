#' Bayesian weighted quantile regression - ALD
#'
#' This function estimates a Bayesian quantile regression model for complex 
#' survey data under informative sampling, where the response variable is 
#' assumed to follow an asymmetric Laplace distribution.
#'
#' @param y n-dimensional vector of responses.
#' @param x Matrix - nxp - of predictors (include the intercept)
#' @param w n-dimensional vector of survey weights.
#' @param tau Quantile of interest (value between 0 and 1)
#' @param n_mcmc Number of iterations.
#' @param burnin_mcmc Number of initial iterations to be discarded.
#' @param thin_mcmc Thinning parameter.
#' 
#' @return A list with the chains of all parameters of interest.

## Packages
require(mvtnorm); require(GIGrvg)

## MCMC
bayesQR_weighted <- function(y,x,w,tau,n_mcmc,burnin_mcmc,thin_mcmc){
  n         <- length(y)
  numcov    <- ncol(x)
  resultado <- list()
  # Normalize weights
  w <- n*(w/sum(w))
  # Create auxiliary objects
  beta  <- matrix(NA, n_mcmc, numcov)
  sigma <- matrix(NA, n_mcmc, 1)
  # Set the initial values
  reg_lm     <- lm(y~-1+x)
  beta[1,]   <- reg_lm$coefficients
  sigma[1,1] <- 1
  v          <- rgamma(n,2,1)
  # Auxiliary constants
  delta2  <- 2/(tau*(1-tau))
  theta   <- (1-2*tau)/(tau*(1-tau))
  # MCMC
  for(k in 2:n_mcmc){
    beta[k,]   <- atualizarBETA(rep(0,numcov),diag(rep(1000,numcov)),y,x,w,sigma[k-1,1],delta2,theta,v)
    v          <- atualizarV(y,x,w,beta[k,],delta2,theta,1,n)
    sigma[k,1] <- atualizarSIGMA(0.001,0.001,y,x,w,beta[k,],delta2,theta,v,n)
  }
  resultado[[1]]  <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[[2]]  <- sigma[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  
  return(resultado)
}

#######################################################################
## Auxiliary functions: Sampling from full conditional distributions ##
#######################################################################

### Full conditional for beta
atualizarBETA<-function(b,B,y,x,w,sigma,delta2,theta,v){
  B.inv  <- chol2inv(chol(B))
  covar  <- chol2inv(chol(B.inv+(1/(delta2*sigma)*(t((w/v)*x)%*%x))))
  media  <- covar%*%((B.inv%*%b)+(1/(delta2*sigma))*(t((w/v)*x)%*%(y-theta*v)))
  beta   <- rmvnorm(1,media,covar)
  return(beta)
}

### Full conditional for sigma
atualizarSIGMA<-function(c,C,y,x,w,beta,tau2,theta,v,n){
  alpha1 <- c + 1.5*n
  beta1  <- C + sum(w*v) + (t(w*(y-x%*%beta-theta*v)/v)%*%(y-x%*%beta-theta*v))/(2*tau2)
  sigma  <- 1/rgamma(1, alpha1, beta1)
  return(sigma)
}

### Full conditional for the latent variable
atualizarV<-function(y,x,w,beta,delta2,theta,sigma,n){
  p1 <- 0.5
  p2 <- w*(y-x%*%beta)^2/(delta2*sigma)
  p3 <- w*(2/sigma + theta^2/(delta2*sigma))
  v  <- NULL
  for(i in 1:n){
    v[i]  <- rgig(1, chi=p2[i], psi=p3[i], lambda=p1)
  }
  return(v)
}
