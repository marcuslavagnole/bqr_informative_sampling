#' Bayesian weighted Quantile Regression - Approximate Likelihood
#'
#' This function estimates a Bayesian quantile regression model for complex 
#' survey data under informative sampling, where an approximate likelihood
#' is employed.
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
require(mvtnorm)

## MCMC
bayesQRAP_weighted <- function(y,x,w,tau,n_mcmc,burnin_mcmc,thin_mcmc){
  n         <- length(y)
  numcov    <- ncol(x)
  resultado <- list()
  # Create auxiliary objects
  beta  <- matrix(NA, n_mcmc, numcov)
  # Set the initial values
  reg_lm     <- lm(y~-1+x)
  beta[1,]   <- reg_lm$coefficients
  const      <- NULL
  const[1]   <- 1
  # MCMC
  for(k in 2:n_mcmc){
    beta_aux   <- atualizarBETA_MH_bio(rep(0,numcov),diag(rep(1000,numcov)),y,x,w,beta[k-1,],tau,const[k-1],k)
    beta[k,]   <- beta_aux[1:numcov] ; const[k] <- beta_aux[(numcov+1)] 
  }
  resultado[[1]]  <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  
  return(resultado)
}


######################################################################
## Auxiliary functions: Sampling from full conditional distribution ##
######################################################################

### Metropolis-Hastings for the coefficients
atualizarBETA_MH_bio <- function(b,B,y,x,w,beta,tau,ct,k){
  valoratual    <- beta
  sigma_MH      <- (tau*(1-tau))*chol2inv(chol((1/length(y))*t(w^2*x)%*%x))
  valorproposto <- as.vector(rmvnorm(1, mean=valoratual, sigma=ct*sigma_MH))
  candidato     <- exp(condicionalBETA_MH_bio(valorproposto,b,B,y,x,w,tau)-
                       condicionalBETA_MH_bio(valoratual,b,B,y,x,w,tau))
  
  chanceaceitar <- min(1,candidato)
  if(runif(1)<chanceaceitar){
    BETAfinal   <- valorproposto
  } else{
    BETAfinal   <- valoratual
  }
  log_ct        <- log(ct)+(1/k^0.8)*(chanceaceitar-0.234)
  
  return(c(BETAfinal,exp(log_ct)))
}

### Full conditional for the coefficients
condicionalBETA_MH_bio <- function(beta,b,B,y,x,w,tau){
  n       <- length(y)
  
  B_inv   <- chol2inv(chol(B))
  priori  <- -0.5*(t(beta-b)%*%B_inv%*%(beta-b))
  
  s_tau   <- t(w*x)%*%(tau-((y-x%*%beta)<0))
  s_tau_i <- x*as.vector(tau-((y-x%*%beta)<0))
  w_cov_aux<-t(((1-(1/w))/(1/w)^2)*s_tau_i)%*%s_tau_i
  w_cov   <- chol2inv(chol(w_cov_aux))
  verossi <- -0.5*t(s_tau)%*%w_cov%*%s_tau - 0.5 * log(det(2*pi*w_cov_aux))
  
  funcao  <- priori+verossi
  
  return(funcao)
}
