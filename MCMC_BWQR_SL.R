#' Bayesian Weighted Quantile Regression - Score-based Likelihood
#'
#' This function estimates a Bayesian quantile regression model for complex 
#' survey data under informative sampling, where a score-based likelihood
#' is assumed.
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

#Packages
require(mvtnorm)

# MCMC
bayesQRSL_weighted <- function(y,x,w,tau,n_mcmc,burnin_mcmc,thin_mcmc){
  n         <- length(y)
  numcov    <- ncol(x)
  resultado <- list()
  # Normalize weights
  w <- n*(w/sum(w))
  # Create auxiliary objects
  beta  <- matrix(NA, n_mcmc, numcov)
  # Set the initial values
  reg_lm     <- lm(y~-1+x)
  beta[1,]   <- reg_lm$coefficients
  const      <- NULL
  const[1]   <- 1
  # MCMC
  for(k in 2:n_mcmc){
    beta_aux   <- atualizarBETA_MH(rep(0,numcov),diag(rep(1000,numcov)),y,x,beta[k-1,],w,tau,const[k-1],k)
    beta[k,]   <- beta_aux[1:numcov] ; const[k] <- beta_aux[(numcov+1)] 
  }
  resultado[[1]]  <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  
  return(resultado)
}

# Metropolis-Hastings for the coefficients
atualizarBETA_MH <- function(b,B,dados,x,beta,w,tau,ct,k){
  valoratual    <- beta
  sigma_MH      <- (tau*(1-tau))*chol2inv(chol((1/length(dados))*t(w^2*x)%*%x))
  valorproposto <- as.vector(rmvnorm(1, mean=valoratual, sigma=ct*sigma_MH))
  candidato     <- exp(condicionalBETA_MH(valorproposto,b,B,dados,x,w,tau)-
                         condicionalBETA_MH(valoratual,b,B,dados,x,w,tau))
  
  chanceaceitar <- min(1,candidato)
  if(runif(1)<chanceaceitar){
    BETAfinal   <- valorproposto
  } else{
    BETAfinal   <- valoratual
  }
  log_ct    <- log(ct)+(1/k^0.8)*(chanceaceitar-0.234)
  
  return(c(BETAfinal,exp(log_ct)))
}

# Full conditional for the coefficients
condicionalBETA_MH <- function(beta,b,B,dados,x,w,tau){
  n       <- length(dados)
  
  B_inv   <- chol2inv(chol(B))
  priori  <- -0.5*(t(beta-b)%*%B_inv%*%(beta-b))
  
  s_tau   <- t(w*x)%*%(tau-((dados-x%*%beta)<0))
  w_cov   <- (n/(tau*(1-tau)))*chol2inv(chol(t(w^2*x)%*%x))
  verossi <- -(1/(2*n))*t(s_tau)%*%w_cov%*%s_tau
  
  funcao  <- priori+verossi
  
  return(funcao)
}
