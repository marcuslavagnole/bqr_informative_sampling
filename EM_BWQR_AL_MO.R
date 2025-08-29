#' Bayesian approach to multiple-output quantile regression analysis under 
#' informative sampling
#'
#' This function employs an EM algorithm for a Bayesian approach to 
#' multiple-output quantile regression analysis under informative sampling, 
#' where the "transformed" response variable is assumed to follow an asymmetric 
#' Laplace distribution.
#'
#' @param y Matrix - nxd - of responses.
#' @param x Matrix - nxp - of predictors (include the intercept)
#' @param w n-dimensional vector of survey weights.
#' @param u direction (d-dimensional vector).
#' @param gamma_u orthogonal basis (dx(d-1)).
#' @param tau Quantile of interest (value between 0 and 1),
#' @param mu0 hyperparameter of the mean of the multivariate normal 
#' distribution for the coefficients.
#' @param sigma0 hyperparameter of the covariance matrix of the multivariate 
#' normal distribution for the coefficients.
#' @param a0 hyperparameter of the inverse gamma distribution for the scale 
#' parameter of the asymmetric Laplace distribution.
#' @param b0 hyperparameter of the inverse gamma distribution for the scale 
#' parameter of the asymmetric Laplace distribution.
#' 
#' @return A list with the estimates of all parameters of interest.

bayesQR_weighted_EM <- function(y,x,w,u,gamma_u,tau,mu0,sigma0,a0,b0){
  p <- ifelse(is.null(ncol(x)),1,ncol(x))
  n <- dim(y)[1]
  resultado <- list()
  
  # Normalize weights
  w <- n*(w/sum(w))
  
  # Directional data
  y_u     <- as.vector(y%*%u)
  y_gamma <- as.vector(y%*%gamma_u)
  x_star  <- cbind(x,y_gamma)
    
  # Create auxiliary objects
  beta  <- NULL
  sigma <- NULL
  
  # Auxiliary constants
  delta2  <- 2/(tau*(1-tau))
  theta   <- (1-2*tau)/(tau*(1-tau))
  
  # Set the initial values
  reg_qr <- lm(y_u~-1+x_star)
  beta   <- rbind(beta,reg_qr$coefficients)
  sigma  <- c(sigma,1)
  
  # Stop criteria
  aux_criteria <- c(rep(1,p+1),1)
  criteria <- 10^(-5)
  i = 2
  
  while(sum(aux_criteria)>criteria){
    
    # Expectation
    a_star  <- as.vector((w*(y_u-x_star%*%beta[i-1,])^2)/(delta2*sigma[i-1]))
    b_star  <- w*(2/sigma[i-1] + theta^2/(delta2*sigma[i-1]))
    a_star[which(a_star==0)] <- 10^(-10)
    
    e_nu_inverse <- sqrt(b_star)/sqrt(a_star) 
    aux_obj      <- sqrt(a_star*b_star)
    e_nu         <- NULL
    for(j in 1:n){
      e_nu[j] <- 1/e_nu_inverse[j]*(besselK(aux_obj[j], 1.5)/besselK(aux_obj[j], 0.5))  
    }
    
    # Maximization
    cov_inv     <- (1/(delta2*sigma[i-1]))*(diag(w*e_nu_inverse))
    y_aux       <- y_u-theta/e_nu_inverse
    
    sigma0.inv  <- chol2inv(chol(sigma0))
    beta_aux    <- chol2inv(chol((t(x_star)%*%cov_inv%*%x_star+sigma0.inv)))%*%(t(x_star)%*%cov_inv%*%y_aux+sigma0.inv%*%mu0)
    beta        <- rbind(beta,as.vector(beta_aux))
    
    numerator   <- sum( (w*e_nu_inverse*(y_aux-x_star%*%beta[i,])^2)/(2*delta2) ) +
      sum( (w*theta^2*(e_nu-1/e_nu_inverse))/(2*delta2) ) +
      sum( w*e_nu ) +
      b0
    denominator <- (3*n+a0+1)/2
    sigma_aux   <- numerator/denominator
    sigma       <- c(sigma,sigma_aux)
    
    aux_criteria<-c(abs(beta[i,]-beta[i-1,]),abs(sigma[i]-sigma[i-1]))  
    i = i+1
  }
  resultado[[1]]  <- beta[i-1,]
  resultado[[2]]  <- sigma[i-1]
  
  return(resultado)
}
