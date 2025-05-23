# Full conditional for the coefficients
condicionalBETA_MH_bio <- function(beta,b,B,dados,x,w,tau){
  n       <- length(dados)
  
  B_inv   <- chol2inv(chol(B))
  priori  <- -0.5*(t(beta-b)%*%B_inv%*%(beta-b))
  
  w_un    <- sum(data_obj[[6]])/n*w
  s_tau   <- t(w_un*x)%*%(tau-((dados-x%*%beta)<0))
  s_tau_i <- x*as.vector(tau-((dados-x%*%beta)<0))
  w_cov_aux<-t(((1-(1/w_un))/(1/w_un)^2)*s_tau_i)%*%s_tau_i
  w_cov   <- chol2inv(chol(w_cov_aux))
  verossi <- -0.5*t(s_tau)%*%w_cov%*%s_tau - 0.5 * log(det(2*pi*w_cov_aux))
  
  funcao  <- priori+verossi
  
  return(funcao)
}

# Metropolis-Hasting for the coefficients
atualizarBETA_MH_bio <- function(b,B,dados,x,beta,w,tau,ct,k,c0,c1){
  valoratual    <- beta
  sigma_MH      <- (tau*(1-tau))*chol2inv(chol((1/length(dados))*t(w^2*x)%*%x))
  valorproposto <- as.vector(rmvnorm(1, mean=valoratual, sigma=ct*sigma_MH))
  candidato     <- exp(condicionalBETA_MH_bio(valorproposto,b,B,dados,x,w,tau)-
                       condicionalBETA_MH_bio(valoratual,b,B,dados,x,w,tau))
  
  chanceaceitar <- min(1,candidato)
  if(runif(1)<chanceaceitar){
    BETAfinal   <- valorproposto
    aceita      <- 1
  } else{
    BETAfinal   <- valoratual
    aceita      <- 0
  }
  log_ct        <- log(ct)+c0*(1/k^c1)*(chanceaceitar-0.234)
  
  return(c(BETAfinal,aceita,exp(log_ct)))
}


# Bayesian Quantile Regression - Approximate Likelihood
bayesQRSL_weighted_bio <- function(y,x,w,tau,n_mcmc,burnin_mcmc,thin_mcmc,cte,c0,c1){
  n         <- length(y)
  numcov    <- ncol(x)
  resultado <- list()
  
  # Create auxiliary objects
  beta  <- matrix(NA, n_mcmc, numcov)
  # Set the initial values
  reg_lm     <- lm(y~-1+x)
  beta[1,]   <- reg_lm$coefficients
  contador   <- NULL
  const      <- NULL
  contador[1]<- 0
  const[1]   <- cte
  # MCMC
  for(k in 2:n_mcmc){
    beta_aux   <- atualizarBETA_MH_bio(rep(0,numcov),diag(rep(1000,numcov)),y,x,beta[k-1,],w,tau,const[k-1],k,c0,c1)
    beta[k,]   <- beta_aux[1:numcov] ; contador[k] <- beta_aux[(numcov+1)] ; const[k] <- beta_aux[(numcov+2)] 
  }
  resultado[[1]]  <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  
  return(resultado)
}
