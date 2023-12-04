# Full conditional for the coefficients
condicionalBETA_MH <- function(beta,b,B,dados,x,w,tau){
  n       <- length(dados)
  
  B_inv   <- chol2inv(chol(B))
  priori  <- -0.5*(t(beta-b)%*%B_inv%*%(beta-b))
  
  verossi <- sum( - w * (dados-x%*%beta) * (tau-((dados-x%*%beta)<0)) )
  
  funcao  <- priori+verossi
  
  return(funcao)
}

# Metropolis-Hasting for the coefficients
atualizarBETA_MH <- function(b,B,dados,x,beta,w,tau,ct,k,t,c0,c1,contador){
  valoratual    <- beta
  #log_ct        <- ifelse(k%%t==1,
  #                        log(ct)+c0*(1/k^c1)*(sum(contador[(k-t):(k-1)])/t-0.234),
  #                        log(ct))
  ##sigma_MH      <- (length(dados)/(tau*(1-tau)))*chol2inv(chol(t(w^2*x)%*%x))
  sigma_MH      <- (tau*(1-tau))*chol2inv(chol((1/length(dados))*t(w^2*x)%*%x))
  #valorproposto <- as.vector(rmvnorm(1, mean=valoratual, sigma=exp(log_ct)*sigma_MH))
  valorproposto <- as.vector(rmvnorm(1, mean=valoratual, sigma=ct*sigma_MH))
  candidato     <- exp(condicionalBETA_MH(valorproposto,b,B,dados,x,w,tau)-
                         condicionalBETA_MH(valoratual,b,B,dados,x,w,tau))
  
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


# Bayesian Quantile Regression - MCMC
bayesQRPL_weighted <- function(y,a_min,b_max,x,w,tau,n_mcmc,burnin_mcmc,thin_mcmc,cte,t,c0,c1){
  n         <- length(y)
  numcov    <- ncol(x)
  resultado <- list()
  
  # Create auxiliary objects
  beta  <- matrix(NA, n_mcmc, numcov)
  # Set the initial values
  reg_lm     <- glm(y~-1+x,family="poisson")
  beta[1,]   <- reg_lm$coefficients
  v          <- rgamma(n,2,1)
  contador   <- NULL
  const      <- NULL
  contador[1]<- 0
  const[1]   <- cte
  # Auxiliary constants
  delta2  <- 2/(tau*(1-tau))
  theta   <- (1-2*tau)/(tau*(1-tau))
  # MCMC
  for(k in 2:n_mcmc){
    u_jit      <- runif(n,0,1)
    y_star     <- NULL
    y_star     <- log((y+u_jit-tau-(a_min+0)))-log((b_max+1)-(y+u_jit-tau))
    #beta[k,]   <- atualizarBETA(rep(0,numcov),diag(rep(1000,numcov)),x,w,1,delta2,theta,v,y_star)
    beta_aux   <- atualizarBETA_MH(rep(0,numcov),diag(rep(1000,numcov)),y_star,x,beta[k-1,],w,tau,const[k-1],k,t,c0,c1,contador)
    beta[k,]   <- beta_aux[1:numcov] ; contador[k] <- beta_aux[(numcov+1)] ; const[k] <- beta_aux[(numcov+2)] 
  }
  resultado[[1]]  <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  
  return(resultado)
}
