# Full conditional for beta
atualizarBETA<-function(b,B,x,w,sigma,delta2,theta,v,dados){
  B.inv  <- chol2inv(chol(B))
  covar  <- chol2inv(chol(B.inv+(1/(delta2*sigma)*(t((w/v)*x)%*%x))))
  media  <- covar%*%((B.inv%*%b)+(1/(delta2*sigma))*(t((w/v)*x)%*%(dados-theta*v)))
  beta   <- rmvnorm(1,media,covar)
  return(beta)
}

# Full conditional for the latent variable
atualizarV<-function(dados,x,w,beta,delta2,theta,sigma,N){
  p1 <- 0.5
  p2 <- w*(dados-x%*%beta)^2/(delta2*sigma)
  p3 <- w*(2/sigma + theta^2/(delta2*sigma))
  v  <- NULL
  for(i in 1:N){
    v[i]  <- rgig(1, chi=p2[i], psi=p3[i], lambda=p1)
  }
  return(v)
}

# Bayesian Quantile Regression - MCMC
bayesQR_weighted <- function(y,a_min,b_max,x,w,tau,n_mcmc,burnin_mcmc,thin_mcmc){
  n         <- length(y)
  numcov    <- ncol(x)
  resultado <- list()
  
  # Create auxiliary objects
  beta  <- matrix(NA, n_mcmc, numcov)
  # Set the initial values
  reg_lm     <- glm(y~-1+x,family="poisson")
  beta[1,]   <- reg_lm$coefficients
  v          <- rgamma(n,2,1)
  # Auxiliary constants
  delta2  <- 2/(tau*(1-tau))
  theta   <- (1-2*tau)/(tau*(1-tau))
  # MCMC
  for(k in 2:n_mcmc){
    u_jit      <- runif(n,0,1)
    y_star     <- NULL
    y_star     <- log((y+u_jit-tau-(a_min+0)))-log((b_max+1)-(y+u_jit-tau))
    beta[k,]   <- atualizarBETA(rep(0,numcov),diag(rep(1000,numcov)),x,w,1,delta2,theta,v,y_star)
    v          <- atualizarV(y_star,x,w,beta[k,],delta2,theta,1,n)
  }
  resultado[[1]]  <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  #resultado[[2]]  <- u_jit[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  
  return(resultado)
}
