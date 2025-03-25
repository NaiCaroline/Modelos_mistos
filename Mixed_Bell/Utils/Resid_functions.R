##############################################################################
#                               RESIDUALS FUNCTION                           #
#                                                                            #
# Code Author: Naiara Santos                                                 #
##############################################################################

# Read Bell functions R code
source("../Utils/Bell_dpqr.R")

##############################################################################
##############################################################################

# Quantile residuals
Res_Quantile = function(y, fit, method = c('bayesian', 'frequentist'),
                        family = c('BELL', 'PO', 'NBI'))
{
  if(method == 'bayesian')
  {
    mu = get_posterior_mean(fit, pars = c("theta"))[,4]
    
    if(family == "NBI")
    {
      Phi = mean(rstan::extract(fit, "phi", permuted = F))
    }
  }
  
  if(method == 'frequentist')
  {
    mu = fitted(fit)
    
    if(family == "NBI")
    {
      Phi = exp(fit$sigma.coefficients)
    }
  }
  
  
  if(family == "PO")
  {
    q = c()
    for(i in 1:length(y))
    {
      a = ppois(y[i] - 1, lambda = mu[i])
      b = ppois(y[i], lambda = mu[i])
      u = runif(n = 1, min = a, max = b)
      
      q[i] = qnorm(u)
    }
  }
    
  if(family == "NBI")
  {
    q = c()
    for(i in 1:length(y))
    {
      a = pnbinom(y[i] - 1, size = 1/Phi, mu = mu[i])
      b = pnbinom(y[i], size = 1/Phi, mu = mu[i])
      u = runif(n = 1, min = a, max = b)
      
      q[i] = qnorm(u)
    }
  }
    
  if(family == "BELL")
  {
    q = c()
    for(i in 1:length(y))
    {
      a = pBELL(y[i] - 1, mu[i])
      b = pBELL(y[i], mu[i])
      u = runif(n = 1, min = a, max = b)
      
      q[i] = qnorm(u)
    }
  }
  
  return(q)
}


##############################################################################
##############################################################################


Res_Pearson = function(y, fit, method = c('bayesian', 'frequentist'),
                       family = c('BELL', 'PO', 'NBI'))
{
  if(method == 'bayesian')
  {
    mu <- get_posterior_mean(fit, pars = c("theta"))[,4] 
    r  <- y - mu
    
    if(family == "BELL")
    {
      V   <- mu * (1 + mu) * exp(mu)
      res <- r / sqrt(V)
    }
    
    if(family == "PO")
    {
      V   <- mu
      res <- r / sqrt(V)
    }
    
    if(family == "NBI")
    {
      Phi <- mean(rstan::extract(fit, "phi", permuted = F))
      V   <- mu + (1/Phi) * (mu^2)
      res <- r / sqrt(V)
    }
  }
  
  if(method == 'frequentist')
  {
    mu <- fitted(fit)
    r  <- y - mu
    
    if(family == "BELL")
    {
      V   <- mu * (1 + LambertW::W(mu)) * exp(LambertW::W(mu))
      res <- r / sqrt(V)
    }
    if(family == "PO")
    {
      V   <- mu
      res <- r / sqrt(V)
    }
    
    if(family == "NBI")
    {
      Phi <- exp(fit$sigma.coefficients)
      V     <- mu + Phi * mu^2
      res   <- r / sqrt(V)
    }
  } 
  
  return(res)
}



























