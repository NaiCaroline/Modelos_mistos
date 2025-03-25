##############################################################################
#                       COMPARISON CRITERIA FUNCTION                         #
#                                                                            #
# Code Author: Naiara C. A. Santos                                           #
##############################################################################

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
{library(loo); library(pscl); library(extraDistr); 
library(numbers); library(LambertW)}


##############################################################################
### AUXILIARY FUNCTIONS NEEDED ###

# FUNCAO MASSA  

dBELL = function(x, theta, log = FALSE) 
{
  Bx <- c()
  for(i in 1:length(x)) 
  {
    Bx[i] <- numbers::bell(x[i])
  }
  
  lf <- x * log(theta) - exp(theta) + 1 + log(Bx) - lgamma(x + 1)
  
  if (log == TRUE) 
  {
    return(lf)
  }
  else 
  {
    return(exp(lf))
  }
}


# FUNCAO PROBABILIDADE  

pBELL = function(q, theta, lower.tail = TRUE, log.p = FALSE) 
{
  prob <- c()
  for(i in 1:length(q)) 
  {
    prob[i] <- sum(dBELL(0:q[i], theta, log = FALSE))
  }
  
  if(lower.tail == FALSE) 
  {
    if(log.p == FALSE) 
    {
      return(1 - prob)
    }
    else 
    {
      return(log(1 - prob))
    }
  }
  else 
  {
    if(log.p == FALSE) 
    {
      return(prob)
    }
    else 
    {
      return(log(prob))
    }
  }
}


# FUNCAO QUANTILICA  

qBELL = function(p, theta, log.p = FALSE) 
{
  n <- length(p)
  k <- rep(0, n)
  for(i in 1:n) 
  {
    q <- 0
    while(q <= p[i]) 
    {
      q <- pBELL(k[i], theta)
      k[i] <- k[i] + 1
    }
  }
  if(log.p == TRUE) 
  {
    return(log(k - 1))
  }
  else 
  {
    return(k - 1)
  }
}


# FUNCAO GERADORA  

rBELL = function(n, theta) 
{
  lambda <- exp(theta) - 1
  x <- c()
  for(i in 1:n) 
  {
    N <- stats::rpois(1, lambda)
    x[i] <- sum(extraDistr::rtpois(n = N, lambda = theta, 
                                   a = 0, b = Inf))
  }
  return(x)
}


##############################################################################
### COMPARISON CRITERIA FUNCTION ###

Criterios = function(y, fit.stan, np, dist = "poisson")
{
  Dbar = mean(-2*rowSums(extract_log_lik(fit.stan, 
                         parameter_name = "log_lik")))

  theta = get_posterior_mean(fit.stan, pars= c("theta"))[,4]
  
  if(dist == "poisson")
  {
    Dhat = -2*sum(dpois(y, lambda = theta, log = T))
  }
  
  if(dist == "nbinom")
  {
    phi = mean(rstan::extract(fit.stan, "phi", permuted = F))
    Dhat = -2*sum(dnbinom(y, mu = theta, size = 1/phi, log = T))
  }
  
  if(dist == "bell")
  {
    Dhat = -2*sum(log(dBELL(y, theta)))
  }
  
  pD   = Dbar - Dhat
  
  # DEVIANCE INFORMATION CRITERION (DIC) 
  dic  = Dbar + pD
  
  # empirical Akaike information criterion (EAIC)  
  eaic = Dbar + 2*np
  
  # empirical Bayesian information criterion (EBIC)
  ebic = Dbar + np*log(length(y))
  
  # information criterion (IC)
  ic = Dbar + 2*pD
  
  # WATANABE-AKAIKE INFORMATION CRITERION (WAIC)
  waic = loo::waic(extract(fit.stan)$log_lik)$estimates[3,1]
  
  # LEAVE-ONE-OUT INFORMATION CRITERION (LOO, LOOIC)
  looic = loo::loo(extract(fit.stan)$log_lik)$estimates[3,1]
  
  saida = c(dic, waic, looic, eaic, ebic, ic)
  return(saida)
}


##############################################################################
##############################################################################

# Criteria Function 
CritWL = function(fit.stan)
{ 
  # WATANABE-AKAIKE INFORMATION CRITERION (WAIC)
  waic = loo::waic(extract(fit.stan)$log_lik)$estimates[3,1]
  
  # LEAVE-ONE-OUT INFORMATION CRITERION (LOO, LOOIC)
  looic = loo::loo(extract(fit.stan)$log_lik)$estimates[3,1]
  
  saida = c(waic, looic)
  return(saida)
}











