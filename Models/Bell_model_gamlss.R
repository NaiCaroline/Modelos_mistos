#############################################################################
#            Inclusion of the Bell model in the gamlss package              #
#############################################################################


BELL <- function (mu.link = "log") 
{
  mstats <- checklink("mu.link", "Bell", substitute(mu.link),
                      c("inverse", "log", "sqrt", "identity")) 
  structure(
    list(family = c("BELL", "Bell"),
         parameters = list(mu = TRUE), # the mean
         nopar = 1, 
         type = "Discrete", 
         mu.link = as.character(substitute(mu.link)), 
         mu.linkfun = mstats$linkfun, 
         mu.linkinv = mstats$linkinv, 
         mu.dr = mstats$mu.eta, 
         dldm = function(y, mu) {
           nd <- numeric.deriv(dBELL(y, mu, log = TRUE), 
                               "mu", delta = 0.001)
           dldv <- as.vector(attr(nd, "gradient"))
           dldv
         },
         d2ldm2 =function(y, mu) {
           nd <- numeric.deriv(dBELL(y, mu,  log = TRUE), 
                               "mu", delta = 0.001)
           dldv <- as.vector(attr(nd, "gradient"))
           d2ldv2 <- -dldv * dldv
           d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2, -1e-15)
           d2ldv2 },
         G.dev.incr = function(y, mu,...)
           -2 * dBELL(y, mu, log = TRUE),
         rqres = expression(rqres(pfun="pBELL", type="Discrete", ymin=0, 
                                  y=y, mu=mu)), 
         mu.initial = expression({mu <- (y + mean(y))/2 } ),
         mu.valid = function(mu) all(mu > 0), 
         y.valid = function(y)  all(y >= 0) 
    ),
    class = c("gamlss.family","family"))
}

#############################################################################
#############################################################################
# Functions dpqr of the Bell model

dBELL = function(x, mu = 1, log = FALSE)
{ 
  if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
  if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))  
  ly <- max(length(x), length(mu))
  x <- rep(x, length = ly)
  mu <- rep(mu, length = ly)
  Bx = c()
  for (i in 1:length(x)) 
  {
    Bx[i] = numbers::bell(x[i])
  }
  
  lf = x * log(LambertW::W(mu)) - exp(LambertW::W(mu)) + 1 + 
       log(Bx) - lgamma(x + 1)
  
  if(log == TRUE) 
  {
    return(lf)
  }
  else 
  {
    return(exp(lf))
  }
}

#------------------------------------------------------------------------------

pBELL = function(q, mu = 1, lower.tail = TRUE, log.p = FALSE) 
{
  if (any(mu <= 0)){stop(paste("mu must be greater than 0 ", "\n", ""))} 
  if (any(q < 0)){return(0)}
  
  prob = c()
  for(i in 1:length(q)) 
  {
    prob[i] = sum(dBELL(0:q[i], mu, log = FALSE))
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

#------------------------------------------------------------------------------

qBELL = function(p, mu = 1, log.p = FALSE) 
{
  if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
  if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", 
                                           "\n", ""))
  
  n = length(p)
  k = rep(0, n)
  for(i in 1:n) 
  {
    q = 0
    while(q <= p[i]) 
    {
      q = pBELL(k[i], mu)
      k[i] = k[i] + 1
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

#------------------------------------------------------------------------------

rBELL = function(n, mu = 1, log.p = FALSE) 
{
  if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
  if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))
  
  lambda = exp(LambertW::W(mu)) - 1
  x = c()
  for(i in 1:n) 
  {
    N = stats::rpois(1, lambda)
    x[i] = sum(extraDistr::rtpois(n = N, lambda = LambertW::W(mu), 
                                  a = 0, b = Inf))
  }
  return(x)
}

#------------------------------------------------------------------------------






























