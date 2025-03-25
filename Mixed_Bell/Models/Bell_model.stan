// GLMM BELL COM LOGVEROSSIMILHANCA INCOMPLETA

functions
{
    real lambertW(real x)
    {
      real y;
      real w;
      y= sqrt( 1 + exp(1) *x);
      w=-1 + 2.036 * log( (1 + 1.14956131 * y)/(1 + 0.45495740*log(1+y)) );
      w = (w/(1+w)) * ( 1 + log(x/w) );
      w = (w/(1+w)) * ( 1 + log(x/w) );
      w = (w/(1+w)) * ( 1 + log(x/w) );
      return(w);
    }

    real bellnumber(int n)
    {
      if(n < 2){ return(1); }
      else
      {
        int k;
        vector[n] B;
        vector[n] Bneu;
        B[1] = 1;
        for (i in 1:(n - 1))
        {
          k = i;
          Bneu[1] = B[i];
          for (j in 2:(i + 1))
          {
            Bneu[j] = B[j - 1] + Bneu[j - 1];
          }
          for(j in 1:n)
          {
            B[j] = Bneu[j];
          }
        }
        return(Bneu[k + 1]);
      }
    }

    real bell_lpmf(int x, real theta)
    {
      real Bx;
      real lprob;
      Bx = bellnumber(x);
      lprob = x*log(theta) - exp(theta) + 1 + log(Bx) - lgamma(x+1);
      return lprob;
    }

    real loglik_bell(int[] x, real[] theta)
    {
      real lprob[num_elements(x)];
      for(i in 1:num_elements(x))
      {
        lprob[i] = x[i]*log(theta[i]) - exp(theta[i]);
      }
      return sum(lprob);
    }
}
data
{
     int<lower=1> n; // observacoes
     int y[n]; // variavel resposta

     int<lower=1> p; // nº de betas
     matrix[n, p] X; // Matriz do modelo (Covariaveis)
     
     int<lower = 0> k; // nº de obs efeito aleatorio 
     int<lower=0, upper=k> id[n]; // efeito aleatorio
}
parameters
{
     vector[p] beta; // parametros fixos
     real<lower = 0> taub; // precisao
	vector[k] b; // parametro aleatorio
}
transformed parameters
{
     vector[n] eta = X*beta; // preditora linear
	real<lower = 0> sigmab = sqrt(1/taub); // DP efeito aleatorio
     real mu[n];
     real theta[n]; // Media
     
     for (i in 1:n)
     {
          mu[i] = exp(eta[i] + b[id[i]]);
          theta[i] = lambertW(mu[i]);
     }
}
model
{
     for(i in 1:p)
     {
         beta[i] ~ normal(0, 100); // priori efeitos fixos
     }

     taub ~ gamma(1, 0.0001);
     b ~ normal(0, sigmab); // priori efeito aleatorio


     target += loglik_bell(y, theta); // Loverossimilhanca incompleta
}
generated quantities
{
  vector[n] log_lik;
  for(i in 1:n)
  {
    log_lik[i] = bell_lpmf(y[i] | theta[i]);
  }
}