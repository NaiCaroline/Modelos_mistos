// GLMM POISSON

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
     real theta[n]; // Media

     for (i in 1:n)
     {
          theta[i] = exp(eta[i] + b[id[i]]);
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
     

     y ~ poisson(theta); // maximaverossimilhanca
}
generated quantities
{
  vector[n] log_lik;
  for(i in 1:n)
  {
    log_lik[i] = poisson_lpmf(y[i] | theta[i]);
  }
}