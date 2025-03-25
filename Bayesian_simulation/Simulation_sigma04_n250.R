#############################################################################
#                  Mixed Bell model simulation study - stan                 #
#                                                                           #
# Code Autor: Naiara Santos                                                 #
#                                                                           #
# Scenario: X1 ~ Normal(0,1)                                                #
#           b ~ Normal(0,sigma)                                             #
#           beta0 = 1.5                                                     #
#           beta1 = 0.5                                                     #
#           sigma = 0.4                                                     #
#                                                                           #
#############################################################################

# Remove objects and clear console
rm(list = ls()); cat('\014')

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
{library(rstan); library(rio); library(TeachingDemos)}

# Import criteria functions
source("../utils/Bayesian_criteria.R")

sigma_b = "04"

#############################################################################
# Global simulation parameters #
################################ 

# Number of replicas
R = 1000

# Sample size
p = 3   
k = 250
n = p * k

# Data structures to store the information generated in the replicas
Matriz_Pars_B = matrix(nrow = R, ncol = 3)
Matriz_Pars_P = matrix(nrow = R, ncol = 3)
Matriz_crit_B  = matrix(nrow = R, ncol = 6)
Matriz_crit_P  = matrix(nrow = R, ncol = 6)
Cont_Cobertura_B = rep(0, 3)
Cont_Cobertura_P = rep(0, 3)
Rhat_B = matrix(nrow = R, ncol = 3)
Rhat_P = matrix(nrow = R, ncol = 3)

## STARTING REPLICAS ##
for(r in 1:R)
{
  ## Generating data ##
  X1 = rep(rnorm(k), each = p)
  indiv = rep(1:k, each = p)
  beta0 = 1.5
  beta1 = 0.5
  sigmab = 0.4
  Pars = c(beta0, beta1, sigmab)
  
  b = rnorm(k, mean = 0, sd = sigmab)
  Y = c()
  aux1 = 1
  aux2 = p
  
  for (w in 1:k) {
    Y[aux1:aux2] = rBELL(p, LambertW::W(exp(beta0 + beta1 * X1[aux1:aux2] + b[w])))
    aux1 = aux1 + p
    aux2 = aux2 + p
  }
  
  dados = data.frame(indiv = indiv, Y = Y, X1 = X1)
  
  ## Estimation process ##
  print(paste("--------------- Starting replica", r, "---------------"))
  
  modMat = model.matrix(~ X1, data = dados)
  data_stan = list(n = length(dados$Y), y = dados$Y, p = dim(modMat)[2], 
                   X = modMat, k = length(unique(dados$indiv)),
                   id = dados$indiv)
  
  fit_Bell <- stan(file = "../Models/Bell_model.stan", data = data_stan, 
                   iter = 20000, chains = 3, thin = 10, 
                   warmup = 5000, seed = 10000003)
  
  fit_Pois <- stan(file = "../Models/Pois_model.stan", data = data_stan, 
                   iter = 20000, chains = 3, thin = 10, 
                   warmup = 5000, seed = 10000003)
  
  
  ####################################
  # Processing and storing results
  
  ##### BELL #####
  beta1_est_B = rowMeans(extract(fit_Bell, "beta[1]", inc_warmup = F, permuted = F))
  beta2_est_B = rowMeans(extract(fit_Bell, "beta[2]", inc_warmup = F, permuted = F))
  sigmab_est_B = rowMeans(extract(fit_Bell, "sigmab", inc_warmup = F, permuted = F))
  
  matriz_amostras_B = data.frame(beta1_est_B, beta2_est_B, sigmab_est_B)
  Matriz_Pars_B[r, ] = colMeans(matriz_amostras_B)
  Matriz_crit_B[r, ] = Criterios(dados$Y, fit_Bell, 3, "bell")
  
  Cont_Cobertura_B = Cont_Cobertura_B + (Pars >= summary(fit_Bell)$summary[c("beta[1]", "beta[2]", "sigmab"), "2.5%"] & 
                                           Pars <= summary(fit_Bell)$summary[c("beta[1]", "beta[2]", "sigmab"), "97.5%"])
  Rhat_B[r, ] = as.numeric(summary(fit_Bell)$summary[c("beta[1]", "beta[2]", "sigmab"), "Rhat"])
  
  ##### POISSON #####
  beta1_est_P = rowMeans(extract(fit_Pois, "beta[1]", inc_warmup = F, permuted = F))
  beta2_est_P = rowMeans(extract(fit_Pois, "beta[2]", inc_warmup = F, permuted = F))
  sigmab_est_P = rowMeans(extract(fit_Pois, "sigmab", inc_warmup = F, permuted = F))
  
  matriz_amostras_P = data.frame(beta1_est_P, beta2_est_P, sigmab_est_P)
  Matriz_Pars_P[r, ] = colMeans(matriz_amostras_P)
  Matriz_crit_P[r, ] = Criterios(dados$Y, fit_Pois, 3, "poisson")
  
  Cont_Cobertura_P = Cont_Cobertura_P + (Pars >= summary(fit_Pois)$summary[c("beta[1]", "beta[2]", "sigmab"), "2.5%"] & 
                                           Pars <= summary(fit_Pois)$summary[c("beta[1]", "beta[2]", "sigmab"), "97.5%"])
  Rhat_P[r, ] = as.numeric(summary(fit_Pois)$summary[c("beta[1]", "beta[2]", "sigmab"), "Rhat"])
  
  
  ###################################
  # Saving progress for each replica
  
  saveRDS(list(Matriz_Pars_B, Matriz_Pars_P, Matriz_crit_B, Matriz_crit_P, Cont_Cobertura_B, Cont_Cobertura_P, Rhat_B, Rhat_P), 
          file = paste0("Resultados/Simulation_sigma", sigma_b, "_n", k, "_replica_", ".rds"))
  
  
  print(paste("--------------- Replica", r, "completed ---------------"))
}


print("------- SIMULATION FINISHED!! -------")





































