#############################################################################
#                Mixed Bell model simulation study - gamlss                 #
#                                                                           #
# Code Autor: Naiara Santos                                                 #
#                                                                           #
# Scenario: X1 ~ Normal(0,1)                                                #
#           b ~ Normal(0,sigma)                                             #
#           beta0 = 1.5                                                     #
#           beta1 = 0.5                                                     #
#           sigma = 0.6                                                     #
#                                                                           #
#############################################################################

# Remove objects and clear console
rm(list = ls()); cat('\014')

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
{library(TeachingDemos); library(gamlss); library(ICglm); library(rio)}

# Import model
source("../Models/Bell_model_gamlss.R")


#############################################################################
# Global Simulation Parameters #
################################

# Number of replicas
R = 1000

# Lists where to store the results of each sample size
Res_amostra = Dados_crit = Dados_pred = list()

# Fixed quantities
p = 3   
amost = c(50, 100, 250, 500) 


for(i in 1:length(amost))
{
  # Sample size
  k = amost[i]
  n = p*k
  
  # Data structures to store the information generated in the replicates
  Matriz_Pars_B = matrix(nrow = R, ncol = 3)
  Matriz_Pars_P = matrix(nrow = R, ncol = 3)
  Matriz_crit_B  = matrix(nrow = R, ncol = 5)
  Matriz_crit_P  = matrix(nrow = R, ncol = 5)
  Matriz_pred_B = matrix(nrow = n, ncol = R)
  Matriz_pred_P = matrix(nrow = n, ncol = R)
  Matriz_time_B = matrix(nrow = R, ncol = 1)
  Matriz_time_P = matrix(nrow = R, ncol = 1)
  
  
  #### STARTING REPLICAS ####
  for(r in 1:R)
  {
    ## Generating data ##
    X1 = rep(rnorm(k), each = p)
    indiv = rep(1:k, each = p)
    beta0 = 1.5
    beta1 = 0.5
    sigmab = 0.6
    Pars = c(beta0, beta1, sigmab)
    
    b = rnorm(k, mean = 0, sd = sigmab)
    Y = c()
    aux1 = 1
    aux2 = p
    
    for(w in 1:k)
    {
      Y[aux1:aux2] = rBELL(p, exp(beta0 + beta1*X1[aux1:aux2] + b[w]))
      aux1 = aux1 + p
      aux2 = aux2 + p
    }
    
    dados = data.frame(indiv = as.factor(indiv), Y = Y, X1 = X1)
    
      
    ## Estimation process ##
    print(paste("Replica:", toString(r),
                "Sample size:", toString(n),
                collapse = " "))
    
      
    tempo_fit_Bell = system.time(
      fit_Bell <- gamlss(Y ~ X1 + random(indiv), data = dados, 
                         family = BELL, gd.tol = Inf)
    ) 
    
    
    tempo_fit_Pois = system.time(
      fit_Pois <- gamlss(Y ~ X1 + random(indiv), data = dados, 
                         family = PO, gd.tol = Inf)
    ) 
    
    
    ####################################
    # Processing and storing results
    
    
    ##### BELL #####
    beta0_est_B = coef(fit_Bell, "mu")[1]
    beta1_est_B = coef(fit_Bell, "mu")[2]
    sigmab_est_B = getSmo(fit_Bell)$sigb 
    
    
    Matriz_Pars_B[r, ] = c(beta0_est_B, beta1_est_B, sigmab_est_B)
    Matriz_time_B[r, 1] = tempo_fit_Bell[[3]]
    Matriz_crit_B[r, ] = round(c(fit_Bell$aic, fit_Bell$sbc, ICglm::HQIC(fit_Bell),
                                 ICglm::CAIC(fit_Bell), ICglm::KICC(fit_Bell)), 2)
    Matriz_pred_B[, r] = predict(fit_Bell, type = "response")
    
    
    ##### POISSON #####
    beta0_est_P = coef(fit_Pois, "mu")[1]
    beta1_est_P = coef(fit_Pois, "mu")[2]
    sigmab_est_P = getSmo(fit_Pois)$sigb 
    
    
    Matriz_Pars_P[r, ] = c(beta0_est_P, beta1_est_P, sigmab_est_P)
    Matriz_time_P[r, 1] = tempo_fit_Pois[[3]]
    Matriz_crit_P[r, ] = round(c(fit_Pois$aic, fit_Pois$sbc, ICglm::HQIC(fit_Pois),
                                 ICglm::CAIC(fit_Pois), ICglm::KICC(fit_Pois)), 2) 
    Matriz_pred_P[, r] = predict(fit_Pois, type = "response")
    
    
  } 
  
  print(paste("------- Replicas finished!! Sample", amost[i]), digits = 15)
  
  
  #############################################################################
  #############################################################################
  # Calculating simulation metrics

  print(paste("------- Starting calculation of metrics: Sample", 
              amost[i], "-------"))
  
  
  ##### BELL #####
  Estimativas_B = SD_B = Vies_B = EQM_B = RMSE_B = NULL
  for(j in 1:ncol(Matriz_Pars_B)) 
  {
    Estimativas_B[j] = mean(Matriz_Pars_B[, j])
    SD_B[j] = sd(Matriz_Pars_B[, j])
    Vies_B[j] = Pars[j] - Estimativas_B[j]
    EQM_B[j] = SD_B[j]^2 + Vies_B[j]^2
    RMSE_B[j] = sqrt(SD_B[j]^2 + Vies_B[j]^2)
  }
  
  Res_B = data.frame(Parametros = c('beta0', 'beta1', 'sigma_b'), 
                     Real = Pars, Estimativas = Estimativas_B, SD = SD_B, 
                     Vies = Vies_B, EQM = EQM_B, RMSE = RMSE_B, 
                     tempo = mean(Matriz_time_B))
  
      
  ##### POISSON #####
  Estimativas_P = SD_P = Vies_P = EQM_P = RMSE_P = NULL
  for(j in 1:ncol(Matriz_Pars_P)) 
  {
    Estimativas_P[j] = mean(Matriz_Pars_P[, j])
    SD_P[j] = sd(Matriz_Pars_P[, j])
    Vies_P[j] = Pars[j] - Estimativas_P[j]
    EQM_P[j] = SD_P[j]^2 + Vies_P[j]^2
    RMSE_P[j] = sqrt(SD_P[j]^2 + Vies_P[j]^2)
  }
  
  Res_P = data.frame(Parametros = c('beta0', 'beta1', 'sigma_b'), 
                     Real = Pars, Estimativas = Estimativas_P, SD = SD_P, 
                     Vies = Vies_P, EQM = EQM_P, RMSE = RMSE_P, 
                     tempo = mean(Matriz_time_P))
  
  
  Preditos = list(Matriz_pred_B, Matriz_pred_P)
  names(Preditos) = c("Bell", "Poisson")
  Dados_pred[[i]] = Preditos
  
  Res = list(Res_B, Res_P)
  names(Res) = c("Bell", "Poisson")
  Res_amostra[[i]] = Res
  
  Criterios = cbind.data.frame(Matriz_crit_B, Matriz_crit_P)
  colnames(Criterios) = rep(c("AIC", "BIC", "HQIC", "CAIC", "KICC"), 2)
  Dados_crit[[i]] = Criterios
    
} 


#############################################################################
# Saving all simulation results

names(Dados_pred) = as.character(amost)
names(Res_amostra) = as.character(amost)
names(Dados_crit) = as.character(amost)

rio::export(Dados_crit, 'Simu_CenarioI_Criterios.xlsx')
rio::export(Res_amostra, 'Simu_CenarioI_Resultados.xlsx')
rio::export(Dados_pred, 'Simu_CenarioI_Preditos.xlsx')

save.image("Simu_CenarioI.RData")

txtStart("ResultadosParametros_CenarioI.txt")

Res_amostra

txtStop()


print("------- SIMULATION FINISHED!! -------")




































      