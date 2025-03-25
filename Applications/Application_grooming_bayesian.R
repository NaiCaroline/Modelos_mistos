#############################################################################
#                        APPLICATION: GROOMING DATA                         #
#                           Bayesian approach                               #
#                                                                           #
# Code Author: Naiara Santos                                                #
#############################################################################

# Remove objects and clear the console
rm(list = ls()); cat('\014')

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
{library(rstan); library(MASS); library(TeachingDemos); 
  library(loo); library(reshape2)}

# Import criteria function
source("../utils/Bayesian_criteria.R")


#############################################################################
#############################################################################
# Data processing


# Dataset
dados_amplo <- data.frame(
  Tempo = rep(c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60), 2),
  Tratamento = c(rep("Saline", 12), rep("Oxytocin", 12)),
  Rato1 = c(0, 9, 6, 0, 9, 11, 16, 0, 0, 0, 0, 0, 
            9, 16, 16, 13, 17, 18, 18, 15, 13, 11, 8, 5),
  Rato2 = c(3, 3, 10, 0, 0, 0, 10, 0, 0, 0, 0, 0, 
            0, 13, 9, 0, 1, 5, 6, 2, 3, 0, 0, 0),
  Rato3 = c(2, 9, 8, 0, 6, 0, 4, 10, 0, 4, 12, 2, 
            15, 19, 18, 14, 9, 13, 1, 16, 15, 4, 13, 0),
  Rato4 = c(2, 6, 8, 5, 2, 0, 0, 2, 0, 0, 0, 0,
            0, 0, 2, 11, 14, 16, 7, 11, 0, 0, 7, 0),
  Rato5 = c(4, 6, 8, 5, 14, 4, 2, 5, 1, 0, 4, 0,
            2, 7, 17, 18, 19, 6, 16, 13, 7, 8, 9, 0),
  Rato6 = c(0, 0, 7, 1, 4, 0, 0, 0, 0, 0, 0, 0,
            2, 0, 3, 10, 11, 14, 9, 0, 1, 8, 3, 7),
  Rato7 = c(5, 9, 8, 7, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 8, 11, 14, 19, 15, 14, 9, 19, 8, 7, 1),
  Rato8 = c(8, 0, 0, 6, 1, 0, 0, 0, 0, 0, 0, 0,
            5, 6, 0, 1, 9, 11, 10, 8, 1, 10, 0, 7)
)

# Adjusting the data
dados <- melt(dados_amplo, id.vars = c("Tempo", "Tratamento"), 
                   variable.name = "ID", value.name = "y")
dados$ID <- as.integer(gsub("Rato", "", dados$ID))


head(dados)
str(dados)


#############################################################################
#############################################################################
# Model fit 
#############################################################################

# Stan data
modMat = model.matrix(~ Tratamento + log(Tempo), data = dados)
head(modMat)

dat_stan = list(n = nrow(modMat), 
                p = ncol(modMat),
                y = dados$y, 
                X = modMat, 
                k = length(unique(dados$ID)), 
                efa = dados$ID)

# Bell model
fit_bell = stan(file = "../Models/Bell.stan", data = dat_stan, chains = 3, 
                iter = 60000, warmup = 20000, thin = 10, seed = 10000003)

# Poisson model
fit_pois = stan(file = "../Models/Poisson.stan", data = dat_stan, chains = 3, 
                iter = 60000, warmup = 20000, thin = 10, seed = 10000003)


#############################################################################
#############################################################################
# Criterias

Crit_B = round(CritWL(fit_bell), 2)
Crit_P = round(CritWL(fit_pois), 2)

Res_Crit = cbind(Crit_B, Crit_P)
colnames(Res_Crit) = c("Bell", "Poisson")
rownames(Res_Crit) = c("WAIC", "LOOIC")

Res_Crit


#############################################################################
#############################################################################













