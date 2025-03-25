#############################################################################
#                        APPLICATION: EPILEPSY DATA                         #
#                           Bayesian approach                               #
#                                                                           #
# Code Author: Naiara Santos                                                #
#############################################################################

# Remove objects and clear the console
rm(list = ls()); cat('\014')

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
{library(rstan); library(MASS); library(TeachingDemos); library(loo)}

# Import criteria function
source("../utils/Bayesian_criteria.R")


#############################################################################
#############################################################################
# Data reading

data(epil)  
dados = transform(epil, visit = (2*epil$period-5)/10)
str(dados)


#############################################################################
#############################################################################
# Model fit 
#############################################################################

# Stan data
modMat = model.matrix(~ trt + lbase, data = dados)
head(modMat)

dat_stan = list(n = nrow(modMat), 
                p = ncol(modMat),
                y = dados$y, 
                X = modMat, 
                k = length(unique(dados$subject)), 
                efa = dados$subject)

# Bell model
fit.bell = stan(file = "../Models/Bell.stan", data = dat_stan, chains = 3, 
                iter = 60000, warmup = 20000, thin = 10, seed = 10000003)

# Poisson model
fit.pois = stan(file = "../Models/Poisson.stan", data = dat_stan, chains = 3, 
                iter = 60000, warmup = 20000, thin = 10, seed = 10000003)


#############################################################################
#############################################################################
# Criterias

Crit.B = round(CritWL(fit.bell), 2)
Crit.P = round(CritWL(fit.pois), 2)

Res_Crit = cbind(Crit.B, Crit.P)
colnames(Res_Crit) = c("Bell", "Poisson")
rownames(Res_Crit) = c("WAIC", "LOOIC")

Res_Crit


#############################################################################
#############################################################################















