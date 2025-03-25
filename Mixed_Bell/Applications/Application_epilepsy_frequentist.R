#############################################################################
#                        APPLICATION: EPILEPSY DATA                         #
#                           Frequentist approach                            #
#                                                                           #
# Code Author: Naiara Santos                                                #
#############################################################################

# Remove objects and clear the console
rm(list = ls()); cat('\014')

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
{library(gamlss); library(MASS); library(TeachingDemos); library(ICglm)}

# Import model
source("../Models/Bell_model_gamlss.R")


#############################################################################
#############################################################################
# Data reading

data(epil)  
dados = transform(epil, visit = (2*epil$period-5)/10, 
                  Subject = factor(epil$subject))
str(dados)


#############################################################################
#############################################################################
# Model fit 
#############################################################################

# Bell model
fit_bell = gamlss(y ~ trt + lbase + re(random = ~1|Subject), 
                  data = dados, family = BELL, gd.tol=Inf)

# Poisson model
fit_pois = gamlss(y ~ trt + lbase + re(random = ~1|Subject),  
                  data = dados, family = PO, gd.tol=Inf) 

####################

# Poisson inverse Gaussian model
fit_pig = gamlss(y ~ trt + lbase + re(random = ~1|Subject),  
                  data = dados, family = PIG, gd.tol=Inf) 


#############################################################################
#############################################################################
# Criterias

Crit_B = round(c(fit_bell$aic, fit_bell$sbc, ICglm::CAIC(fit_bell)), 2)
Crit_P = round(c(fit_pois$aic, fit_pois$sbc, ICglm::CAIC(fit_pois)), 2) 
Crit_Pig = round(c(fit_pig$aic, fit_pig$sbc, ICglm::CAIC(fit_pig)), 2) 


Res_Crit = cbind(Crit_P, Crit_Pig, Crit_B)
colnames(Res_Crit) = c("Poisson", "Poisson IG", "Bell")
rownames(Res_Crit) = c("AIC", "SBC", "CAIC")

Res_Crit


#############################################################################
#############################################################################
#                            RESIDUALS ANALYSIS                             #
#############################################################################
#############################################################################

# Read residuals functions R code
source("../Utils/Resid_functions.R")


#############################################################################
# Bell model

Qres = Res_Quantile(dados$y, fit_bell, method = 'frequentist',
                    family = 'BELL')

obs = seq(1:length(Qres))
ResBell = data.frame(obs, Qres) 
head(ResBell)

summary(ResBell)

###########################################

par(mar = c(4.5, 4.5, 1, 1))
plot_quantile = ggplot(data = ResBell, mapping = aes(x = obs, y = Qres)) + 
  geom_point() + 
  geom_hline(yintercept = -2, colour = "black", linetype = 2, linewidth = 0.9) + 
  geom_hline(yintercept = 2, colour = "black", linetype = 2, linewidth = 0.9) + 
  labs(title = "Bell model", x = "Sample", y = "Quantile residuals") +
  coord_cartesian(ylim = c(-3, 3.5)) +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  scale_x_continuous(limits = c(0, 240), expand = c(0.01, 0.1)) +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(-3,3.5, 1))
plot_quantile


#############################################################################
# Poisson model

QresPO = Res_Quantile(dados$y, fit_pois, method = 'frequentist',
                      family = 'PO')

obs = seq(1:length(QresPO))
ResPO = data.frame(obs, QresPO) 
head(ResPO)

summary(ResPO)

###########################################

par(mar = c(4.5, 4.5, 1, 1))
plot_quantilePO = ggplot(data = ResPO, mapping = aes(x = obs, y = QresPO)) + 
  geom_point() + 
  geom_hline(yintercept = -2, colour = "black", linetype = 2, linewidth = 0.9) + 
  geom_hline(yintercept = 2, colour = "black", linetype = 2, linewidth = 0.9) + 
  labs(title = "Poisson model", x = "Sample", 
       y = "Quantile residuals") +
  coord_cartesian(ylim = c(-4, 7)) +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  scale_x_continuous(limits = c(0, 240), expand = c(0.01, 0.1)) +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(-4,7, 1))
plot_quantilePO






























