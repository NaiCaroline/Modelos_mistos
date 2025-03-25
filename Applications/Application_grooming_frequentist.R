#############################################################################
#                       APPLICATION: GROOMING DATA                          #
#                           Frequentist approach                            #
#                                                                           #
# Code Author: Naiara Santos                                                #
#############################################################################

# Remove objects and clear the console
rm(list = ls()); cat('\014')

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
{
  library(gamlss); library(TeachingDemos); 
  library(ICglm); library(ggplot2);
  library(tidyr); library(dplyr); library(reshape2)
}

# Import model
source("../Models/Bell_model_gamlss.R")


#############################################################################
#############################################################################
# Data processing


# # Dataset
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

# Poisson
fit_pois = gamlss(y ~ Tratamento + log(Tempo) + re(random = ~1|ID),  
                  data = dados, family = PO) 

# Bell
fit_bell = gamlss(y ~ Tratamento + log(Tempo) + re(random = ~1|ID),  
                  data = dados, family = BELL) 


#############################################################################
#############################################################################
# Criterias

Crit_pois = round(c(fit_pois$aic, fit_pois$sbc, ICglm::HQIC(fit_pois),
                    ICglm::CAIC(fit_pois), ICglm::KICC(fit_pois)), 2)

Crit_bell = round(c(fit_bell$aic, fit_bell$sbc, ICglm::HQIC(fit_bell),
                    ICglm::CAIC(fit_bell), ICglm::KICC(fit_bell)), 2)


Res_Crit = cbind(Crit_pois, Crit_bell)
colnames(Res_Crit) = c("Poisson", "Bell")
rownames(Res_Crit) = c("AIC", "SBC", "HQIC", "CAIC", "KICC")

Res_Crit = as.data.frame(Res_Crit)
colnames_menor_valor <- apply(Res_Crit, 1, function(row) {
  names(Res_Crit)[which.min(row)]
})

Res_Crit$MelhorModelo <- colnames_menor_valor
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
  geom_hline(yintercept = -3, colour = "black", linetype = 2, linewidth = 0.9) + 
  geom_hline(yintercept = 3, colour = "black", linetype = 2, linewidth = 0.9) + 
  labs(title = "Bell model", x = "Sample", y = "Quantile residuals") +
  coord_cartesian(ylim = c(-4, 4.5)) +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  scale_x_continuous(limits = c(0, 193), expand = c(0.01, 0.1)) +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(-3,3, 1.5))
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
  geom_hline(yintercept = -3, colour = "black", linetype = 2, linewidth = 0.9) + 
  geom_hline(yintercept = 3, colour = "black", linetype = 2, linewidth = 0.9) + 
  labs(title = "Poisson model", x = "Sample", y = "Quantile residuals") +
  coord_cartesian(ylim = c(-5, 5)) +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  scale_x_continuous(limits = c(0, 193), expand = c(0.01, 0.1)) +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(-3,3, 1.5))
plot_quantilePO




