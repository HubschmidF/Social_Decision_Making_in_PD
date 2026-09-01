#==============================================================================# 
# SocialDBS: Posterior Parameters and parameter inference (SM)
#==============================================================================#

# FH

#==============================================================================#

#===============#
# Path, packages and data
#===============#

# Software packages
suppressMessages(suppressWarnings(library(rstudioapi)))  # path management
suppressMessages(suppressWarnings(library(HDInterval)))  # functions for highest density intervals
suppressMessages(suppressWarnings(library(ggplot2)))     # plots
suppressMessages(suppressWarnings(library(ggpubr)))      # figure generation
suppressMessages(suppressWarnings(library(readxl)))      # to read the data
suppressMessages(suppressWarnings(library(ggridges)))    # annotating plots 
suppressMessages(suppressWarnings(library(tibble)))      # data manipulation
suppressMessages(suppressWarnings(library(grid)))        # plot annotation
suppressMessages(suppressWarnings(library(plyr)))        # summarizing functions

# Path
analysis_path = getActiveDocumentContext()$path
setwd(dirname(analysis_path))
analysis_dir <- getwd()
  
# Data files --
# HC
data_healthy <- file.path(analysis_dir, "SM_posterior_draws_healthy.xlsx")
data_HC <- read_xlsx(data_healthy)
DF_HC <- as.data.frame(data_HC)
# DBS ON 
data_PATON <- file.path(analysis_dir, "SM_posterior_draws_DBS_ON.xlsx")
data_ON <- read_xlsx(data_PATON)
DF_ON <- as.data.frame(data_ON)
# DBS OFF 
data_PATOFF <- file.path(analysis_dir, "SM_posterior_draws_DBS_OFF.xlsx")
data_OFF <- read_xlsx(data_PATOFF)
DF_OFF <- as.data.frame(data_OFF)
  
#===============#
# Plots hyper parameters depending on groups
#===============#
  
# Utility
n_draws <- 8000
group_HC <- rep(1,n_draws)
group_ON <- rep(2,n_draws)
group_OFF <- rep(3,n_draws)
#
Group <- c(group_HC, group_ON, group_OFF)
Group <- as.factor(Group)
levels(Group) <- c('HC', 'ON', 'OFF')
  
# Annotation / significance
sig <- grobTree(textGrob("*", x = 0.9 , y = 0.7 , gp=gpar(col = "black", fontsize = 50)))
  
# inverse temperature (beta)
betaHC <- DF_HC$hyper_beta
betaON <- DF_ON$hyper_beta
betaOFF <- DF_OFF$hyper_beta
#
betas <- c(betaHC, betaON, betaOFF)
  
BETA <- data.frame(Group, betas)
mu1 <- ddply(BETA, "Group", summarise, grp.mean=mean(betas))
  
pl1 <- ggplot(BETA, aes(x = betas, color = Group, fill = Group)) +
              geom_density(alpha = 0.2) +
              geom_vline(data=mu1, aes(xintercept=grp.mean, color=Group), linetype = 'dashed') +
              labs(y = 'Density', x = 'Inverse temperature β') +
              scale_fill_manual(values=c("grey", "orange1", "purple1")) +
              scale_color_manual(values=c("grey", "orange1", "purple1")) +
              theme_pubr() +
              theme(legend.position='right')
pl1
  
# Discount
k_HC <- DF_HC$log_k
k_ON <- DF_ON$log_k
k_OFF <- DF_OFF$log_k
#
discounts <- c(k_HC, k_ON, k_OFF)
  
DISCOUNT <- data.frame(Group, discounts)
mu2 <- ddply(DISCOUNT, "Group", summarise, grp.mean=mean(discounts))
  
pl2 <- ggplot(DISCOUNT, aes(x = discounts, color = Group, fill = Group)) +
              geom_density(alpha = 0.2) +
              geom_vline(data=mu2, aes(xintercept=grp.mean, color=Group), linetype = 'dashed') +
              labs(y = 'Density', x = 'Discount log(k)') +
              scale_fill_manual(values=c("grey", "orange1", "purple1")) +
              scale_color_manual(values=c("grey", "orange1", "purple1")) +
              theme_pubr() +
              theme(legend.position='right')
pl2
  
#===============#
# Difference densities 
#===============#
  
# Inv Temp (difference distribution) ---
  
beta_HC_ON <- DF_HC$hyper_beta - DF_ON$hyper_beta
beta_HC_OFF <- DF_HC$hyper_beta - DF_OFF$hyper_beta
beta_OFF_ON <-  DF_OFF$hyper_beta - DF_ON$hyper_beta
# Utility
DF_beta_HC_ON <- tibble(beta_HC_ON)
DF_beta_HC_OFF <- tibble(beta_HC_OFF)
DF_beta_OFF_ON <- tibble(beta_OFF_ON)
  
# 95% HDI
# comparison 1
mean(beta_HC_ON)                  # 0.02405311
hdi(beta_HC_ON, credMass = .95)   # 0.005235014 0.044185076  ***
# comparison 2
mean(beta_HC_OFF)                 # 0.02843222
hdi(beta_HC_OFF, credMass = .95)  # 0.01028485 0.04759475 ***
# comparison 3
mean(beta_OFF_ON)                 # -0.004379116
hdi(beta_OFF_ON, credMass = .95)  # -0.018332267  0.008745081  
  
# corresponding plots ---
  
# comp 1
pl1A <- ggplot(DF_beta_HC_ON, aes(x = beta_HC_ON, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - ON") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr() + annotation_custom(sig)
pl1A
  
# comp 2
pl1B <- ggplot(DF_beta_HC_OFF, aes(x = beta_HC_OFF, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - OFF") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr() + annotation_custom(sig)
pl1B
  
# comp 3
pl1C <- ggplot(DF_beta_OFF_ON , aes(x = beta_OFF_ON, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "OFF - ON") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr()
pl1C
  
# Log(k) (difference distribution) ---
  
k_HC_ON <- DF_HC$log_k - DF_ON$log_k
k_HC_OFF <- DF_HC$log_k - DF_OFF$log_k
k_ON_OFF <-  DF_ON$log_k - DF_OFF$log_k
# Utility
DF_k_HC_ON <- tibble(k_HC_ON)
DF_k_HC_OFF <- tibble(k_HC_OFF)
DF_k_ON_OFF <- tibble(k_ON_OFF)

# comparison 1
mean(k_HC_ON)                  # -0.1504336
hdi(k_HC_ON, credMass = .95)   # -0.9327508  0.6915514 
# comparison 2
mean(k_HC_OFF)                 # -0.1303948
hdi(k_HC_OFF, credMass = .95)  # -0.9130716  0.6699845 
# comparison 3
mean(k_ON_OFF)                 # 0.02003888
hdi(k_ON_OFF, credMass = .95)  # -0.7466482  0.7839538 

# corresponding plots ---
  
# comp 1
pl2A <- ggplot(DF_k_HC_ON, aes(x = k_HC_ON, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - ON") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr()
pl2A
  
# comp 2
pl2B <- ggplot(DF_k_HC_OFF, aes(x = k_HC_OFF, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - OFF") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr()
pl2B
  
# comp 3
pl2C <- ggplot(DF_k_ON_OFF, aes(x = k_ON_OFF, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "ON - OFF") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr()
pl2C
  
# Diff densities for supplementary
  
Figure_it <- ggarrange(pl1A,  pl1B,  pl1C, labels = c('','',''), legend = "top", common.legend = TRUE ,nrow = 1, ncol = 3)
Figure_ks <- ggarrange(pl2A,  pl2B,  pl2C, labels = c('','',''), legend = "top", common.legend = TRUE ,nrow = 1, ncol = 3)
Figure_it
Figure_ks
  
#==============================================================================#
#------------------------------------------------------------------------------#
#==============================================================================#