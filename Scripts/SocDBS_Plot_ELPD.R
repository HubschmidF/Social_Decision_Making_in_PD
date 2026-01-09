#==============================================================================# 
# SocialDBS: ELPD differences plots
#==============================================================================#

# FH

# Description:
#  - Plots ELPD difference to best model in all groups
#  - Calculates sigma heuristic

#==============================================================================#

#===============#
# Utility
#===============#

suppressMessages(suppressWarnings(library(rstudioapi))) # path management
suppressMessages(suppressWarnings(library(dplyr)))      # data manipulation
suppressMessages(suppressWarnings(library(plyr)))       # data manipulation
suppressMessages(suppressWarnings(library(tidyr)))      # R backbone package
suppressMessages(suppressWarnings(library(car)))        # R backbone package
suppressMessages(suppressWarnings(library(ggplot2)))    # plotting
suppressMessages(suppressWarnings(library(ggpubr)))     # ggplot add on
suppressMessages(suppressWarnings(library(ggridges)))   # ggplot add on
suppressMessages(suppressWarnings(library(grid)))       # plot anotation
suppressMessages(suppressWarnings(library(readxl)))     # data importation
suppressMessages(suppressWarnings(library(ggcorrplot))) # ggplot add on
analysis_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(analysis_path))
analysis_dir <- getwd()

#===============#
# Barplots, relative model fits
#===============#

#
Group <- c("HC", "HC", "HC",
           "ON", "ON", "ON",
           "OFF","OFF","OFF")

Models <- c("Mod AA", "Mod AB", "Mod AC",                        # Model A = HDDM  (group A = HC...)
            "Mod BA", "Mod BB", "Mod BC",                        # Model B = HDDDM
            "Mod CA", "Mod CB", "Mod CC")                        # Model C = EWHDDDM

# (previously acquired in SocDBS_Comparison_HDDM.r)
ELPD_diff <- c(-308.765, -339.066, 0,
               -182.697, -213.155, 0,
               -142.551, -163.934, 0)

SE_diff <- c(68.609, 61.556, 0,
             54.724, 60.464, 0,
             59.219, 62.893, 0)

LowBounds <- ELPD_diff - SE_diff
HighBounds <- ELPD_diff + SE_diff
#
DF_fits <- data.frame(Group, Models, ELPD_diff, SE_diff, LowBounds, HighBounds)
DF_fits
# annotations
sig_B_OF <- grobTree(textGrob("σ = 2.61", x=0.55,  y=0.88, hjust=0, gp=gpar(col="black", fontsize=10)))   
sig_A_OF <- grobTree(textGrob("σ = 2.41", x=0.59,  y=0.77, hjust=0, gp=gpar(col="black", fontsize=10)))
sig_B_ON <- grobTree(textGrob("σ = 3.52", x=0.43,  y=0.56, hjust=0, gp=gpar(col="black", fontsize=10)))
sig_A_ON <- grobTree(textGrob("σ = 3.33", x=0.49,  y=0.45, hjust=0, gp=gpar(col="black", fontsize=10)))
sig_B_HC <- grobTree(textGrob("σ = 5.51", x=0.14,  y=0.235, hjust=0, gp=gpar(col="black", fontsize=10)))
sig_A_HC <- grobTree(textGrob("σ = 4.50", x=0.21,  y=0.12, hjust=0, gp=gpar(col="black", fontsize=10)))
#
pl0 <- ggplot(data = DF_fits, aes(x=ELPD_diff, y = Models, fill= Group, color = Group)) +
              geom_point(size=5) +
              geom_errorbar(aes(xmin  = LowBounds,xmax  = HighBounds), width = 0.0, size  = 2) +
              scale_fill_manual(values = c("grey", "purple1", "orange1")) +
              scale_color_manual(values = c("grey", "purple1", "orange1")) +
              geom_vline(aes(xintercept=0), size=1) +
              geom_hline(aes(yintercept=3.5), size=0.5, linetype="dashed") +
              geom_hline(aes(yintercept=6.5), size=0.5, linetype="dashed") +
              labs(y="Models", x="ELPD Difference") +
              annotation_custom(sig_B_OF) +
              annotation_custom(sig_A_OF) +
              annotation_custom(sig_B_ON) +
              annotation_custom(sig_A_ON) +
              annotation_custom(sig_B_HC) +
              annotation_custom(sig_A_HC) +
              theme_pubr()  + theme(legend.position = c(0.1, 0.85))
pl0

#===============#
# Correlograms of posterior param values
#===============#

params <- c("k", "f", "a", "ter", "vcoeff")

# HC
data_healthy <- file.path(analysis_dir, "posterior_draws_healthy.xlsx")
data_HC <- read_xlsx(data_healthy)
DF <- as.data.frame(data_HC)

hyperparameters_HC <- data.frame(DF$hyper_k, DF$hyper_f, DF$hyper_a, DF$hyper_ter, DF$hyper_v)
colnames(hyperparameters_HC) <- params

cormat <- round(cor(hyperparameters_HC), 3)
head(cormat[, 1:5])
ggcorrplot(cormat, method = "circle", colors = c("purple1", "white", "orange1"))

# ON
data_ON <- file.path(analysis_dir, "posterior_draws_DBS_ON.xlsx")
data_ON <- read_xlsx(data_ON)
DF <- as.data.frame(data_ON)

hyperparameters_ON <- data.frame(DF$hyper_k, DF$hyper_f, DF$hyper_a, DF$hyper_ter, DF$hyper_v)
colnames(hyperparameters_ON) <- params

cormat <- round(cor(hyperparameters_ON), 3)
head(cormat[, 1:5])
ggcorrplot(cormat, method = "circle", colors = c("purple1", "white", "orange1"))

# OFF
data_OFF <- file.path(analysis_dir, "posterior_draws_DBS_OFF.xlsx")
data_OFF <- read_xlsx(data_OFF)
DF <- as.data.frame(data_OFF)

hyperparameters_OFF <- data.frame(DF$hyper_k, DF$hyper_f, DF$hyper_a, DF$hyper_ter, DF$hyper_v)
colnames(hyperparameters_OFF) <- params

cormat <- round(cor(hyperparameters_OFF), 3)
head(cormat[, 1:5])
ggcorrplot(cormat, method = "circle", colors = c("purple1", "white", "orange1"))

#==============================================================================#
#------------------------------------------------------------------------------#
#==============================================================================#