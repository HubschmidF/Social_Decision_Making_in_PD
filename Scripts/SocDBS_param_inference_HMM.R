#==============================================================================# 
# SocialDBS: Posterior Parameters and parameter inference (HMM)
#==============================================================================#

# FH

#==============================================================================#

#===============#
# Path, packages and data
#===============#

# Software packages
suppressMessages(suppressWarnings(library(rstudioapi)))     # path management
suppressMessages(suppressWarnings(library(HDInterval)))     # functions for highest density intervals
suppressMessages(suppressWarnings(library(ggplot2)))        # plots
suppressMessages(suppressWarnings(library(ggpubr)))         # figure generation (probably not final used figure)
suppressMessages(suppressWarnings(library(readxl)))         # to read the data
suppressMessages(suppressWarnings(library(ggridges)))       # annotating plots 
suppressMessages(suppressWarnings(library(tibble)))         # data manipulation
suppressMessages(suppressWarnings(library(grid)))           # plot annotation
suppressMessages(suppressWarnings(library(plyr)))           # data manipulation
suppressMessages(suppressWarnings(library(BayesFactor)))    # BF estimation

# Path
analysis_path = getActiveDocumentContext()$path
setwd(dirname(analysis_path))
analysis_dir <- getwd()

# Data files --
# HC
data_file <- file.path(analysis_dir, "000_HMM_results.csv")
data <- read.csv(data_file, header = T, sep = ",")
DF <- as.data.frame(data)

# Log transformation of parameters
DF$log_k     <- log(DF$k)
DF$log_lapse <- log(DF$lapse)
DF$log_recov <- log(DF$recov)

#===============#
# Posterior parameters
#===============#

# Discount rate
log_k <- ggplot(DF, aes(x = group, y = log_k, fill = group)) +
                geom_point(position = position_jitter(width = 0.05), size = 1.5, alpha = 0.5) +
                geom_violin(width = 0.6, alpha = 0.4) +
                stat_summary(fun = median,
                             fun.min = function(x) quantile(x, 0.25),
                             fun.max = function(x) quantile(x, 0.75),
                              geom = "crossbar", width = 0.05, fill = "white", alpha = 0.8) +
                scale_fill_manual(values = c("grey", "purple1", "orange1")) +
                theme_pubr() +
                labs(x = "Group", y = "Discount log(k)") +
                theme(legend.position = "none")
log_k


# Lapse probability
p_lapse <- ggplot(DF, aes(x = group, y = log_lapse, fill = group)) +
                  geom_point(position = position_jitter(width = 0.05), size = 1.5, alpha = 0.5) +
                  geom_violin(width = 0.6, alpha = 0.4) +
                  stat_summary(fun = median,
                               fun.min = function(x) quantile(x, 0.25),
                               fun.max = function(x) quantile(x, 0.75),
                               geom = "crossbar", width = 0.05, fill = "white", alpha = 0.8) +
                  scale_fill_manual(values = c("grey", "purple1", "orange1")) +
                  theme_pubr() +
                  labs(x = "Group", y = "log(T(lapse))") +
                  theme(legend.position = "none")
p_lapse

# Recovery probability
p_recov <- ggplot(DF, aes(x = group, y = log_recov, fill = group)) +
                  geom_point(position = position_jitter(width = 0.05), size = 1.5, alpha = 0.5) +
                  geom_violin(width = 0.6, alpha = 0.4) +
                  stat_summary(fun = median,
                               fun.min = function(x) quantile(x, 0.25),
                               fun.max = function(x) quantile(x, 0.75),
                               geom = "crossbar", width = 0.05, fill = "white", alpha = 0.8) +
                  scale_fill_manual(values = c("grey", "purple1", "orange1")) +
                  theme_pubr() +
                  labs(x = "Group", y = "log(T(recover))") +
                  theme(legend.position = "none")
p_recov

# Figure with additional plot
Util_f <- ggarrange(log_k,  p_lapse,  p_recov, pl5,   # pl 5 comes from another script
                    labels = c('','','',''),
                    legend = "none", common.legend = TRUE, nrow = 1, ncol = 4)
Util_f

# ============================================================================ #
# ---------------------------------------------------------------------------- #
# ============================================================================ # 