#==============================================================================# 
# SocialDBS: Posterior Parameters and parameter inference
#==============================================================================#

# FH

# Description:

#  - plots posterior draws of hyper parameters of EWHDDDM
#  - inference on posterior draws based on difference densities (95% HDI)

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
suppressMessages(suppressWarnings(library(plyr)))           # summarizing functions

# Path
analysis_path = getActiveDocumentContext()$path
setwd(dirname(analysis_path))
analysis_dir <- getwd()

# Data files --
# HC
data_healthy <- file.path(analysis_dir, "posterior_draws_healthy.xlsx")
data_HC <- read_xlsx(data_healthy)
DF_HC <- as.data.frame(data_HC)
# DBS ON 
data_PATON <- file.path(analysis_dir, "posterior_draws_DBS_ON.xlsx")
data_ON <- read_xlsx(data_PATON)
DF_ON <- as.data.frame(data_ON)
# DBS OFF 
data_PATOFF <- file.path(analysis_dir, "posterior_draws_DBS_OFF.xlsx")
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

# non-decision times (ter)
terHC <- DF_HC$hyper_ter
terON <- DF_ON$hyper_ter
terOFF <- DF_OFF$hyper_ter
#
ters <- c(terHC, terON, terOFF)

TER <- data.frame(Group, ters)
mu1 <- ddply(TER, "Group", summarise, grp.mean=mean(ters))

pl1 <- ggplot(TER, aes(x = ters, color = Group, fill = Group)) +
              geom_density(alpha = 0.2) +
              geom_vline(data=mu1, aes(xintercept=grp.mean, color=Group), linetype = 'dashed') +
              labs(y = 'Density', x = 'Non-Decision Time (ter)') +
              scale_fill_manual(values=c("grey", "orange1", "purple1")) +
              scale_color_manual(values=c("grey", "orange1", "purple1")) +
              theme_pubr() +
              theme(legend.position='right')
pl1

# boundaries (a)
boundHC <- DF_HC$hyper_a
boundON <- DF_ON$hyper_a
boundOFF <- DF_OFF$hyper_a
#
bounds <- c(boundHC, boundON, boundOFF)

BOUND <- data.frame(Group, bounds)
mu2 <- ddply(BOUND, "Group", summarise, grp.mean=mean(bounds))

pl2 <- ggplot(BOUND, aes(x = bounds, color = Group, fill = Group)) +
              geom_density(alpha = 0.2) +
              geom_vline(data=mu2, aes(xintercept=grp.mean, color=Group), linetype = 'dashed') +
              labs(y = 'Density', x = 'Boundary Separation (a)') +
              scale_fill_manual(values=c("grey", "orange1", "purple1")) +
              scale_color_manual(values=c("grey", "orange1", "purple1")) +
              theme_pubr() +
              theme(legend.position='right') 
pl2

# Drift rate scalling (v)
driftHC <- DF_HC$hyper_v
driftON <- DF_ON$hyper_v
driftOFF <- DF_OFF$hyper_v
#
drifts <- c(driftHC, driftON, driftOFF)

DRIFT <- data.frame(Group, drifts)
mu3 <- ddply(DRIFT, "Group", summarise, grp.mean=mean(drifts))

pl3 <- ggplot(DRIFT, aes(x = drifts, color = Group, fill = Group)) +
              geom_density(alpha = 0.2) +
              geom_vline(data=mu3, aes(xintercept=grp.mean, color=Group), linetype = 'dashed') +
              labs(y = 'Density', x = 'Drift Rate Scalling (vcoeff)') +
              scale_fill_manual(values=c("grey", "orange1", "purple1")) +
              scale_color_manual(values=c("grey", "orange1", "purple1")) +
              theme_pubr() +
              theme(legend.position='right')
pl3

# Discount
k_HC <- DF_HC$log_k
k_ON <- DF_ON$log_k
k_OFF <- DF_OFF$log_k
#
discounts <- c(k_HC, k_ON, k_OFF)

DISCOUNT <- data.frame(Group, discounts)
mu4 <- ddply(DISCOUNT, "Group", summarise, grp.mean=mean(discounts))

pl4 <- ggplot(DISCOUNT, aes(x = discounts, color = Group, fill = Group)) +
              geom_density(alpha = 0.2) +
              geom_vline(data=mu4, aes(xintercept=grp.mean, color=Group), linetype = 'dashed') +
              labs(y = 'Density', x = 'Discount (log(k))') +
              scale_fill_manual(values=c("grey", "orange1", "purple1")) +
              scale_color_manual(values=c("grey", "orange1", "purple1")) +
              theme_pubr() +
              theme(legend.position='right')
pl4

# Discount
f_HC <- DF_HC$hyper_f
f_ON <- DF_ON$hyper_f
f_OFF <- DF_OFF$hyper_f
#
weights <- c(f_HC, f_ON, f_OFF)

WEIGHTS <- data.frame(Group, weights)
mu5 <- ddply(WEIGHTS, "Group", summarise, grp.mean=mean(weights))

pl5 <- ggplot(WEIGHTS, aes(x = weights, color = Group, fill = Group)) +
              geom_density(alpha = 0.2) +
              geom_vline(data=mu5, aes(xintercept=grp.mean, color=Group), linetype = 'dashed') +
              labs(y = 'Density', x = 'Empathy Weight (f)') +
              scale_fill_manual(values=c("grey", "orange1", "purple1")) +
              scale_color_manual(values=c("grey", "orange1", "purple1")) +
              theme_pubr() +
              theme(legend.position='right')
pl5

#===============#
# Difference densities 
#===============#

# Non-decision times (difference distribution) ---

ter_HC_ON <- DF_HC$hyper_ter - DF_ON$hyper_ter
ter_HC_OFF <- DF_HC$hyper_ter - DF_OFF$hyper_ter 
ter_OFF_ON <-  DF_OFF$hyper_ter - DF_ON$hyper_ter
# Utility
DF_ter_HC_ON <- tibble(ter_HC_ON)
DF_ter_HC_OFF <- tibble(ter_HC_OFF)
DF_ter_OFF_ON <- tibble(ter_OFF_ON)

# 95% HDI
# comparison 1
mean(ter_HC_ON)                  # 0.3597574
hdi(ter_HC_ON, credMass = .95)   # 0.2782701 0.4375111     ***
# comparison 2
mean(ter_HC_OFF)                 # 0.2643057
hdi(ter_HC_OFF, credMass = .95)  # 0.1745258 0.3565558     ***
# comparison 3
mean(ter_OFF_ON)                 # 0.09545173
hdi(ter_OFF_ON, credMass = .95)  # 0.003538324 0.198896887 ***

# corresponding plots ---

# comp 1
pl1A <- ggplot(DF_ter_HC_ON, aes(x = ter_HC_ON, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - ON") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr() + annotation_custom(sig)
pl1A

# comp 2
pl1B <- ggplot(DF_ter_HC_OFF, aes(x = ter_HC_OFF, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - OFF") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr() + annotation_custom(sig)
pl1B

# comp 3
pl1C <- ggplot(DF_ter_OFF_ON , aes(x = ter_OFF_ON, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "OFF - ON") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr() + annotation_custom(sig)
pl1C

# boundaries (difference distribution) ---

a_HC_ON <- DF_HC$hyper_a - DF_ON$hyper_a
a_HC_OFF <- DF_HC$hyper_a - DF_OFF$hyper_a
a_ON_OFF <-  DF_ON$hyper_a - DF_OFF$hyper_a
# Utility
DF_a_HC_ON <- tibble(a_HC_ON)
DF_a_HC_OFF <- tibble(a_HC_OFF)
DF_a_ON_OFF <- tibble(a_ON_OFF)

# comparison 1
mean(a_HC_ON)                  # 0.355299
hdi(a_HC_ON, credMass = .95)   # 0.001294053 0.705616284     ***
# comparison 2
mean(a_HC_OFF)                 # 0.3144333
hdi(a_HC_OFF, credMass = .95)  # -0.09035118  0.72915745
# comparison 3
mean(a_ON_OFF)                 # -0.04086573
hdi(a_ON_OFF, credMass = .95)  # -0.4671197  0.3474664 

# corresponding plots ---

# comp 1
pl2A <- ggplot(DF_a_HC_ON, aes(x = a_HC_ON, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - ON") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr() + annotation_custom(sig)
pl2A

# comp 2
pl2B <- ggplot(DF_a_HC_OFF, aes(x = a_HC_OFF, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - OFF") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr() 
pl2B

# comp 3
pl2C <- ggplot(DF_a_ON_OFF, aes(x = a_ON_OFF, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "ON - OFF") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr()
pl2C

# Drift rates (difference distribution) ---

v_HC_ON <- DF_HC$hyper_v - DF_ON$hyper_v
v_HC_OFF <- DF_HC$hyper_v - DF_OFF$hyper_v
v_ON_OFF <-  DF_ON$hyper_v - DF_OFF$hyper_v
# Utility
DF_v_HC_ON <- tibble(v_HC_ON)
DF_v_HC_OFF <- tibble(v_HC_OFF)
DF_v_ON_OFF <- tibble(v_ON_OFF)

# comparison 1
mean(v_HC_ON)                  # 0.003271509
hdi(v_HC_ON, credMass = .95)   # 7.767546e-05 6.550787e-03    ***
# comparison 2
mean(v_HC_OFF)                 # 0.005028542
hdi(v_HC_OFF, credMass = .95)  # 0.001920854 0.008007923      ***
# comparison 3
mean(v_ON_OFF)                 # 0.001757033
hdi(v_ON_OFF, credMass = .95)  # -0.0005505653  0.0040603868 

# corresponding plots ---

# comp 1
pl3A <- ggplot(DF_v_HC_ON, aes(x = v_HC_ON, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - ON") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr() + annotation_custom(sig)
pl3A

# comp 2
pl3B <- ggplot(DF_v_HC_OFF, aes(x = v_HC_OFF, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - OFF") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr() + annotation_custom(sig)
pl3B

# comp 3
pl3C <- ggplot(DF_v_ON_OFF, aes(x = v_ON_OFF, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "ON - OFF") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr()
pl3C

# Log(k) (difference distribution) ---

k_HC_ON <- DF_HC$log_k - DF_ON$log_k
k_HC_OFF <- DF_HC$log_k - DF_OFF$log_k
k_ON_OFF <-  DF_ON$log_k - DF_OFF$log_k
# Utility
DF_k_HC_ON <- tibble(k_HC_ON)
DF_k_HC_OFF <- tibble(k_HC_OFF)
DF_k_ON_OFF <- tibble(k_ON_OFF)

# comparison 1
mean(k_HC_ON)                  # -0.6379484
hdi(k_HC_ON, credMass = .95)   # -1.3730487  0.1018019
# comparison 2
mean(k_HC_OFF)                 # -0.8232134
hdi(k_HC_OFF, credMass = .95)  # -1.53898398 -0.06225657    ***
# comparison 3
mean(k_ON_OFF)                 # -0.1852649
hdi(k_ON_OFF, credMass = .95)  # -0.9761111  0.6290493

# corresponding plots ---

# comp 1
pl4A <- ggplot(DF_k_HC_ON, aes(x = k_HC_ON, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - ON") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr()
pl4A

# comp 2
pl4B <- ggplot(DF_k_HC_OFF, aes(x = k_HC_OFF, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - OFF") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr() + annotation_custom(sig)
pl4B

# comp 3
pl4C <- ggplot(DF_k_ON_OFF, aes(x = k_ON_OFF, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "ON - OFF") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr()
pl4C

# EMpathy weight (difference distribution) ---

f_HC_ON <- DF_HC$hyper_f - DF_ON$hyper_f
f_HC_OFF <- DF_HC$hyper_f - DF_OFF$hyper_f
f_ON_OFF <-  DF_ON$hyper_f - DF_OFF$hyper_f
# Utility
DF_f_HC_ON <- tibble(f_HC_ON)
DF_f_HC_OFF <- tibble(f_HC_OFF)
DF_f_ON_OFF <- tibble(f_ON_OFF)

# comparison 1
mean(f_HC_ON)                  # -0.0949109
hdi(f_HC_ON, credMass = .95)   # -0.9606738  0.7552719
# comparison 2
mean(f_HC_OFF)                 # -0.7551205
hdi(f_HC_OFF, credMass = .95)  # -1.6773834  0.1298565
# comparison 3
mean(f_ON_OFF)                 # -0.6602096
hdi(f_ON_OFF, credMass = .95)  # -1.6441830  0.4321178 

# corresponding plots ---

# comp 1
pl5A <- ggplot(DF_f_HC_ON, aes(x = f_HC_ON, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - ON") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr()
pl5A

# comp 2
pl5B <- ggplot(DF_f_HC_OFF, aes(x =f_HC_OFF, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "HC - OFF") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr() 
pl5B

# comp 3
pl5C <- ggplot(DF_f_ON_OFF, aes(x = f_ON_OFF, y = 0, fill = stat(quantile))) +
               geom_density_ridges_gradient(quantile_lines = TRUE, quantile_fun = hdi, vline_linetype = 2) +
               geom_vline(xintercept = 0, color = 'red', size = 1) +
               labs(y = "Diff. Density", x = "ON - OFF") +
               scale_fill_manual(values = alpha(c("transparent", "navy", "transparent"), 0.5), guide = "none") +
               theme_pubr()
pl5C

#===============#
# Difference densities
#===============#

Figure <- ggarrange(pl1,  pl2,  pl3,  pl4,  pl5,
                    pl1A, pl2A, pl3A, pl4A, pl5A,
                    pl1B, pl2B, pl3B, pl4B, pl5B,
                    pl1C, pl2C, pl3C, pl4C, pl5C,
                    labels = c('','','','','',
                               '','','','','',
                               '','','','','',
                               '','','','',''),
                    legend = "top", common.legend = TRUE,
                    nrow = 4, ncol = 5)
Figure

#==============================================================================#
#------------------------------------------------------------------------------#
#==============================================================================#