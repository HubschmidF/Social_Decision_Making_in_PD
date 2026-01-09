#==============================================================================# 
# SocialDBS: Mixed effect models group comparisons
#==============================================================================#

# FH

# Description:
#  - Fitting and results of the Mixed effects models
#  - Generates plots for the LMEs

#==============================================================================#

#===============#
# Path, Packages and data curation
#===============#

# Path and software packages
suppressMessages(suppressWarnings(library(rstudioapi))) # utility for path managment
suppressMessages(suppressWarnings(library(tidyverse)))  # utility for data manipulation
suppressMessages(suppressWarnings(library(lme4)))       # used for fitting of LMMs
suppressMessages(suppressWarnings(library(effects)))    # used for quick visualisation only
suppressMessages(suppressWarnings(library(lmerTest)))   # used for significance of fixed effects
suppressMessages(suppressWarnings(library(emmeans)))    # used for contrast post-hoc comparisons
suppressMessages(suppressWarnings(library(ggplot2)))    # used for plotting
suppressMessages(suppressWarnings(library(ggpubr)))     # ggplot add on
suppressMessages(suppressWarnings(library(grid)))       # ggplot add on
analysis_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(analysis_path))
analysis_dir <- getwd()

# Data file
data_file <- file.path(analysis_dir, "Task_Data_Full.csv")
data <- read.csv(data_file, header = T, sep = ",")
DF <- as.data.frame(data)

# To account for inherently skewed RTs
DF$RT <- DF$rt/1000             # to have rt in seconds
DF$log_rt <- log(DF$RT)         # Log and the base 10 transformation both make the RT normal so log transformation kept for further analyses

# to make LMMs more interpretables
DF$stim_choice <- as.factor(DF$response) 
print(levels(DF$stim_choice))

# it seems like one response was logged as "-1" wich should be impossible 
# given the PTB script. For convinience we just do not account for that one trial
# ordrered as  -1 , 0, 1, 2 

levels(DF$stim_choice) <- c("miss", "missed", "selfish", "generous")

# Exclude missed responses and the one weird trial that was logged
all_trials <- subset(DF, stim_choice %in% c("selfish","generous"))  
excluded_trials <- subset(DF, stim_choice %in% c("whatever","missed")) 
print(nrow(excluded_trials)) 

# Reward difference "rd"
# indication of the difference the social and selfish amount, i.e. it increases relative inequity 
# if one person gets 120 and the other nothing.
all_trials$RD <- all_trials$selfish_amount - 75

# Re-organises the mod_choice so that 1 corresponds to generocity and 0 to selfish
all_trials$bin_choice <- as.factor(all_trials$response)
levels(all_trials$bin_choice)  <- c("selfish", "generous")   #
all_trials$bin_choice <- relevel(all_trials$bin_choice, ref = "selfish") #
levels(all_trials$bin_choice) <- c("0", "1")
all_trials$mod_choice <- as.integer(all_trials$bin_choice)
all_trials$mod_choice <- as.integer(all_trials$mod_choice) -1

#===============#
# Linear Mixed Model: log(RT)
#===============#

# Fit and visualisation
LMM_rt <- lmer(log_rt ~ group + social_distance + RD + (1 | Part_idx),
               data = all_trials)

confint(LMM_rt)
summary(LMM_rt)
plot(allEffects(LMM_rt))       # visualisation only

# Contrast comparisons between groups
(LMM_group <- pairs(emmeans(LMM_rt, ~ group,type="log_rt"))%>%test())

#===============#
# Generalized Mixed Logistic Model: P(generous)
#===============#

# Fit and visualisation
GLM_choice <- glmer(mod_choice ~ group + social_distance + RD + (1 | Part_idx),
                    data = all_trials,
                    family = binomial(link='logit'))
confint(GLM_choice)
summary(GLM_choice)
plot(allEffects(GLM_choice))  # visualisation only

# Contrast comparisons between groups
(GLM_group <- pairs(emmeans(GLM_choice, ~ group,type="mod_choice"))%>%test())

#===============#
# Plots (model expected values)
#===============#

# Significance annotations
sig1 <- grobTree(textGrob("*", x=0.7,  y=0.8, hjust=0, gp=gpar(col="black", fontsize=30, fontface="bold")))
sig2 <- grobTree(textGrob("**", x=0.7,  y=0.8, hjust=0, gp=gpar(col="black", fontsize=30, fontface="bold")))
sig3 <- grobTree(textGrob("***", x=0.7,  y=0.8, hjust=0, gp=gpar(col="black", fontsize=30, fontface="bold")))
ns <- grobTree(textGrob("n.s", x=0.7,  y=0.8, hjust=0, gp=gpar(col="black", fontsize=20, fontface="bold"))) 

# Fixed effect of social distance on P(generous)
effect_sd <- effects::effect(term= "social_distance", mod= GLM_choice)
summary(effect_sd) 
df_sd<- as.data.frame(effect_sd)

pl1 <- ggplot() +
              #geom_point(data=df_sd, aes(x=social_distance, y=fit), color="black") +   # to add points
              geom_line(data=df_sd, aes(x=social_distance, y=fit), size = 2, color="navy") +
              geom_ribbon(data= df_sd, aes(x=social_distance, ymin=lower, ymax=upper), alpha= 0.3, fill="lightblue") +
              scale_color_manual(values = c("navy")) +
              labs(y = 'P(Generous)', x = 'Social Distance') +
              theme_pubr() +
              annotation_custom(sig3)
pl1

# Fixed effect of the reward difference on P(generous)
effect_rd <- effects::effect(term= "RD", mod= GLM_choice)
summary(effect_rd) 
df_rd<- as.data.frame(effect_rd)

pl2 <- ggplot() +
              #geom_point(data=df_sd, aes(x=social_distance, y=fit), color="black") +   # to add points
              geom_line(data=df_rd, aes(x=RD, y=fit), size = 2, color="navy") +
              geom_ribbon(data= df_rd, aes(x=RD, ymin=lower, ymax=upper), alpha= 0.3, fill="lightblue") +
              scale_color_manual(values = c("navy")) +
              labs(y = 'P(Generous)', x = 'Δ Reward') +
              theme_pubr() +
              annotation_custom(sig3)
pl2

# Fixed effect of social distance on log(rt)
eff_sd <- effects::effect(term= "social_distance", mod= LMM_rt)
summary(eff_sd) 
d_sd<- as.data.frame(eff_sd)

pl3 <- ggplot() +
              #geom_point(data=df_sd, aes(x=social_distance, y=fit), color="black") +   # to add points
              geom_line(data=d_sd, aes(x=social_distance, y=fit), size = 2, color="navy") +
              geom_ribbon(data= d_sd, aes(x=social_distance, ymin=lower, ymax=upper), alpha= 0.3, fill="lightblue") +
              scale_color_manual(values = c("navy")) +
              labs(y = 'Log(RT)', x = 'Social Distance') +
              theme_pubr() +
              annotation_custom(sig2)
pl3

# Fixed effect of reward difference on log(rt)
eff_rd <- effects::effect(term= "RD", mod= LMM_rt)
summary(eff_rd) 
d_rd<- as.data.frame(eff_rd)

pl4 <- ggplot() +
              #geom_point(data=df_sd, aes(x=social_distance, y=fit), color="black") +   # to add points
              geom_line(data=d_rd, aes(x=RD, y=fit), size = 2, color="navy") +
              geom_ribbon(data= d_rd, aes(x=RD, ymin=lower, ymax=upper), alpha= 0.3, fill="lightblue") +
              scale_color_manual(values = c("lightblue")) +
              labs(y = 'Log(RT)', x = 'Δ Reward') +
              theme_pubr() +
              annotation_custom(ns)
pl4

#===============#
# Group effect plots (model expected values)
#===============#

# Predicted generocity proportions 
effect_gr <- effects::effect(term= "group", mod= GLM_choice)
print(summary(effect_gr)) # to get values

group_plot_df <- data.frame(group = c("HC", "OFF", "ON"),
                            value = c(0.7551868,0.5666841,0.5593935), # said values
                            upper = c(0.8349619,0.6752946,0.6687772),
                            lower = c(0.6528815,0.4512649,0.4439226))

# annotations
bar1 <- grobTree(textGrob("_______",            x=0.20,  y=0.90, hjust=0, gp=gpar(col="black", fontsize=20, fontface="italic")))
bar2 <- grobTree(textGrob("________________",   x=0.20,  y=0.97, hjust=0, gp=gpar(col="black", fontsize=20, fontface="italic")))
bar3 <- grobTree(textGrob("________",           x=0.50,  y=0.63, hjust=0, gp=gpar(col="black", fontsize=20, fontface="italic")))
sig1 <- grobTree(textGrob("*",                  x=0.45,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=30, fontface="bold")))
sig2 <- grobTree(textGrob("*",                  x=0.31,  y=0.87, hjust=0, gp=gpar(col="black", fontsize=30, fontface="bold")))
sig3 <- grobTree(textGrob("n.s",                x=0.60,  y=0.63, hjust=0, gp=gpar(col="black", fontsize=18, fontface="bold")))

pl5 <- ggplot(group_plot_df, aes(x=group, y=value, color=group)) + 
              geom_errorbar(aes(ymin=lower, ymax=upper), size = 2, width=.1) +
              geom_point(size = 6) +
              scale_color_manual(values = c("darkgrey", "purple1", "orange1")) +
              labs(y = 'P(Generous)', x = 'Group') +
              ylim(0.4,1.0) +
              theme_pubr() +
              annotation_custom(bar1) +     # significance from contrast comparisons
              annotation_custom(bar2) +
              annotation_custom(bar3) +
              annotation_custom(sig1) +
              annotation_custom(sig2) +
              annotation_custom(sig3) +
              theme(legend.position='none')
pl5

# Predicted log(RTs)
effect_gr <- effects::effect(term= "group", mod= LMM_rt)
print(summary(effect_gr)) # to get values

group_plot_df <- data.frame(group = c("HC", "OFF", "ON"),
                            value = c(0.8543679,0.7982794,0.6882923), # said values
                            upper = c(1.0518383,0.9865709,0.8765945),
                            lower = c(0.6568976,0.6099878,0.4999902))

# annotations
bar1 <- grobTree(textGrob("________",           x=0.20,  y=0.90, hjust=0, gp=gpar(col="black", fontsize=20, fontface="italic")))
bar2 <- grobTree(textGrob("________________",   x=0.20,  y=0.97, hjust=0, gp=gpar(col="black", fontsize=20, fontface="italic")))
bar3 <- grobTree(textGrob("________",           x=0.51,  y=0.80, hjust=0, gp=gpar(col="black", fontsize=20, fontface="italic")))
sig1 <- grobTree(textGrob("n.s",                x=0.42,  y=0.98, hjust=0, gp=gpar(col="black", fontsize=20, fontface="bold")))
sig2 <- grobTree(textGrob("n.s",                x=0.30,  y=0.90, hjust=0, gp=gpar(col="black", fontsize=20, fontface="bold")))
sig3 <- grobTree(textGrob("***",                x=0.60,  y=0.77, hjust=0, gp=gpar(col="black", fontsize=30, fontface="bold")))

pl6 <- ggplot(group_plot_df, aes(x=group, y=value, color=group)) + 
              geom_errorbar(aes(ymin=lower, ymax=upper), size = 2, width=.1) +
              geom_point(size = 6) +
              scale_color_manual(values = c("darkgrey", "purple1", "orange1")) +
              labs(y = 'Log(RT)', x = 'Group') +
              ylim(0.4,1.2) +
              theme_pubr() +
              annotation_custom(bar1) +     # significance from contrast comparisons
              annotation_custom(bar2) +
              annotation_custom(bar3) +
              annotation_custom(sig1) +
              annotation_custom(sig2) +
              annotation_custom(sig3) +
              theme(legend.position='none') 
pl6

#===============#
# Manipulation check (splitting group variable)  # Effect stays significant
#===============#

# Manipulation check (if split the within group from between...)
HC_ON <- subset(all_trials, group %in% c("ON","NO"))  

GLM_choice <- glmer(mod_choice ~ group + social_distance + RD + (1 | Part_idx),
                    data = HC_ON,
                    family = binomial(link='logit'))
summary(GLM_choice)
plot(allEffects(GLM_choice))


HC_OFF <- subset(all_trials, group %in% c("OFF","NO"))  


GLM_choice <- glmer(mod_choice ~ group + social_distance + RD + (1 | Part_idx),
                    data = HC_OFF,
                    family = binomial(link='logit'))
summary(GLM_choice)
plot(allEffects(GLM_choice))

# ============================================================================ #
# ---------------------------------------------------------------------------- #
# ============================================================================ #