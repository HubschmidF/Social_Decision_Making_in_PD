#==============================================================================# 
# SocialDBS: Model extracted quantities: Q values and conflict
#==============================================================================#

# FH

                           # !!! NOTE !!! #

# Note: To plot model estimated Qvalues, you need to extract modelled values
# from stan summary output file and add them to the original dataset. Finally,
# adapt filename line 41.

# Description:
#  - plots modeled Q values and modeled conflict
#  - Situations of increased conflict (GLM / LMM) (not kept in manuscript)

#==============================================================================#

#===============#
# Path, packages and data
#===============#

suppressMessages(suppressWarnings(library(rstudioapi)))     # path management
suppressMessages(suppressWarnings(library(tictoc)))         # utility to estimate sampling times
suppressMessages(suppressWarnings(library(dplyr)))          # utility for data frame manipulations
suppressMessages(suppressWarnings(library(plyr)))           # utility for data frame manipulations
suppressMessages(suppressWarnings(library(ggplot2)))        # used for plotting
suppressMessages(suppressWarnings(library(ggpubr)))         # ggplot add on
suppressMessages(suppressWarnings(library(HDInterval)))     # used for HDI estimation
suppressMessages(suppressWarnings(library(lme4)))           # used for fitting of LMMs
suppressMessages(suppressWarnings(library(effects)))        # used for quick visualisation only
suppressMessages(suppressWarnings(library(lmerTest)))       # used for significance of fixed effects
suppressMessages(suppressWarnings(library(emmeans)))        # used for contrast post-hoc comparisons

# Path
analysis_path = getActiveDocumentContext()$path
setwd(dirname(analysis_path))
analysis_dir <- getwd()

# Data file
data_file <- file.path(analysis_dir, "FileName.csv")    # ADAPT 
data <- read.csv(data_file, header = T, sep = ",")
DF <- as.data.frame(data)

DF$RT <- DF$rt/1000
DF$stim_choice <- as.factor(DF$response)
print(levels(DF$stim_choice))
levels(DF$stim_choice) <- c("selfish", "generous")

# for utility
all_trials <- subset(DF, stim_choice %in% c("selfish","generous"))
all_trials$log_q_generous <- log(all_trials$Q_generous)
all_trials$log_selfish <-  log(all_trials$selfish_amount)

all_trials$conflict <- all_trials$Q_generous - all_trials$selfish_amount
all_trials$conflict_load <- all_trials$log_q_generous - all_trials$log_selfish
all_trials$absolute_load <- abs(all_trials$conflict_load)
all_trials$RD <- all_trials$selfish_amount - 75

# for legends
all_trials$Group <- as.factor(all_trials$group)
levels(all_trials$Group) <- c("HC", "OFF", "ON")

# for LMM analyses
# Re-organises the mod_choice so that 1 corresponds to generocity and 0 to selfish
all_trials$bin_choice <- as.factor(all_trials$response)
levels(all_trials$bin_choice)  <- c("selfish", "generous")   #
all_trials$bin_choice <- relevel(all_trials$bin_choice, ref = "selfish") #
levels(all_trials$bin_choice) <- c("0", "1")
all_trials$mod_choice <- as.integer(all_trials$bin_choice)
all_trials$mod_choice <- as.integer(all_trials$mod_choice) - 1

one_ammount <- subset(all_trials, selfish_amount %in% c(85))
one_ammount <- subset(one_ammount, button_choice %in% c(1))

#===============#
# Plots Q values and value differences
#===============#

# Q generous as a function of the social distance
pl1 <- ggplot(one_ammount, aes(x = social_distance, y = Q_generous, color = Group)) +
             geom_line(aes(x = social_distance, y =  Q_generous, group = log_idx, color = Part_idx), color = c("lightgrey"), linewidth = 0.5) +
             geom_point(alpha = 0.9, size = 3) +
             labs(y = 'Q generous', x = 'Social Distance') +
             scale_color_manual(values = c("darkgrey", "purple1", "orange1")) +
             #ylim(0,500) +
             theme_pubr() 
pl1

# Distribution of decision conflict (not kept in final version)
pl2 <- ggplot(all_trials, aes(x = conflict_load, color = Group, fill = Group)) + 
              geom_density(alpha = 0.2) +
              geom_vline(xintercept = 0, color = 'black', size = 1, linetype = "dashed") +
              scale_color_manual(values = c("darkgrey", "purple1", "orange1")) +
              scale_fill_manual(values = c("darkgrey", "purple1", "orange1")) +
              labs(y = 'Density', x = 'Decisional Conflict (a.u)') +
              ylim(0, 1) +     # to add legend on plot
              theme_pubr() +
              theme(legend.position='none')
pl2

# Trial sequence as a function of the conflict load (not kept in final version)
pl3 <- ggplot(all_trials, aes(x = trial_id, y = conflict, color = Group)) +
              geom_line(aes(x = trial_id, y = conflict, group = subj_idx, color = Group), color = c("lightgrey"), linewidth = 0.5) +
              geom_point(alpha = 0.6, size = 2) +
              geom_hline(aes(yintercept = 0), size = 1) + 
              #geom_hline(aes(yintercept = 148.2328), size = 1, color = "lightblue") + 
              labs(y = 'Q Generous - Q Selfish', x = 'Trial n°') +
              scale_color_manual(values = c("darkgrey", "purple1", "orange1")) +
              #xlim(0, 113) +
              theme_pubr() +
              theme(legend.position='none')
pl3

mean(all_trials$conflict)
median(all_trials$conflict)

#===============#
# Decision as a function of conflict (sanity check not kept in manuscript)
#===============#

# Effect of conflict on Prob(gen)
GLM_full <- glmer(mod_choice ~ group * conflict_load + (1 | subj_idx),
                  data =all_trials,
                  family = binomial(link='logit'))
summary(GLM_full)
plot(allEffects(GLM_full))
(GLM_high <- pairs(emmeans(GLM_full, ~ group,type="mod_choice"))%>%test())

# Effect of conflict on RT
LMM_full <- lmer(RT ~ group * conflict_load + (1 | subj_idx),
                  data =all_trials)
summary(LMM_full)
plot(allEffects(LMM_full))
(GLM_high <- pairs(emmeans(LMM_full, ~ group,type="mod_choice"))%>%test())

# ============================================================================ #
# ---------------------------------------------------------------------------- #
# ============================================================================ #