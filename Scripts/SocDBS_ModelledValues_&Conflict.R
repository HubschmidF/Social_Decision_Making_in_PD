#==============================================================================# 
# SocialDBS: Model extracted quantities: Q values and conflict
#==============================================================================#

# > FH

# Description:
#  - plots modeled Q values and highlights modeled conflict
#  - Situations of increased conflict (GLM / LMM)

#==============================================================================#

#===============#
# Path, packages and data
#===============#

suppressMessages(suppressWarnings(library(rstudioapi))) # path management
suppressMessages(suppressWarnings(library(tictoc)))     # utility to estimate sampling times
suppressMessages(suppressWarnings(library(dplyr)))      # utility for data frame manipulations
suppressMessages(suppressWarnings(library(plyr)))       # utility for data frame manipulations
suppressMessages(suppressWarnings(library(ggplot2)))  
suppressMessages(suppressWarnings(library(ggpubr)))  
suppressMessages(suppressWarnings(library(HDInterval))) 
suppressMessages(suppressWarnings(library(lme4)))       #used for fitting of LMMs
suppressMessages(suppressWarnings(library(effects)))    #used for quick visualisation only
suppressMessages(suppressWarnings(library(lmerTest)))   #used for significance of fixed effects
suppressMessages(suppressWarnings(library(emmeans)))    #used for contrast post-hoc comparisons

# Path
analysis_path = getActiveDocumentContext()$path
setwd(dirname(analysis_path))
analysis_dir <- getwd()

# Data file
data_file <- file.path(analysis_dir, "Task_Data_Modelling_Added_Qvalues.csv")
data <- read.csv(data_file, header = T, sep = ",")
DF <- as.data.frame(data)

DF$RT <- DF$rt/1000
DF$stim_choice <- as.factor(DF$response)
print(levels(DF$stim_choice))
levels(DF$stim_choice) <- c("selfish", "generous")

# for utility
all_trials <- subset(DF, stim_choice %in% c("selfish","generous"))
all_trials$log_q_generous <- log(all_trials$Q_gen)
all_trials$log_selfish <-  log(all_trials$selfish_amount)

all_trials$conflict <- all_trials$Q_gen - all_trials$selfish_amount
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
pl1 <- ggplot(one_ammount, aes(x = social_distance, y = Q_gen, color = Group)) +
             geom_line(aes(x = social_distance, y =  Q_gen, group = log_idx, color = Part_idx), color = c("lightgrey"), linewidth = 0.5) +
             geom_point(alpha = 0.9, size = 3) +
             labs(y = 'Q generous', x = 'Social Distance') +
             scale_color_manual(values = c("darkgrey", "purple1", "orange1")) +
             #ylim(0,500) +
             theme_pubr() +
             theme(legend.position='none')
pl1

# Distribution of decision conflict (not kept in final version)
pl2 <- ggplot(all_trials, aes(x = conflict_load, color = Group, fill = Group)) + 
              geom_density(alpha = 0.2) +
              geom_vline(xintercept = 0, color = 'black', size = 1, linetype = "dashed") +
              scale_color_manual(values = c("darkgrey", "purple1", "orange1")) +
              scale_fill_manual(values = c("darkgrey", "purple1", "orange1")) +
              labs(y = 'Density', x = 'Decisional Conflict (a.u)') +
              #ylim(0, 1) +     # to add legend on plot
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
# Decision as a function of conflict (sanity check not kept)
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

#===============#
# Engagement states over sequence
#===============#

new_trials <- subset(all_trials, subj_idx %in% c(7,13,15)) # are all HCs
new_trials$subj <- as.factor(new_trials$subj_idx)

pl4 <- ggplot(new_trials, aes(x = trial_id, y = P_eng, color = subj)) +
              geom_line(aes(x = trial_id, y = P_eng, group = subj, color = Group), color = c("darkgrey"), linewidth = 0.8) +
              geom_point(alpha = 0.8, size = 2) +
              geom_hline(aes(yintercept = 0.5), size = 0.5, color = "red", linetype = "dashed") + 
              labs(y = 'Estimated P(engaged)', x = 'Trial n°') +
              scale_color_manual(values = c("blue", "darkblue", "lightblue")) +
              theme_pubr() +
              theme(legend.position='none')
pl4

# ============= #
# Some GLMs
# ============= #

# implementational check
all_trials$bin_choice <- as.factor(all_trials$response)
levels(all_trials$bin_choice)  <- c("selfish", "generous")   #
all_trials$bin_choice <- relevel(all_trials$bin_choice, ref = "selfish") #
levels(all_trials$bin_choice) <- c("0", "1")
all_trials$mod_choice <- as.integer(all_trials$bin_choice)
all_trials$mod_choice <- as.integer(all_trials$mod_choice) -1

# trial pools
eng_trials <- subset(all_trials, P_eng > 0.5)
ran_trials <- subset(all_trials, P_eng < 0.5)
eng_trials$bin_eng <- 1
ran_trials$bin_eng <- 0

# Choice probs in engaged state
GLM_eng <- glmer(mod_choice ~ group + social_distance + RD + (1 | Part_idx),
                    data = eng_trials ,
                    family = binomial(link='logit'))
summary(GLM_eng)

# Choice probs in random state (effects shrink to zero)
GLM_ran <- glmer(mod_choice ~ group + social_distance + RD + (1 | Part_idx),
                 data = ran_trials ,
                 family = binomial(link='logit'))
summary(GLM_ran)

# Reconstructed dataset
recon_trials <- rbind(eng_trials,ran_trials)

GLM_eng <- glmer(bin_eng ~ group + (1 | Part_idx),
                 data = recon_trials ,
                 family = binomial(link='logit'))
summary(GLM_eng)
confint(GLM_eng)

LMM_eng <- lmer(P_eng ~ group + (1 | Part_idx),
                 data = recon_trials)
summary(LMM_eng)
confint(LMM_eng)

#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
#Formula: P_eng ~ group + (1 | Part_idx)
#   Data: recon_trials

#REML criterion at convergence: -19132.7

#Scaled residuals: 
#    Min      1Q  Median      3Q     Max 
#-4.4115 -0.5079 -0.0044  0.4151  5.1475 

#Random effects:
# Groups   Name        Variance Std.Dev.
# Part_idx (Intercept) 0.042153 0.20531 
# Residual             0.002678 0.05175 
#Number of obs: 6306, groups:  Part_idx, 39

#Fixed effects:
#            Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)  0.50410    0.04841 37.00968  10.414  1.5e-12 ***
#groupOFF    -0.21328    0.06597 37.00993  -3.233  0.00258 ** 
#groupON     -0.17905    0.06597 37.01014  -2.714  0.01003 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#         (Intr) grpOFF
#groupOFF -0.734       
#groupON  -0.734  1.000

#2.5 %      97.5 %
#  .sig01       0.16262088  0.25414811
#.sigma       0.05085138  0.05266346
#(Intercept)  0.40936605  0.59883400
#groupOFF    -0.34238215 -0.08418047
#groupON     -0.30815350 -0.0499514

(LMM_group <- pairs(emmeans(LMM_eng, ~ group,type="P_eng"))%>%test())

#contrast estimate      SE  df z.ratio p.value
#NO - OFF   0.2133 0.06600 Inf   3.233  0.0035
#NO - ON    0.1791 0.06600 Inf   2.714  0.0182
#OFF - ON  -0.0342 0.00156 Inf -21.910  <.0001

# ============= #
# Boxplot
# ============= #

all_trials$group <- as.factor(all_trials$group)
levels(all_trials$group) <- c("HC", "OFF", "ON")

pl5 <- ggplot(all_trials, aes(x = group, y = P_eng, fill = group)) +
              geom_point(position = position_jitter(width = 0.05), size = 1.5, alpha = 0.005) +
              geom_boxplot(width = 0.3, notch = TRUE, fill = c("darkgrey", "purple1", "orange1"), alpha = 0.5) +
              theme_pubr() +
              labs(x = "Group", y = "Estimated P(engaged)") +
              theme(legend.position = "none")
pl5

# ============================================================================ #
# ---------------------------------------------------------------------------- #
# ============================================================================ #