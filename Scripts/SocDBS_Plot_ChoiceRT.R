#==============================================================================# 
# SocialDBS: Plot Choice-Reactions
#==============================================================================#

# > FH

# Description:
#  - Script used for plotting the histogram of choices and RTs 
#  - not kept in final version of manuscript

#==============================================================================#

#===============#
# Path, Packages and data curation
#===============#

# Path and software packages
suppressMessages(suppressWarnings(library(rstudioapi))) # utility for path management
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
DF$log_rt <- log(DF$RT)

# Factorization
DF$stim_choice <- as.factor(DF$response) 
DF$group <- as.factor(DF$group) 
print(levels(DF$stim_choice))
levels(DF$stim_choice) <- c("miss", "missed", "selfish", "generous")
levels(DF$group) <- c("HC", "OFF", "ON")

# Exclude missed responses and the one weird trial that was logged
all_trials <- subset(DF, stim_choice %in% c("selfish","generous"))  
excluded_trials <- subset(DF, stim_choice %in% c("whatever","missed")) 
print(length(excluded_trials)) 

# Reward difference "rd"
all_trials$RD <- all_trials$selfish_amount - 75

# Re-organises the mod_choice so that 1 corresponds to generocity and 0 to selfish
all_trials$bin_choice <- as.factor(all_trials$response)
levels(all_trials$bin_choice)  <- c("selfish", "generous")   #
all_trials$bin_choice <- relevel(all_trials$bin_choice, ref = "selfish") #
levels(all_trials$bin_choice) <- c("0", "1")
all_trials$mod_choice <- as.integer(all_trials$bin_choice)
all_trials$mod_choice <- as.integer(all_trials$mod_choice) -1

#===============#
# Choice reaction time distributions
#===============#

# Flipping the reaction time for the selfish choice to be negative 
generous <- subset(all_trials, mod_choice %in% 1)
selfish  <- subset(all_trials, mod_choice %in% 0)
neg_selfish <- selfish
neg_selfish[sapply(neg_selfish, is.double)] <- neg_selfish[sapply(neg_selfish, is.double)] * -1
DFplot <- rbind(generous, neg_selfish)
DFplot$Group <- DFplot$group     # just for capitalization in the plot legend
levels(DFplot$group)

# Annotations
an_gen <- grobTree(textGrob("Generous", x=0.8,  y=0.8, hjust=0, gp=gpar(col="black", fontsize=12, fontface="bold")))
an_sel <- grobTree(textGrob("Selfish", x=0.1,  y=0.8, hjust=0, gp=gpar(col="black", fontsize=12, fontface="bold")))

# Histogram
pl1 <- ggplot(DFplot, aes(x = RT, color = Group, fill = Group)) +
              geom_histogram( bins =150, alpha = 0.1, position = 'identity') +
              #geom_density(alpha=0.1) +
              labs(y = 'Counts', x = 'Reaction Time (s)') +
              scale_color_manual(values=c("darkgrey", "purple1", "orange1")) +
              scale_fill_manual(values=c("darkgrey", "purple1", "orange1")) +
              theme_pubr() +
              theme(legend.position='bottom') +
              annotation_custom(an_gen) +
              annotation_custom(an_sel)

pl1

#==============================================================================#
#------------------------------------------------------------------------------#
#==============================================================================#