#==============================================================================# 
# SocialDBS: Single model fit
#==============================================================================#

# FH

# Description:
#  - fits one model to all three data samples (HC - ON - OFF)
#  - Saves model output for potential future use
#  - Adapt the path at which model outputs need to be saved below !!!

# Model
model_name <- c("Model_HMM_discount.stan")  # Model name

# MCMC parameterization
n_chs <- 4                           # n° of chains and cores used for sampling
iters <- 4000                        # n° of sampling iterations
burns <- 2000                        # n° of initial discarded samples
n_draws <- (iters - burns)* n_chs    # n° of posterior samples used for inference

# Groups to fit 
wantedgroup <- c(1:3)                # 1 = HC, 2 = DBS-ON, 3 = DBS-OFF

# Removal of outlying trials or not
remove_outlying_trials <- FALSE     # if activated will remove the trials
                                    # in the 5% pre and tail ends of the 
                                    # RT distribution

#==============================================================================#

#===============#
# Path, packages and curation
#===============#

# Software packages
suppressMessages(suppressWarnings(library(rstudioapi))) # path management
suppressMessages(suppressWarnings(library(tictoc)))     # utility to estimate sampling times
suppressMessages(suppressWarnings(library(dplyr)))      # utility for data frame manipulations
suppressMessages(suppressWarnings(library(plyr)))       # utility for data frame manipulations
suppressMessages(suppressWarnings(library(writexl)))    # for saving matrix of posterior draws
suppressMessages(suppressWarnings(library(tidyr)))      # backbone R software package
suppressMessages(suppressWarnings(library(rstan)))      # for interfacing stan models
suppressMessages(suppressWarnings(library(loo)))        # functions for Bayesian Leave-One-Out validation
suppressMessages(suppressWarnings(library(HDInterval))) # highest density interval estimations

# Path
analysis_path = getActiveDocumentContext()$path
setwd(dirname(analysis_path))
analysis_dir <- getwd()

# Data file
data_file <- file.path(analysis_dir, "Task_Data_Modelling.csv")
data <- read.csv(data_file, header = T, sep = ",")
DF <- as.data.frame(data)

DF$RT <- DF$rt/1000
DF$stim_choice <- as.factor(DF$response)
print(levels(DF$stim_choice))

levels(DF$stim_choice) <- c("miss", "missed", "selfish", "generous")

# Exclude missed responses and the one weird trial that was logged as -1
all_trials <- subset(DF, stim_choice %in% c("selfish","generous"))
excluded_trials <- subset(DF, stim_choice %in% c("miss","missed"))
print(nrow(excluded_trials))

# Utility for manipulation of posterior draws
group_HC <- rep(1,n_draws)
group_ON <- rep(2,n_draws)
group_OFF <- rep(3,n_draws)
#
Group <- c(group_HC, group_ON, group_OFF)
Group <- as.factor(Group)
levels(Group) <- c('HC', 'ON', 'OFF')

#===============#
# Fitting
#===============#

# initial pre-compilation
model_path <- file.path(analysis_dir, model_name)
rstan_options("auto_write" = T) 
model <- rstan::stan_model(model_path)

for (thatgroup in wantedgroup) {
  
  print(thatgroup)
  
  # Gets relevant data
  if (thatgroup == 1) {    
    hec_trials <- subset(all_trials, group %in% c("NO"))
  } else if (thatgroup == 2) {
    hec_trials <- subset(all_trials, group %in% c("ON"))
  } else if (thatgroup == 3) {
    hec_trials <- subset(all_trials, group %in% c("OFF"))
  }
  
  if (remove_outlying_trials == TRUE) {    
    # removes outlying rt values                                                                              
    RT_Bounds <- hdi(hec_trials$RT, credMass = .95)         # gets the 95% HDI in the RT distribution
    low <- RT_Bounds[1]                                     # corresponds to lower bound
    hig <- RT_Bounds[2]                                     # corresponds to upper bound
    hec_trials <- subset(hec_trials, RT >= low & RT <= hig) # removes trials with outlying values
  }

  # Organizes the columns needed in the dataset
  hec_trials$bin_choice <- as.factor(hec_trials$response)
  levels(hec_trials$bin_choice)  <- c("selfish", "generous")   #
  hec_trials$bin_choice <- relevel(hec_trials$bin_choice, ref = "generous") #
  levels(hec_trials$bin_choice) <- c("0", "1")
  hec_trials$mod_choice <- as.integer(hec_trials$bin_choice)
  hec_trials$mod_choice <- as.integer(hec_trials$mod_choice)
  # generous amount 
  hec_trials$generous_amount <- rep(75, length(hec_trials$trial_id))
  hec_trials$generous_amount <- as.integer(hec_trials$generous_amount)
  
  # Creates the data list passed onto stan 
  n_subj <-                                 
    length(unique(hec_trials$Part_idx))
  
  n_trials_subjs <-                          
    hec_trials %>%
    group_by(Part_idx) %>%
    dplyr::summarise(n=n()) %>%
    pull(n)
  
  n_trials_max <- max(n_trials_subjs)        
  
  choices <-                                   
    hec_trials %>%
    dplyr::select(Part_idx, mod_choice) %>%    
    group_by(Part_idx) %>%
    dplyr::mutate(trial_id= 1:n()) %>%
    ungroup()%>%
    pivot_wider(names_from = trial_id, values_from = mod_choice) %>%
    dplyr::select(-Part_idx) %>%
    replace(is.na(.), 0) %>%
    simplify2array()
  
  distances <-                                   
    hec_trials %>%
    dplyr::select(Part_idx, social_distance) %>% 
    group_by(Part_idx) %>%
    dplyr::mutate(trial_id= 1:n()) %>%
    ungroup()%>%
    pivot_wider(names_from = trial_id, values_from = social_distance) %>%
    dplyr::select(-Part_idx) %>%
    replace(is.na(.), 0) %>%
    simplify2array()
  
  shell_fishes <-                                   
    hec_trials %>%
    dplyr::select(Part_idx, selfish_amount) %>% 
    group_by(Part_idx) %>%
    dplyr::mutate(trial_id= 1:n()) %>%
    ungroup()%>%
    pivot_wider(names_from = trial_id, values_from = selfish_amount) %>%
    dplyr::select(-Part_idx) %>%
    replace(is.na(.), 0) %>%
    simplify2array()
  
  genes_n_roses <-                                   
    hec_trials %>%
    dplyr::select(Part_idx, generous_amount) %>% 
    group_by(Part_idx) %>%
    dplyr::mutate(trial_id= 1:n()) %>%
    ungroup()%>%
    pivot_wider(names_from = trial_id, values_from = generous_amount) %>%
    dplyr::select(-Part_idx) %>%
    replace(is.na(.), 0) %>%
    simplify2array()
  
  reactions <-                                   
    hec_trials %>%
    dplyr::select(Part_idx, RT) %>%    
    group_by(Part_idx) %>%
    dplyr::mutate(trial_id= 1:n()) %>%
    ungroup()%>%
    pivot_wider(names_from = trial_id, values_from = RT) %>%
    dplyr::select(-Part_idx) %>%
    replace(is.na(.), 0) %>%
    simplify2array()
  
  # minimum reaction times per subject
  vec <- as.integer(c())
  
  if (thatgroup == 1) {   # Healthy
    
    ids <- c(23:27,29:40,42)           # Participant 28,and 41 excluded from those analyses  

    for (subject in ids) {
      sub_string <- as.character(subject)
      that_one_subject <- subset(hec_trials, Part_idx %in% c(sub_string))
      vec <- cbind(vec, min(that_one_subject$RT))
    }
    
  } else if (thatgroup == 2) {    # Stimulation ON
    
    ids <- c(1:2,4:22)
    
    for (subject in ids) {
      sub_string <- as.character(subject)
      that_one_subject <- subset(hec_trials, Part_idx %in% c(sub_string))
      vec <- cbind(vec, min(that_one_subject$RT))
    }
    
  } else if (thatgroup == 3) {     # Stimulation OFF
    
    ids <- c(1:2,4:22)
    
    for (subject in ids) {
      sub_string <- as.character(subject)
      that_one_subject <- subset(hec_trials, Part_idx %in% c(sub_string))
      vec <- cbind(vec, min(that_one_subject$RT))
    }
  }
  
  MinRTS <- vec[1,]    # recovers vector with min RT for each subject
  # Minimum reaction in dataset   
  RT_min <- 0.000000000001 #(needs a global minimum bound that is not 0)   
  
  # Final List
  stan_list <- 
    list(
      nS = n_subj,
      nT_max = n_trials_max,
      nT_subj = n_trials_subjs,
      social_distance = distances,
      amount_selfish = shell_fishes,    
      amount_generous = genes_n_roses,  
      choice = choices,
      rt = reactions,
      min_rt = MinRTS,
      rt_bound = RT_min
    )

  # Fits the model to all groups
  
  if (thatgroup == 1) {
    
    Model_Healthy <- 
      sampling(
        model,
        data = stan_list,
        chains = n_chs,
        warmup = burns,
        iter = iters,
        cores = n_chs,
        refresh = 100
      )
    
  } else if (thatgroup == 2) {
    
    Model_PD_ON <- 
      sampling(
        model,
        data = stan_list,
        chains = n_chs,       
        warmup = burns,     
        iter = iters,
        cores = n_chs,
        refresh = 100
      )
    
  } else if (thatgroup == 3) {
    
    Model_PD_OFF <- 
      sampling(
        model,
        data = stan_list,
        chains = n_chs,
        warmup = burns,
        iter = iters,
        cores = n_chs,
        refresh = 100
      )
  }
}
print(warnings())

# Matrices of posterior draws
DrawsHC <- as.data.frame(Model_Healthy)
DrawsON <- as.data.frame(Model_PD_ON)
DrawsOFF <- as.data.frame(Model_PD_OFF)

#===============#
# Saving posterior draws for further use
#===============#

# Healthy data
write_xlsx(DrawsHC,
           "PATH\\TO\\YOUR\\ANALYSIS\\HMM_posterior_draws_healthy.xlsx",  # ADAPT PATH!
           col_names = TRUE)

# Stimulation ON data
write_xlsx(DrawsON,
           "PATH\\TO\\YOUR\\ANALYSIS\\HMM_posterior_draws_DBS_ON.xlsx",   # ADAPT PATH!
           col_names = TRUE)

# Stimulation OFF data
write_xlsx(DrawsOFF,
           "PATH\\TO\\YOUR\\ANALYSIS\\HMM_posterior_draws_DBS_OFF.xlsx",  # ADAPT PATH!
           col_names = TRUE)

#===============#
# Diagnostics
#===============#

if (model_name == "SM_discount.stan") {
  params <- c("log_k", "hyper_beta")
} else if (model_name == "DDM_discount_nosim.stan") {
  params <- c("hyper_a", "hyper_ter", "hyper_vcoeff", "log_k")
} else if (model_name == "HMM_discount_alt.stan") {
  params <- c("log_k","hyper_beta")
}

# Pairs Plots
pairs(Model_Healthy, pars = params)
pairs(Model_PD_ON, pars = params)
pairs(Model_PD_OFF, pars = params)

# chainplots
traceplot(Model_Healthy ,pars = params, inc_warmup = FALSE)
traceplot(Model_PD_ON ,pars = params, inc_warmup = FALSE)
traceplot(Model_PD_OFF ,pars = params, inc_warmup = FALSE)

#===============#
# Model Summaries
#===============#

summ_HC <- summary(Model_Healthy)
summ_ON <- summary(Model_PD_ON)
summ_OF <- summary(Model_PD_OFF)
#
overview_HC <- summ_HC$summary
overviewDF_HC <- as.data.frame(overview_HC)
#
overview_ON <- summ_ON$summary
overviewDF_ON <- as.data.frame(overview_ON)
#
overview_OF <- summ_OF$summary
overviewDF_OF <- as.data.frame(overview_OF)

# Save to xcel sheet

# HC
colHC <- row.names(overviewDF_HC)
overviewDF_HC$parameters <- colHC
write_xlsx(overviewDF_HC,
           "PATH\\TO\\YOUR\\ANALYSIS\\HMM_Model_Summary_HC.xlsx",  # ADAPT PATH!
           col_names = TRUE)

# ON
colON <- row.names(overviewDF_ON)
overviewDF_ON$parameters <- colON
write_xlsx(overviewDF_ON,
           "PATH\\TO\\YOUR\\ANALYSIS\\HMM_Model_Summary_ON.xlsx",  # ADAPT PATH!
           col_names = TRUE)

# OFF
colOF <- row.names(overviewDF_OF)
overviewDF_OF$parameters <- colOF
write_xlsx(overviewDF_OF,
           "PATH\\TO\\YOUR\\ANALYSIS\\HMM_Model_Summary_OFF.xlsx", # ADAPT PATH!
           col_names = TRUE)

#==============================================================================#
#------------------------------------------------------------------------------#
#==============================================================================#