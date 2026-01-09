#==============================================================================# 
# SocialDBS: Comparison of models with a DDM mapping function
#==============================================================================#

# FH

# Description:
#  - fits all models to the three experimental groups
#  - insures adequate convergence of MCMC chains
#  - recovers matrix of log-likelihoods
#  - compares model predictive performance based on LOO

# Models to fit
model_names <- c("Model_HDDM_final.stan",      # Baseline HDDM model
                 "Model_HDDDM_final.stan",     # Hyperbolic discount DDM model
                 "Model_EWHDDDM_final.stan")   # Empathy weighted hyperbolic discount DDM model

# Groups to fit it to
wantedgroup <- c(1:3)  # fits to all experimental groups (1 = HC, 2 = DBS-ON, 3 = DBS-OFF)

# MCMC parametrization
n_chs <- 4                          # n° of chains (and cores) used for sampling
iters <- 4000                       # n° of sampling iterations
burns <- 2000                       # n° of initial discarded samples

# Removal of outlying trials or not
remove_outlying_trials <- TRUE       # if activated will remove the trials
                                     # in the 5% upper and lower tail ends of
                                     # the RT distribution

#==============================================================================#

#===============#
# Path, packages and curation
#===============#

# Software packages (probably not all in use here)
suppressMessages(suppressWarnings(library(rstudioapi)))     # path management
suppressMessages(suppressWarnings(library(tictoc)))         # utility to estimate sampling times
suppressMessages(suppressWarnings(library(dplyr)))          # utility for data frame manipulations
suppressMessages(suppressWarnings(library(plyr)))           # utility for data frame manipulations
suppressMessages(suppressWarnings(library(tidyr)))          # backbone R software package
suppressMessages(suppressWarnings(library(rstan)))          # for interfacing STAN code
suppressMessages(suppressWarnings(library(loo)))            # functions for Bayesian Leave-One-Out validation
suppressMessages(suppressWarnings(library(HDInterval)))     # functions for highest density intervals

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

# Exclude missed responses
all_trials <- subset(DF, stim_choice %in% c("selfish","generous"))
excluded_trials <- subset(DF, stim_choice %in% c("miss","missed"))
print(nrow(excluded_trials))

#===============#
# Iterative fitting
#===============#

# initial pre-compilation of all models to save time
tic()                                                   
for (c in 1:length(model_names)) {                        
  model_path <- file.path(analysis_dir, model_names[c])
  rstan_options("auto_write" = T) 
  model <- rstan::stan_model(model_path)
}
toc()

tic()
for (n in 1:length(model_names)) {   
  
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
      RT_Bounds <- hdi(hec_trials$RT, credMass = .95)           
      low <- RT_Bounds[1]                                      
      hig <- RT_Bounds[2]                                       
      hec_trials <- subset(hec_trials, RT >= low & RT <= hig)
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
      dplyr::summarise(n=n())  %>%
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
      
      ids <- c(1:2,4:22)                 # Participant 3 excluded from those analyses  
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
    # Minimum reaction in dataset (needs a global minimum bound that is not 0)   
    RT_min <- 0.000000000001
    
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
    
    # Getting the right model
    model_path <- file.path(analysis_dir, model_names[n])
    rstan_options("auto_write" = T) 
    model <- rstan::stan_model(model_path)
    
    # Fitting statement
    Model <- 
      sampling(
        model,
        data = stan_list,
        chains = n_chs,       
        warmup = burns,     
        iter = iters,
        init = 0,
        cores = n_chs,
        refresh = 100
      )
    
    # information recovery
    if      (n == 1 & thatgroup == 1) {  
      
      # Gets Log-likelyhood matrix
      LL_m1_HC = loo::extract_log_lik(Model, 'log_lik')                       
      print("1/9 now done")                                                       
      
      # Estimates chain convergence
      summ <- summary(Model)
      overview <- summ$summary
      overviewDF <- as.data.frame(overview)
      Rhat <- na.omit(overviewDF$Rhat)
      print(min(Rhat))
      print(max(Rhat))
      if (max(Rhat) > 1.1) {
        stop()   # clause only works if ran in windows
      }
      
      # Clears the workspace of unecessery objects
      rm("summ")
      rm("overview")
      rm("overviewDF")
      rm("Rhat")
    }
    else if (n == 1 & thatgroup == 2) {   
      
      # Gets Log-likelyhood matrix
      LL_m1_ON = loo::extract_log_lik(Model, 'log_lik')
      print("2/9 now done")
      
      # Estimates chain convergence
      summ <- summary(Model)
      overview <- summ$summary
      overviewDF <- as.data.frame(overview)
      Rhat <- na.omit(overviewDF$Rhat)
      print(min(Rhat))
      print(max(Rhat))
      if (max(Rhat) > 1.1) {
        stop()
      }
      
      # Clears the workspace of unecessery objects
      rm("summ")
      rm("overview")
      rm("overviewDF")
      rm("Rhat")
    }
    else if (n == 1 & thatgroup == 3) {   
      
      # Gets Log-likelyhood matrix
      LL_m1_OFF = loo::extract_log_lik(Model, 'log_lik')
      print("3/9 now done")
      
      # Estimates chain convergence
      summ <- summary(Model)
      overview <- summ$summary
      overviewDF <- as.data.frame(overview)
      Rhat <- na.omit(overviewDF$Rhat)
      print(min(Rhat))
      print(max(Rhat))
      if (max(Rhat) > 1.1) {
        stop()
      }
      
      # Clears the workspace of unecessery objects
      rm("summ")
      rm("overview")
      rm("overviewDF")
      rm("Rhat")
    }
    else if (n == 2 & thatgroup == 1) {   
      
      # Gets Log-likelyhood matrix
      LL_m2_HC = loo::extract_log_lik(Model, 'log_lik')
      print("4/9 now done")
      
      # Estimates chain convergence
      summ <- summary(Model)
      overview <- summ$summary
      overviewDF <- as.data.frame(overview)
      Rhat <- na.omit(overviewDF$Rhat)
      print(min(Rhat))
      print(max(Rhat))
      if (max(Rhat) > 1.1) {
        stop()
      }
      
      # Clears the workspace of unecessery objects
      rm("summ")
      rm("overview")
      rm("overviewDF")
      rm("Rhat")
    }
    else if (n == 2 & thatgroup == 2) {   
      
      # Gets Log-likelyhood matrix
      LL_m2_ON = loo::extract_log_lik(Model, 'log_lik')
      print("5/9 now done")
      
      # Estimates chain convergence
      summ <- summary(Model)
      overview <- summ$summary
      overviewDF <- as.data.frame(overview)
      Rhat <- na.omit(overviewDF$Rhat)
      print(min(Rhat))
      print(max(Rhat))
      if (max(Rhat) > 1.1) {
        stop()
      }
      
      # Clears the work space of unnecessary objects
      rm("summ")
      rm("overview")
      rm("overviewDF")
      rm("Rhat")
    }
    else if (n == 2 & thatgroup == 3) {   
      
      # Gets Log-likelyhood matrix
      LL_m2_OFF = loo::extract_log_lik(Model, 'log_lik')
      print("6/9 now done")
      
      # Estimates chain convergence
      summ <- summary(Model)
      overview <- summ$summary
      overviewDF <- as.data.frame(overview)
      Rhat <- na.omit(overviewDF$Rhat)
      print(min(Rhat))
      print(max(Rhat))
      if (max(Rhat) > 1.1) {
        stop()
      }
      
      # Clears the workspace of unecessery objects
      rm("summ")
      rm("overview")
      rm("overviewDF")
      rm("Rhat")
    }
    else if (n == 3 & thatgroup == 1) {   
      
      # Gets Log-likelyhood matrix
      LL_m3_HC = loo::extract_log_lik(Model, 'log_lik')
      print("7/9 now done")
      
      # Estimates chain convergence
      summ <- summary(Model)
      overview <- summ$summary
      overviewDF <- as.data.frame(overview)
      Rhat <- na.omit(overviewDF$Rhat)
      print(min(Rhat))
      print(max(Rhat))
      if (max(Rhat) > 1.1) {
        stop()
      }
      
      # Clears the workspace of unecessery objects
      rm("summ")
      rm("overview")
      rm("overviewDF")
      rm("Rhat")
    }
    else if (n == 3 & thatgroup == 2) {
      
      # Gets Log-likelyhood matrix
      LL_m3_ON = loo::extract_log_lik(Model, 'log_lik')
      print("8/9 now done")
      
      # Estimates chain convergence
      summ <- summary(Model)
      overview <- summ$summary
      overviewDF <- as.data.frame(overview)
      Rhat <- na.omit(overviewDF$Rhat)
      print(min(Rhat))
      print(max(Rhat))
      if (max(Rhat) > 1.1) {
        stop()
      }
      
      # Clears the workspace of unecessery objects
      rm("summ")
      rm("overview")
      rm("overviewDF")
      rm("Rhat")
    }
    else if (n == 3 & thatgroup == 3) {   
      
      # Gets Log-likelyhood matrix
      LL_m3_OFF = loo::extract_log_lik(Model, 'log_lik')
      print("9/9 now done")
      
      # Estimates chain convergence
      summ <- summary(Model)
      overview <- summ$summary
      overviewDF <- as.data.frame(overview)
      Rhat <- na.omit(overviewDF$Rhat)
      print(min(Rhat))
      print(max(Rhat))
      if (max(Rhat) > 1.1) {
        stop()
      }
      
      # Clears the workspace of unecessery objects
      rm("summ")
      rm("overview")
      rm("overviewDF")
      rm("Rhat")
    }
  }        # closes group loop
}        # closes model loop
toc()

#===============#
# Model Log likelihoods
#===============#

# Baseline HDDM
loo_HC_HDDM <- loo(LL_m1_HC)
loo_ON_HDDM <- loo(LL_m1_ON)
loo_OFF_HDDM <- loo(LL_m1_OFF)

# Discount HDDM
loo_HC_HypDDM <- loo(LL_m2_HC)
loo_ON_HypDDM <- loo(LL_m2_ON)
loo_OFF_HypDDM <- loo(LL_m2_OFF)

# Weighted-Discount HDDM
loo_HC_WeiDDM <- loo(LL_m3_HC)
loo_ON_WeiDDM <- loo(LL_m3_ON)
loo_OFF_WeiDDM <- loo(LL_m3_OFF)

#===============#
# Comparisons of goodness of fit
#===============#

# In Healthy individuals
comp_HC <- loo_compare(loo_HC_HDDM, loo_HC_HypDDM, loo_HC_WeiDDM)
HC_elpd_df <- as.data.frame(comp_HC)
HC_elpd_df$sig <- HC_elpd_df$elpd_diff / HC_elpd_df$se_diff
print(HC_elpd_df)

# In Patients with DBS-ON
comp_ON <- loo_compare(loo_ON_HDDM, loo_ON_HypDDM, loo_ON_WeiDDM)
ON_elpd_df <- as.data.frame(comp_ON)
ON_elpd_df$sig <- ON_elpd_df$elpd_diff / ON_elpd_df$se_diff
print(ON_elpd_df)

# In Patients with DBS-OFF
comp_OFF <- loo_compare(loo_OFF_HDDM, loo_OFF_HypDDM, loo_OFF_WeiDDM)
OFF_elpd_df <- as.data.frame(comp_OFF)
OFF_elpd_df$sig <- OFF_elpd_df$elpd_diff / OFF_elpd_df$se_diff
print(OFF_elpd_df)

#==============================================================================#
#------------------------------------------------------------------------------#
#==============================================================================#