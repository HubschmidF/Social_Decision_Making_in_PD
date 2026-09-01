//============================================================================//
// Hierarchical Bayesian Hyperbolic Discount Model
//============================================================================//

// Author: Fabrice Hubschmid
// Contact: hubschmidfabrice@gmail.com

// Value function: Hyperbolic (k)
// Link function: SoftMax (beta)

// Disclaimers / Context:

//  - Hierarchical implementation of the model of Sellito et al. (2021)
//    https://doi.org/10.1371/journal.pcbi.1007615
//  - Code structured as the models of the hBayesDM package (Ahn & al. 2017)
//  - using non-centered parametrization (Matt Trick)
//  - Parameter prior for the discount rate come from Ahn & al. (2017)
//    https://doi.org/10.1162/cpsy_a_00002

//============================================================================//

//-----------//
// Data 
//-----------//
data {
  int <lower=1> nS;                           // number of subjects
  int <lower=1> nT_max;                       // max number of trials
  int <lower=1> nT_subj[nS];                  // number of trial per subect
  int <lower=0> social_distance[nS,nT_max];   // social distance in a trial
  int <lower=0> amount_selfish[nS,nT_max];    // amount that can be won with the selfish choice
  int <lower=0> amount_generous[nS,nT_max];   // amount that can be won with the generous choice
  int <lower=0,upper=2> choice[nS,nT_max];    // stimulus choice (1 for generous and 2 for selfish)
  real <lower=0> rt[nS, nT_max];              // reaction time in a trial (not used)
  real <lower=0> min_rt[nS];                  // minimum reaction time for each subject (not used)
  real <lower=0> rt_bound;                    // global minumum bound for non-decision time (not used)
}
//-----------//
// Model Parameters
//-----------//
// Raw parameters
parameters {
  
  real k;                   // discount rate
  real <lower=0> k_sd;      // SD of discount rate
  real beta;                // inverse temperature
  real <lower=0> beta_sd;   // SD of inverse temperature
  
  // Individual level parameters
  vector [nS] k_s;          // individual raw discount rate
  vector [nS] beta_s;       // individual raw inverse temperature

}
//-----------//
// Transformed parameters
//-----------//
// Effective parameters used in the sampling statement
transformed parameters {
  vector <lower=0, upper=1> [nS] k_effective;
  vector <lower=0, upper=20> [nS] beta_effective;

  for (i in 1:nS) {
    k_effective[i]    = Phi_approx(k + k_sd * k_s[i]);
    beta_effective[i] = Phi_approx(beta + beta_sd * beta_s[i]) * 20;
  }
}
//-----------//
// Model
//-----------//
model {
  // Model priors for hyper parameters
  k    ~ normal(0,1);
  beta ~ normal(0,1);

  // Model priors for SD parameters
  k_sd    ~ normal(0,0.2);
  beta_sd ~ normal(0,0.2);
  
  // Model priors for individual level parameters 
  k_s    ~ normal(0,1);
  beta_s ~ normal(0,1);
  
  { // local section for fitting
  // iterates over subjects
    for (s in 1:nS) {
      
      int nT_s = nT_subj[s];
      vector[2] Q_values;
      
      // iterates over trials
      for (t in 1:nT_s) {
        
        // estimated choice values
        Q_values[1] = amount_generous[s,t] + (amount_generous[s,t] / (1 + k_effective[s] * social_distance[s,t]));
        Q_values[2] = amount_selfish[s,t];
        
        // sampling statement
        choice[s,t] ~ categorical_logit(beta_effective[s] * Q_values);
      }
    }
  }
}
//-----------//
// Model extracted quantities
//-----------//
generated quantities {
  // Hyper-parameters
  real <lower=0,upper=1> hyper_k;
  real <lower=0,upper=20> hyper_beta;
  // k in log space
  real log_k;
  
  // forward sampled quantities
  real log_lik[nS];
  real y_pred[nS,nT_max];
  real Q_gen[nS,nT_max];
  
  // Back transformed hyper parameters into interpretable space
  hyper_k = Phi_approx(k);
  log_k = log(hyper_k);
  hyper_beta = Phi_approx(beta)*20;
  
  // local section for forward sampled quantities 
     { 
      for (s in 1:nS) {
        
        int nT_s = nT_subj[s];
        vector[2] Q_vals;
        log_lik[s]=0;

       for (t in 1:nT_s) {
         
          // estimated action values
          Q_vals[1] = amount_generous[s,t] + (amount_generous[s,t] / (1 + k_effective[s] * social_distance[s,t]));
          Q_vals[2] = amount_selfish[s,t];
          Q_gen[s,t] = Q_vals[1];
          
          // log-likelyhood and model forward sampled predicted choice
          log_lik[s] += categorical_logit_lpmf(choice[s,t] | beta_effective[s] * Q_vals);
          y_pred[s,t] = categorical_logit_rng(beta_effective[s] * Q_vals);
       } 
     }
   }
}
