//============================================================================//
// Hierarchical discount drift diffusion model
//============================================================================//

// Author: Fabrice Hubschmid
// Contact: hubschmidfabrice@gmail.com

// Value function: Hyperbolic (k)
// Link function: DDM (a, ter, vcoeff) (fixed bias (0.5))

// Disclaimers / Context:

//  - uses non-centered parametrization (Matt Trick)
//  - Corresponds approx. to the model of Peters & D'Esposito (2020)
//    https://doi.org/10.1371/journal.pcbi.1007615
//  - Structured as the models of the hBayesDM package (Ahn & al. 2017)
//  - Parameter prior for the discount come from Ahn & al. (2017)
//    https://doi.org/10.1162/cpsy_a_00002
//  - Parameter priors for the DDM come from Matzke & Wagenmaker. (2009)
//    https://doi.org/10.3758/PBR.16.5.798 
//  - Adapted some code snippets from Desch & al. (2023)
//    https://doi.org/10.7554/eLife.81436

//============================================================================//

// INFO: in order for the script to simulate based on wiener first passage time
// distributions, you need to insert the function of  that can be found online
// @: 
// https://discourse.mc-stan.org/t/alternative-to-wiener-rng-in-stan/28443/11

// -- Here! --

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
  real <lower=0> rt[nS, nT_max];              // reaction time in a trial
  real <lower=0> min_rt[nS];                  // minimum reaction time for each subject
  real <lower=0> rt_bound;                    // global minumum bound for non-decision time
}
transformed data {
  // Normalization of monetary ammounts: 
  // this needed to keep vcoeff away from 0 as the sampler struggles with v = 0
  real a_selfish_norm[nS,nT_max];
  real a_generous_norm[nS,nT_max];
  
  for (s in 1:nS) {
    for (t in 1:nT_subj[s]) {
      a_selfish_norm[s,t] = amount_selfish[s,t] * 1.0 / amount_generous[s,t];
      a_generous_norm[s,t] = 1;
    }
  }
}
//-----------//
// Model Parameters
//-----------//
// Raw parameters
parameters {
  
  real k;                   // discount rate
  real <lower=0> k_sd;      // SD of discount rate
  real a;                   // boundary separation
  real <lower=0> a_sd;      // SD of boundary separation
  real ter;                 // non-decision times
  real <lower=0> ter_sd;    // SD of non-decision times
  real vcoeff;              // drift rate (scale factor)
  real <lower=0> vcoeff_sd; // SD of drift rate
  
  // Individual level parameters
  vector [nS] k_s;          // individual raw discount rate
  vector [nS] a_s;          // individual raw boundaries
  vector [nS] ter_s;        // individual raw non-decisions
  vector [nS] vcoeff_s;     // individual raw drift rate coefficient
}
//-----------//
// Transformed parameters
//-----------//
// Effective parameters used in the sampling statement
transformed parameters {
  vector <lower=0, upper=1> [nS] k_effective;
  vector <lower=0> [nS] a_effective;
  vector <lower=rt_bound, upper=max(min_rt)> [nS] ter_effective;
  vector <lower=0> [nS] vcoeff_effective;
  
  for (i in 1:nS) {
    k_effective[i]      = Phi_approx(k + k_sd * k_s[i]);
    a_effective[i]      = exp(a + a_sd * a_s[i]);
    ter_effective[i]    = Phi_approx(ter + ter_sd * ter_s[i]) * (min_rt[i] - rt_bound) + rt_bound;
    vcoeff_effective[i] = exp(vcoeff + vcoeff_sd * vcoeff_s[i]);
  }
}
//-----------//
// Model
//-----------//
model {
  // Model priors for hyper parameters
  k      ~ normal(0,1);
  a      ~ normal(1,1);
  ter    ~ normal(0.435,0.12);  
  vcoeff ~ normal(1.5, 1);

  // Model priors for SD parameters
  k_sd      ~ normal(0,0.2);
  a_sd      ~ normal(0,0.2);
  ter_sd    ~ normal(0,0.2);
  vcoeff_sd ~ normal(0,0.2);
  
  // Model priors for individual level parameters 
  k_s      ~ normal(0,1);
  a_s      ~ normal(0,1);
  ter_s    ~ normal(0,1);
  vcoeff_s ~ normal(0,1);
  
  { // local section for fitting
  // iterates over subjects
    for (s in 1:nS) {
      
      int nT_s = nT_subj[s];
      real Q_selfish;
      real Q_generous;
      real Drift;
      
      // iterates over trials
      for (t in 1:nT_s) {
        
        // estimated choice values
        Q_selfish = a_selfish_norm[s,t];                                                                                                     
        Q_generous = a_generous_norm[s,t] + (a_generous_norm[s,t] / (1 + k_effective[s] * social_distance[s,t]));
        
        // upper or lower bound choice
        if (choice[s,t]==1) {                                                                      
          Drift = vcoeff_effective[s] * (Q_generous - Q_selfish);
          
        } else { 
          Drift = vcoeff_effective[s] * (Q_selfish - Q_generous);
        }
        
        // sampling statement
        rt[s,t] ~ wiener(a_effective[s], ter_effective[s], 0.5, Drift);
      }
    }
  }
}
//-----------//
// Model extracted quantities
//-----------//
generated quantities {
  // Hyper-parameters
  real <lower=0, upper=1> hyper_k;
  real <lower=0> hyper_a;
  real <lower=rt_bound, upper=max(min_rt)> hyper_ter;
  real <lower=0> hyper_vcoeff;
  // k in log space
  real log_k;
  
  // Log likelihood
  real log_lik[nS];
  
  // Forward sampled quantities
  real Q_selfish[nS,nT_max];
  real Q_generous[nS,nT_max];
  real y_pred[nS,nT_max];
  real rt_pred[nS, nT_max];
  
  // Back transformed hyper parameters into interpretable space
  hyper_k = Phi_approx(k);
  log_k = log(hyper_k);
  hyper_a = exp(a);
  hyper_ter = Phi_approx(ter) * (mean(min_rt) - rt_bound) + rt_bound;
  hyper_vcoeff = exp(vcoeff);
  
  // local section for forward sampled quantities 
     { 
      for (s in 1:nS) {
        
        log_lik[s]=0;
        int nT_s = nT_subj[s];
        vector[2] tmp_wiener;
        real Trial_Drift;   // used for log-lik estimation
        real Sim_Drift;     // used for synthetic datasets
        
       for (t in 1:nT_s) {
         
          // estimated action values
          Q_selfish[s,t] = a_selfish_norm[s,t];                                                                                            
          Q_generous[s,t] = a_generous_norm[s,t] + (a_generous_norm[s,t] / (1 + k_effective[s] * social_distance[s,t]));
          
          // upper or lower bound choice
          if (choice[s,t]==1) {                                                                                 
            Trial_Drift = vcoeff_effective[s] * (Q_generous[s,t] - Q_selfish[s,t]);
          } else { 
            Trial_Drift = vcoeff_effective[s] * (Q_selfish[s,t] - Q_generous[s,t]);
          }
          
          // Simulated choice is independent from empirical choice
          Sim_Drift = vcoeff_effective[s] * (Q_generous[s,t] - Q_selfish[s,t]);
          
          // log-likelyhood and model forward sampled predicted choice
          log_lik[s] += wiener_lpdf(rt[s,t] | a_effective[s], ter_effective[s], 0.5, Trial_Drift);
          tmp_wiener = wiener_rng(a_effective[s], ter_effective[s], 0.5, Sim_Drift);
          rt_pred[s, t] = tmp_wiener[1];
          y_pred[s,t] = (tmp_wiener[2] == 1) ? 1 : 2;
          
       }
     }
   }
}
