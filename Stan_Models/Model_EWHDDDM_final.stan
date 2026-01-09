//============================================================================//
// Empathy-Weighted Hyperbolic-Discount Drift-Diffusion Model 
//============================================================================//

// Author: Fabrice Hubschmid
// Contact: fabrice.hubschmid@hhu.de

// Value function: Hyperbolic (k) + empathy weighting (f)
// Link function: Drift-Diffusion (a,ter,v(coeff)) linear mapping

// Disclaimers / Context:

//  - Resembles the model of Peters & D'Esposito (2020)
//    https://doi.org/10.1371/journal.pcbi.1007615
//  - with an additonal scaling parameter that weights the value of the reward
//    obtained by the other person.

//  - Code structured as the models of the hBayesDM package (Ahn & al. 2017)
//  - Parameter prior for the discount rate comes from Ahn & al. (2017)
//    https://doi.org/10.1162/cpsy_a_00002

//  - Parameter priors for the DDM come from Matzke & Wagenmakers (2009)
//    https://doi.org/10.3758/PBR.16.5.798

//  - Bias is fixed to 0.5 to account for trades-offs of explained
//    variance with the empathy weight (f)

//============================================================================//

// INFO: in order for the script to simulate based on wiener first passage time
// distributions, you need to insert the function that can be found online
// @: 
// https://discourse.mc-stan.org/t/alternative-to-wiener-rng-in-stan/28443/11

// Additionally activate lines indicated bellow

// ***insert function here***

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
  int <lower=0,upper=2> choice[nS,nT_max];    // stimulus choice (0 for missing 1 for generous and 2 for selfish)
  
  // Configuration of DDM boundaries
  //real <lower=0,upper=10> rt[nS, nT_max];   // reaction time in a trial in case of empirical data
  real <lower=0> rt[nS, nT_max];              // reaction time in a trial in case of simulated data
  real <lower=0> min_rt[nS];                  // minimum reaction time for each subject
  real <lower=0> rt_bound;                    // global minumum bound for non-decision time
}
//-----------//
// Model Parameters
//-----------//
// Raw parameters
parameters {
  // Valuation Model
  real k;                   // discount rate
  real <lower=0> k_sd;      // SD of discount rate
  real f;                   // empathy weight
  real <lower=0> f_sd;      // SD of the empathy weight
  
  // Response Model
  real a;                   // boundary separation
  real <lower=0> a_sd;      // SD of boundary separation
  real ter;                 // non-decision times
  real <lower=0> ter_sd;    // SD of non-decision times
  real v;                   // drift rate (scale factor)
  real <lower=0> v_sd;      // SD of drift rate
  
  // Individual level parameters
  vector [nS] k_s;          // individual raw discount rate
  vector [nS] f_s;          // individual raw empathy weight
  vector [nS] a_s;          // individual raw boundaries
  vector [nS] ter_s;        // individual raw non-decisions
  vector [nS] v_s;          // individual raw scalling drift rate
}
//-----------//
// Transformed parameters
//-----------//
// Effective parameters used in the sampling statement
transformed parameters {
  vector <lower=0, upper=1> [nS] k_effective;
  vector <lower=0> [nS] f_effective;                               
  vector <lower=0> [nS] a_effective;
  vector <lower=rt_bound, upper=max(min_rt)> [nS] ter_effective;
  vector <lower=0> [nS] v_effective;
  
  for (i in 1:nS) {
    k_effective[i] = Phi_approx(k + k_sd * k_s[i]);
    f_effective[i] = exp(f + f_sd * f_s[i]);
    a_effective[i] = exp(a + a_sd * a_s[i]);
    ter_effective[i] = Phi_approx(ter + ter_sd * ter_s[i]) * (min_rt[i] - rt_bound) + rt_bound;
    v_effective[i] = exp(v + v_sd * v_s[i]);
  }
}
//-----------//
// Model
//-----------//
model {
  // Model Priors for hyper parameters
  k ~ normal(0,1);
  f ~ normal(1,1);                          
  a ~ normal(1,1);
  ter ~ normal(0.435,0.12);    
  v ~ normal(2,3);                                      

  // Model Priors for SD parameters
  k_sd ~ normal(0,0.2);                   
  f_sd ~ normal(0,0.2);
  a_sd ~ normal(0,0.2);
  ter_sd ~ normal(0,0.2);
  v_sd ~ normal(0,0.2);
  
  // Model Priors for individual level parameters 
  k_s ~ normal(0,1);
  f_s ~ normal(0,1);
  a_s ~ normal(0,1);
  ter_s ~ normal(0,1);
  v_s ~ normal(0,1);
  
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
        Q_selfish = amount_selfish[s,t];                                                                                                     
        Q_generous = amount_generous[s,t] + (amount_generous[s,t] * f_effective[s]) / (1 + k_effective[s] * social_distance[s,t]);      
        
        // upper or lower bound choice
        if (choice[s,t]==1) {                                                                      
          Drift = v_effective[s] * (Q_generous - Q_selfish);
        } else { 
          Drift = v_effective[s] * (Q_selfish - Q_generous);
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
  real <lower=0,upper=1> hyper_k;
  real <lower=0> hyper_f;
  real <lower=0> hyper_a;
  real <lower=rt_bound, upper=max(min_rt)> hyper_ter;
  real <lower=0> hyper_v;
  real log_k;
  
  // Log likelihood
  real log_lik[nS];
  
  // Forward sampled quantities
  real Q_selfish[nS,nT_max];
  real Q_generous[nS,nT_max];
  //real y_pred[nS,nT_max];                         // ACTIVATE FOR SIMULATION
  //real rt_pred[nS, nT_max];                       // ACTIVATE FOR SIMULATION
  
  // Back transformed hyper parameters into interpretable space
  hyper_k = Phi_approx(k);
  hyper_f = exp(f);                                               
  log_k = log(hyper_k);
  hyper_a = exp(a);
  hyper_ter = Phi_approx(ter) * (mean(min_rt) - rt_bound) + rt_bound;
  hyper_v = exp(v);
  
  // local section for forward sampled quantities 
     { 
      for (s in 1:nS) {
        
        int nT_s = nT_subj[s];
        //vector[2] tmp_wiener;                         // ACTIVATE FOR SIMULATION
        log_lik[s]=0;
        
        real Trial_Drift;
        
       for (t in 1:nT_s) {
         
          // estimated action values
          Q_selfish[s,t] = amount_selfish[s,t];                                                                                                            
          Q_generous[s,t] = amount_generous[s,t] + (amount_generous[s,t] * f_effective[s]) / (1 + k_effective[s] * social_distance[s,t]);             
          
          // upper or lower bound choice
          if (choice[s,t]==1) {                                                                                 
            Trial_Drift = v_effective[s] * (Q_generous[s,t] - Q_selfish[s,t]);
          } else { 
            Trial_Drift = v_effective[s] * (Q_selfish[s,t] - Q_generous[s,t]);                                                       
          }
          
          // log-likelyhood and model forward sampled predicted choice
          log_lik[s] += wiener_lpdf(rt[s,t] | a_effective[s], ter_effective[s], 0.5, Trial_Drift);    
          //tmp_wiener = wiener_rng(a_effective[s], ter_effective[s], 0.5, Trial_Drift);                         // ACTIVATE FOR SIMULATION
          //rt_pred[s, t] = tmp_wiener[1];                                                                       // ACTIVATE FOR SIMULATION
          //y_pred[s, t] = to_int(tmp_wiener[2]);                                                                // ACTIVATE FOR SIMULATION
       } 
     }
   }
}


