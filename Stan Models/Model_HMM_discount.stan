//============================================================================//
// Hyperbolic discount model with dynamic noise estimation (HMM)
//============================================================================//

// Author: Fabrice Hubschmid
// Contact: hubschmidfabrice@gmail.com

// Value function: Hyperbolic (k)
// Link function: SoftMax (fixed inverse temperature (0.4))
// Dynamic noise estimation: P(lapse), P(recover)

// Disclaimers / Context:

//  - Traditional implementation of a hyperbolic social discounting model
//    f.e.: https://doi.org/10.1371/journal.pcbi.1007615
//  - Hidden Markov modelling of engagment states over choices sequences as 
//    suggested in Li & al. (2024)
//    https://doi.org/10.1016/j.jmp.2024.102842 
//  - Adapted some unpublished code very generously shared by Tim Angelike
//    for the mixed log target density estimation
//    https://scholar.google.de/citations?user=alXNO-gAAAAJ&hl=de&oi=ao

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
// Transformed Data                       
//-----------//
transformed data {
  vector[2] init_States = rep_vector(0.5,2);   // initial Probability of engagment [1] = engaged, [2] = random
}
//-----------//
// Model Parameters
//-----------//
// Raw parameters
parameters {
  // Individual level parameters
  vector <lower=0> [nS] k_s;              // individual discount rate
  vector <lower=0, upper=1> [nS] P_R2E_s; // individual probably of switching policy from random to engaged state (recover)
  vector <lower=0, upper=1> [nS] P_E2R_s; // individual probably of switching policy from engaged to random state (lapse)
}
//-----------//
// Model
//-----------//
model {
  
  // Model priors for individual level parameters 
  k_s     ~ normal(0,1);
  P_R2E_s ~ normal(0,1);
  P_E2R_s ~ normal(0,1);
  
  {
  // iterates over subjects
    for (s in 1:nS) {
      
      int nT_s = nT_subj[s];
      vector[2] Q_values;
      // state probabily over sequence
      vector[2] P_states;  // P[1] P being in engaged state - P[2] being in random state
      vector[2] theta;     // mixted target log density of actions
      real P_eng;          // used as denominator for state update
      
      // state initialization
      P_states = init_States;
      
      // iterates over trials
      for (t in 1:nT_s) {
        
        // action probabilities
        vector[2] A_pro;                      // [1] action probability for 1 - [2] action probability for 2
        vector[2] A_ran = rep_vector(0.5,2);  // action probabilities for random state
        
        // estimated choice values
        Q_values[1] = amount_generous[s,t] + (amount_generous[s,t] / (1 + k_s[s] * social_distance[s,t]));
        Q_values[2] = amount_selfish[s,t];
        
        // increment target log density for choice given both policies 
        A_pro = softmax(0.4 * Q_values);
        theta = A_pro * P_states[1] + A_ran * P_states[2];
        choice[s,t] ~ categorical(theta);
        
        // updating probabilities of being in engaged or random state
        P_eng = A_pro[choice[s,t]] * P_states[1] * (1-P_E2R_s[s]) + A_ran[choice[s,t]] * P_states[2] * P_R2E_s[s];
        P_states[1] = P_eng / theta[choice[s,t]];
        P_states[2] = 1 - P_states[1];
        
      }
    }
  }
}
//-----------//
// Model extracted quantities
//-----------//
generated quantities {
  // forward sampled quantities
  real log_lik[nS];
  real y_pred[nS,nT_max];
  real pred_P_S[nS, nT_max];
  real Q_gen[nS, nT_max];

  // local section for forward sampled quantities 
     { 
      for (s in 1:nS) {
        
        int nT_s = nT_subj[s];
        vector[2] Q_vals;
        vector[2] P_states;  // [1] Pr of being in engaged state - [2] Pr of being in random state
        vector[2] theta;     // mixted target log density of actions
        real P_eng;   
        
        // state initialization
        P_states = init_States;
        log_lik[s]=0;

       for (t in 1:nT_s) {
         
          // records predicted task engagement
          pred_P_S[s,t] = P_states[1];
         
          vector[2] A_pro;                      // [1] action probability for 1 - [2] action probability for 2
          vector[2] A_ran = rep_vector(0.5,2);  // action probabilities for random state
         
          // estimated action values
          Q_vals[1] = amount_generous[s,t] + (amount_generous[s,t] / (1 + k_s[s] * social_distance[s,t]));
          Q_vals[2] = amount_selfish[s,t];
          Q_gen[s,t] = Q_vals[1];
          
          // increment target log density for choice given both policies 
          A_pro = softmax(0.4 * Q_vals);
          theta = A_pro * P_states[1] + A_ran * P_states[2];
          log_lik[s] += categorical_lpmf(choice[s,t] | theta);
          y_pred[s,t] = categorical_rng(theta);
          
          // updating probabilities of being in engaged or random state
          P_eng = A_pro[choice[s,t]] * P_states[1] * (1-P_E2R_s[s]) + A_ran[choice[s,t]] * P_states[2] * P_R2E_s[s];
          P_states[1] = P_eng / theta[choice[s,t]];
          P_states[2] = 1 - P_states[1];
          
          // n.b: in order for the script to simulate data for recovery 
          // simulation you need to update the state probabilities over 
          // simulated sequences (y_pred[s,t]) and not choices (choice[s,t])
          
       } 
     }
   }
}
