# Rinit -------------------------------------------------------------------

buildRinit <- function(version = "gen") {
  
  seis_rinit_epi <- Csnippet(
    "
    double *S = &S1;
    double *E = &E1;
    double *Pone = &Pone1;
    double *Ptwo = &Ptwo1;
    double *Pthree = &Pthree1;
    double *I = &I1;
    double *Sprime = &Sprime1;
    double *C = &C1;
    double *X = &X1;
    double *Y = &Y1;
    const double *S_0 = &S_01;
    int u;
    for(u = 0; u < U; u++) {
      if (t_vacc[u] > 0 && t_vacc[u] < t) {
        S[u] = 0;
        Sprime[u] = nearbyint(S_0[u*S_0_unit]);
        E[u] = 1 - Sprime[u];
        C[u] = 1;
        X[u] = t_vacc_fake;
      } else {
        S[u] = nearbyint(S_0[u*S_0_unit]);
        Sprime[u] = 0;
        E[u] = 1 - S[u];
        C[u] = E[u];
        X[u] = 0;
      }
      
      Pone[u] = 0;
      Ptwo[u] = 0;
      Pthree[u] = 0;
      I[u] = 0;
      Y[u] = X[u];
    }
  "
  )
  
  seis_rinit_epi_v2 <- Csnippet(
    "
    double *S = &S1;
    double *E = &E1;
    double *Pone = &Pone1;
    double *Ptwo = &Ptwo1;
    double *Pthree = &Pthree1;
    double *I = &I1;
    double *Sprime = &Sprime1;
    double *C = &C1;
    double *X = &X1;
    double *Y = &Y1;
    const double *S_0 = &S_01;
    const double *C_0 = &C_01;

    int u;
    for(u = 0; u < U; u++) {
      int s0 = nearbyint(S_0[u*S_0_unit]);
      int c0 = nearbyint(C_0[u*C_0_unit]);
      
      if (t_vacc[u] > 0 && t_vacc[u] < t) {
        S[u] = 0;
        Sprime[u] = s0;
        C[u] = 1;
        X[u] = t_vacc_fake;
      } else {
        if (s0 == 1 && c0 == 1) {
          S[u] = 0;
          Sprime[u] = 1;
          C[u] = 1;
          X[u] = t_vacc_fake;
        } else if (s0 == 1 && c0 == 0) {
          S[u] = 1;
          Sprime[u] = 0;
          C[u] = 0;
          X[u] = 0;
        } else {
          S[u] = 0;
          Sprime[u] = 0;
          C[u] = 1;
          X[u] = t_vacc_fake;
        }
      }
      
      
      E[u] = 1 - s0;
      Pone[u] = 0;
      Ptwo[u] = 0;
      Pthree[u] = 0;
      I[u] = 0;
      Y[u] = X[u];
    }
  "
  )
  
  # Version allowing for Sprime to be initialized
  seis_rinit_gen_v2 <- Csnippet(
    "
    double *S = &S1;
    double *E = &E1;
    double *Pone = &Pone1;
    double *Ptwo = &Ptwo1;
    double *Pthree = &Pthree1;
    double *I = &I1;
    double *Sprime = &Sprime1;
    double *C = &C1;
    double *X = &X1;
    double *Y = &Y1;
    double *genlik = &genlik1;
    double *gencnt = &gencnt1;
    double *hazard = &hazard1;
    double *ratio = &ratio1;
    double *C_block = &C_block1;


    const double *S_0 = &S_01;
    const double *C_0 = &C_01;

    int u;
    for(u = 0; u < U; u++) {
      int s0 = nearbyint(S_0[u*S_0_unit]);
      int c0 = nearbyint(C_0[u*C_0_unit]);
      
      if (t_vacc[u] > 0 && t_vacc[u] < t) {
        S[u] = 0;
        Sprime[u] = s0;
        C[u] = 1;
        X[u] = t_vacc_fake;
      } else {
        if (s0 == 1 && c0 == 1) {
          S[u] = 0;
          Sprime[u] = 1;
          C[u] = 1;
          X[u] = t_vacc_fake;
        } else if (s0 == 1 && c0 == 0) {
          S[u] = 1;
          Sprime[u] = 0;
          C[u] = 0;
          X[u] = 0;
        } else {
          S[u] = 0;
          Sprime[u] = 0;
          C[u] = 1;
          X[u] = t_vacc_fake;
        }
        //Rprintf(\"Done Initialization, %f, u: %i\\n\",t,u);
      }
      
      E[u] = 1 - s0;
      Pone[u] = 0;
      Ptwo[u] = 0;
      Pthree[u] = 0;
      I[u] = 0;
      genlik[u] = -999;
      gencnt[u] = 0;
      Y[u] = X[u];
      hazard[u] = 0;
      ratio[u] = 0;
    }
    
    {
    double block_counts[B][2];    // temporary storage of block counts
  
    // Initialize
    for (int b=0; b<B; b++) {
      block_counts[b][0] = 0;
      block_counts[b][1] = 0;
    }
    
    // Compute
    for (u=0; u<U; u++) { 
      block_counts[map_participant_block[u]][map_participant_age_cat[u]] += C[u];
    }
    
    // Assign
    for (u=0; u<U; u++) {
      C_block[u] = block_counts[map_participant_block[u]][map_participant_age_cat[u]];
    }
  }
  "
  )
  
  if (version == "epi") {
    return(seis_rinit_epi_v2)
  } else if (version == "gen") {
    return(seis_rinit_gen_v2)
  }
} 

# RProcess ----------------------------------------------------------------
buildRProcess <- function(version,
                          ...) {
  
  seis_rprocess_gen <- Csnippet(
    "
  const double *lambda_wt=&lambda_wt1;
  const double *lambda_a=&lambda_a1;
  const double *lambda_d=&lambda_d1;
  const double *lambda_o=&lambda_o1;
  const double *lambdaA_wt=&lambdaA_wt1;
  const double *lambdaA_a=&lambdaA_a1;
  const double *lambdaA_d=&lambdaA_d1;
  const double *lambdaA_o=&lambdaA_o1;
  const double *lambdaB_wt=&lambdaB_wt1;
  const double *lambdaB_a=&lambdaB_a1;
  const double *lambdaB_d=&lambdaB_d1;
  const double *lambdaB_o=&lambdaB_o1;
  double *S = &S1;
  double *E = &E1;
  double *Pone = &Pone1;
  double *Ptwo = &Ptwo1;
  double *Pthree = &Pthree1;
  double *I = &I1;
  double *Sprime = &Sprime1;
  double *C = &C1;
  double *X = &X1;
  double *Y = &Y1;
  double *genlik = &genlik1;
  double *gencnt = &gencnt1;
  double *hazard = &hazard1;
  double *ratio = &ratio1;
  double *C_block = &C_block1;

  double foi, tot_inf_rates, vacc;
  double rate[6], trans[8];
  int u,v;
  double day = (t-floor(t))*365;
  int year;    // which year are we in
  int variant_ind = get_variant(t);
  double scov2_gamma = scov2_gamma_vec[variant_ind];
  double scov2_eta = scov2_eta_vec[variant_ind];
  double scov2_phi = scov2_phi_vec[variant_ind];
  double scov2_delta[4];
  double cnt_inf = 0;
  double lambda, lambdaA, lambdaB;
  double lambda_vec[3], lambdaA_vec[3], lambdaB_vec[3];
  double variant_probs[3];
  //const double *prob_alpha = &prob_alpha1;
  //const double *prob_delta = &prob_delta1;
  //const double *prob_omicron = &prob_omicron1;

 
  for (int i=0; i<4; i++) {
    scov2_delta[i] = scov2_delta_mat[i][variant_ind];
  //  Rprintf(\"scov2_delta[%i]: %f\\n\",i,scov2_delta_mat[i][variant_ind]);
  }
 
  //Rprintf(\"Time: %f, variant: %i\\n\",t,variant_ind);


  if (t < 0.6656393) {
    year = 0;
  } else {
    year = 1;
  }

  variant_probs[0] = prob_alpha;
  variant_probs[1] = prob_delta;
  variant_probs[2] = prob_omicron;
    
  for (u=0; u<U; u++) {
    double prev_genlik = genlik[u];    // store previous genlik
     // initialize variant probs
    //variant_probs[0] = prob_alpha[u];
    //variant_probs[1] = prob_delta[u];
    //variant_probs[2] = prob_omicron[u];
    
    // initialize vectors
    lambda_vec[0] = lambda_a[u*lambda_unit];
    lambda_vec[1] = lambda_d[u*lambda_unit];
    lambda_vec[2] = lambda_o[u*lambda_unit];
    lambdaA_vec[0] = lambdaA_a[u*lambdaA_unit];
    lambdaA_vec[1] = lambdaA_d[u*lambdaA_unit];
    lambdaA_vec[2] = lambdaA_o[u*lambdaA_unit];
    lambdaB_vec[0] = lambdaB_a[u*lambdaB_unit];
    lambdaB_vec[1] = lambdaB_d[u*lambdaB_unit];
    lambdaB_vec[2] = lambdaB_o[u*lambdaB_unit];
    
    // compute mean lambda
    lambda = mean_lambda(lambda_vec, variant_probs);
    //if (ISNA(lambda) || isnan(lambda)) {
    //  Rprintf(\"Time: %f, variant: %i, lambda: %f\\n\",t,variant_ind, lambda);
    //  }
    //if (u == 0 && variant_ind == 2) {
    //  Rprintf(\"Time: %f, variant: %i, lambda: %f\\n\",t,variant_ind, lambda);
    //  Rprintf(\" prob: [%f, %f, %f]\\n\", variant_probs[0], variant_probs[1], variant_probs[2]);
    //  Rprintf(\" lambda: [%f, %f, %f]\\n\", lambda_vec[0], lambda_vec[1], lambda_vec[2]);
    //}
    
    for (int i=0; i<8; i++) {
      trans[i] = 0;
    }
    
    // Define vaccination rate for participant
    if (t_vacc[u] > 0) {
      if (t >= (t_vacc[u] - 1/365)) {
        vacc = 1e3;
      } else {
        vacc = 0;
      }
    } else {
      vacc = 0;
    }
       //Rprintf(\"Done vacc, %f, u: %i\\n\",t,u);

    // FOI from community
    // By variant
    if (variant_ind == 0) {
      foi = lambda_wt[u*lambda_unit];
    } else {
      foi = lambda;
    }
    
    //if (u == 0){
    //       Rprintf(\"Done community, time %f, variant %i, lambda %f, u: %i\\n\",t, variant_ind, foi,u);
    //}

    
    // loop over adjacency
    if (start_list[u][year] != 999) {
      for (int i = start_list[u][year]; i <= end_list[u][year]; i++) {
        double P_tmp, I_tmp, mask;
        int l = adj_list[i][1];
        I_tmp = I[l];
        P_tmp = (Pone[l] + Ptwo[l] + Pthree[l]);
        
        
        if (I_tmp > 0 || P_tmp > 0) {
          int v = get_variant(X[l]);
          mask = get_mask(t, n_masks[l], start_mask[l], end_mask[l]);
          
          // Keep adjacency in households
          if (adj_type[i] == -1) {
            mask = 1;
            //Rprintf(\"HH, %i - %i, link %i year %i, start %i, end %i\\n\",u, l, i, year, start_list[u][year], end_list[u][year]);
          }
          
          if (adj_type[i] < 1) {
            // FOI from same group or hh
            if (v == 0) {
              foi += lambdaA_wt[u*lambdaA_unit] * mask * (alpha[v] * P_tmp + I_tmp);
            } else {
              foi += lambdaA_vec[v-1] * mask * (alpha[v] * P_tmp + I_tmp);
            }
          } else {
            // FOI from same school
            if (v == 0) {
              foi += lambdaB_wt[u*lambdaB_unit] * mask * (alpha[v] * P_tmp + I_tmp);
            } else {
              foi += lambdaB_vec[v-1] * mask * (alpha[v] * P_tmp + I_tmp);
            }
          }
        }
      }
    }
    
    // define rates
    rate[0] = foi;    // S->E: infection
    rate[1] = vacc;   // S-Sprime due to vaccination
    rate[2] = scov2_gamma;   // E->P: rate from exposed to pre-symptomatic infectious
    rate[3] = scov2_eta;  // P->I: rate from pre-symptomatic to symptomatic infectious
    rate[4] = scov2_phi;     // I->Sprime: rate from infectious to susceptible
    rate[5] = foi * scov2_delta[get_variant(X[u])];  //Sprime->E rate from susceptible with previous infection to exposed
     
     // save FOI
    if (C[u] > 0) {
      hazard[u] = foi * scov2_delta[get_variant(X[u])];
    } else {
      hazard[u] = foi;
    }
    
    //if (u == 0){
    //  Rprintf(\"Done rates,time %f, u: %i, %f %f %f %f %f %f\\n\",t,u, rate[0], rate[1], rate[2], rate[3], rate[4], rate[5]);
    //  Rprintf(\"Reinfection X %f var %i delta %f foi %f\\n\",X[u], get_variant(X[u]), scov2_delta[get_variant(X[u])], foi * scov2_delta[get_variant(X[u])]);
    //}
         
    // transitions between classes
    reulermultinom(2,S[u],&rate[0],dt,&trans[0]);
    reulermultinom(1,E[u],&rate[2],dt,&trans[2]);
    reulermultinom(1,Pone[u],&rate[3],dt,&trans[3]);
    reulermultinom(1,Ptwo[u],&rate[3],dt,&trans[4]);
    reulermultinom(1,Pthree[u],&rate[3],dt,&trans[5]);
    reulermultinom(1,I[u],&rate[4],dt,&trans[6]);
    reulermultinom(1,Sprime[u],&rate[5],dt,&trans[7]);
    
    //if (u == 0){
    //  Rprintf(\"States before: S[u] = %f\\n\", S[u]);
    //  Rprintf(\"States before: E[u] = %f\\n\", E[u]);
    //  Rprintf(\"States before: Pone[u] = %f\\n\", Pone[u]);
    //  Rprintf(\"States before: Ptwo[u] = %f\\n\", Ptwo[u]);
    //  Rprintf(\"States before: Pthree[u] = %f\\n\", Pthree[u]);
    //  Rprintf(\"States before: I[u] = %f\\n\", I[u]);
    //  Rprintf(\"States before: Sprime[u] = %f\\n\", Sprime[u]);


    //  for (int i=0;i<8;i++) {
    //    Rprintf(\"Transitions: trans[%i] = %f\\n\",i, trans[i]);
    //  }
    //}

    S[u] += -trans[0] - trans[1];
    E[u] += trans[0] + trans[7] - trans[2];
    Pone[u] += trans[2] - trans[3];
    Ptwo[u] += trans[3] - trans[4];
    Pthree[u] += trans[4] - trans[5];
    I[u] += trans[5] - trans[6];
    Sprime[u] += trans[1] + trans[6] - trans[7];
    C[u] += trans[0] + trans[1]; 
    cnt_inf += trans[7] + trans[0];
    
   // if (isnan(Pone[u])) {
  //    for (int i=0;i<8;i++) {
  //      Rprintf(\"Transitions %i time %f: trans[%i] = %f\\n\",u, t, i, trans[i]);
  //    }
  //  }
    
    
    //if (u == 0){
    //  Rprintf(\"States after: S[u] = %f\\n\", S[u]);
    //  Rprintf(\"States after: E[u] = %f\\n\", E[u]);
    //  Rprintf(\"States after: Pone[u] = %f\\n\", Pone[u]);
    //  Rprintf(\"States after: Ptwo[u] = %f\\n\", Ptwo[u]);
    //  Rprintf(\"States after: Pthree[u] = %f\\n\", Pthree[u]);
    //  Rprintf(\"States after: I[u] = %f\\n\", I[u]);
    //  Rprintf(\"States after: Sprime[u] = %f\\n\", Sprime[u]);
    //  
    //  Rprintf(\"Check Sprime,time %f, u: %i, trans1 %f, trans6 %f,  trans7 %f, Sprime %f\\n\",t,u, trans[1], trans[6], trans[7], Sprime[u]);
    //  Rprintf(\"Sum states,time %f, u: %i, sum %f\\n\",t,u, S[u] + E[u] + Pone[u] + Ptwo[u] + Pthree[u] + I[u] + Sprime[u]);
    //}

    // Time of earliest infection
    if ((trans[0] + trans[7]) > 0 && X[u] == 0) {
       Y[u] = t;
    } else if ((trans[1]) > 0 && X[u] == 0) {
       Y[u] = t_vacc_fake;
    }
    
    // Time of latest infection
    if ((trans[0] + trans[7]) > 0) {
       X[u] = t;
    } else if ((trans[1]) > 0) {
       X[u] = t_vacc_fake;
    }
    
    
    // compute ratio
    ratio[u] = (foi - lambda)/lambda;
  
    
    //Rprintf(\"Done update, %f, u: %i\\n\",t,u);

    //masked[u] = get_mask(t, n_masks[u], start_mask[u], end_mask[u]);
    
    // Trick pomp into computing the genetic likelihood here to be able to use
    // the block particle filter. We compute the likelihood only at the time of
    // transition from the S or Sprime classes to the exposed class
    
    // Rest genlink 
    //genlik[u] = 0;
    
    if (use_gen_lik == 1) {
    double new_genlik = 0;
    
      if (time_seq[u] == -999) {
        // Case in which no sequence data is available
        genlik[u] = 0;
      // Attempt to alaways compute the likelihood
      } else if (t >= (time_seq[u] - seq_delay) && t < time_seq[u]  && gencnt[u] == 0.0) {
        // If the sequence is positive
        foi = 0;                                  // foi experience by particiapnt
        tot_inf_rates = 0;                        // total infection rate at time of infection
        int i = seq_covar[u];          // what sequence is it
        int ancestors[n_max_ancestors[i]];        // vector of ancestor sequences
        int infectors[n_max_ancestors[i]];        // vector of participants corresponding to each ancestor sequence. This points to element in adj_list to reuse relation type.
        int adj_type_tmp[n_max_ancestors[i]];     // vector of participants corresponding to each ancestor sequence. This points to element in adj_list to reuse relation type.
        double rel_rates[n_max_ancestors[i]+1];   // relative rates for prior probability
        double marglik_gen_vec[max_kappa];        // relative rates for prior probability
        int j = 0;                                // counter for number of ancestors
  
        // Update status of counted gen likelihood infection occurred
        if ((trans[0] + trans[7]) == 1) {
         gencnt[u] = 1;
        }
        
        // If no ancestors, then the genetic likelihood is only composed of distance to
        // community sequences.
  
        if (n_max_ancestors[i] == 0 || start_list[u][year] == 999) {
          for (int kappa = 0; kappa < max_kappa; kappa++) {
            // Don't update if we consider that mu is known, if not could uncomment line
            marglik_gen_vec[kappa] = dbinom(dist_closest_comm_seq[i], closest_comm_seq_len[i], (kappa+1)*mu, 1) +
              log(w_mat[kappa][closest_comm_seq_temp_diff[i]]);
          }
  
          new_genlik = log_sum_exp(marglik_gen_vec, max_kappa);
  
        } else {
  
          // Initialize relative rates
          for (int k = 0; k<(n_max_ancestors[i]+1); k++) {
            rel_rates[k] = 0;
          }
  
          // First accumulate infected participants
          //Rprintf(\"Part: %i; seq: %i; start %i, end %i \\n\",u, i, start_list_seq[i], end_list_seq[i]);

          for (int k = start_list_seq[i]; k <= end_list_seq[i]; k++) {
            int l = adj_list_infectors[k];
            //Rprintf(\"\t k: %i, l: %i \\n\",k, l);
            if (C[l] > 0 && X[l] > t_vacc_fake) {
              // Add to vector of ancestors if ever infected
              ancestors[j] = adj_list_ancestors[k];
              infectors[j] = l;
              adj_type_tmp[j] = adj_type_infectors[k];
              j++;
            }
          }
  
          // loop over adjacency
          for (int k = start_list[u][year]; k <= end_list[u][year]; k++) {
            double P_tmp, I_tmp, mask;
            int l = adj_list[k][1];
            I_tmp = I[l];
            P_tmp = (Pone[l] + Ptwo[l] + Pthree[l]);
            
            if (I_tmp >0 || P_tmp >0) {
              int v = get_variant(X[l]);
              mask = get_mask(t, n_masks[l], start_mask[l], end_mask[l]);
              
              // Keep adjacency in households
              if (adj_type[k] == -1) {
                mask = 1;
              }
              
              if (adj_type[k] < 1) {
                // FOI from same group or hh
                if (v == 0) {
                  foi += lambdaA_wt[u*lambdaA_unit] * mask * (alpha[v] * P_tmp + I_tmp);
                } else {
                  foi += lambdaA_vec[v-1] * mask * (alpha[v] * P_tmp + I_tmp);
                }
               } else {
                // FOI from same school
                if (v == 0) {
                  foi += lambdaB_wt[u*lambdaB_unit] * mask * (alpha[v] * P_tmp + I_tmp);
                } else {
                  foi += lambdaB_vec[v-1] * mask * (alpha[v] * P_tmp + I_tmp);
                }
              }
            }
          }
  
          // Compute relative importance of infection sources
          if (j>0) {
            double n_non_inf = 0;    // number of participants with history of infection but non infected
  
            for (int k = 0; k < j; k++) {
              double P_tmp, I_tmp, mask;
              int l = infectors[k];
              I_tmp = I[l];
              P_tmp = (Pone[l] + Ptwo[l] + Pthree[l]);
  
              if ((I_tmp + P_tmp) == 1.0) {
                int v = get_variant(X[l]);
                mask = get_mask(t, n_masks[l], start_mask[l], end_mask[l]);
                
                // Keep adjacency in households
                if (adj_type_tmp[k] == -1) {
                  mask = 1;
                }
  
                if (adj_type_tmp[k] < 1) {
                  // FOI from same group
                  if (v == 0) {
                    rel_rates[k] = lambdaA_wt[u*lambdaA_unit] * mask * (alpha[v] * P_tmp + I_tmp);
                  } else if (variant_ind == 1) {
                    rel_rates[k] = lambdaA_vec[v-1] * mask * (alpha[v] * P_tmp + I_tmp);
                  }
                } else {
                  // FOI from same school
                  if (v == 0) {
                    rel_rates[k] = lambdaB_wt[u*lambdaB_unit] * mask * (alpha[v] * P_tmp + I_tmp);
                  } else if (variant_ind == 1) {
                    rel_rates[k] = lambdaB_vec[v] * mask * (alpha[v] * P_tmp + I_tmp);
                  }
               }
              } else {
                n_non_inf++;
              }
            }
  
            // Set total infection hazard coming from participants with sequences
            for (int k = 0; k < j; k++) {
              tot_inf_rates += rel_rates[k];
            }
  
            // now set relative rates for participants not currently infected
            for (int k = 0; k < j; k++) {
              int m = infectors[k];
              double mask = get_mask(t, n_masks[m], start_mask[m], end_mask[m]);
              
              // Keep adjacency in households
              if (adj_type_tmp[k] == -1) {
                mask = 1;
              }
              
              if ((mask * (Pone[m] + Ptwo[m] + Pthree[m] + I[m])) == 0.0) {
                // this assumes equal prior probability
                rel_rates[k] = (foi-tot_inf_rates)/n_non_inf;
              }
            }
          }
  
          // Add community transmission to vector
          if (variant_ind == 0) {
            rel_rates[j] = lambda_wt[u*lambda_unit];
          } else {
            rel_rates[j] = lambda;
          }
          j++;
  
          // compute likelihoods
          if (foi == 0) {
            // if none infected then must come from community
            for (int kappa = 0; kappa < max_kappa; kappa++) {
              marglik_gen_vec[kappa] = dbinom(dist_closest_comm_seq[i], closest_comm_seq_len[i], (kappa+1)*mu, 1) +
              log(w_mat[kappa][closest_comm_seq_temp_diff[i]]);
            }
  
            new_genlik = log_sum_exp(marglik_gen_vec, max_kappa);
            //Rprintf(\"No FOI: infection genlik %i, %f \\n\",u, new_genlik);
  
          } else {
            // Total foce of infection
            double tot_foi = 0;
  
            for (int k=0;k<j;k++) {
              tot_foi += rel_rates[k];
            }
  
            // Normalize force of infection
            for (int k=0;k<j;k++) {
              rel_rates[k] = rel_rates[k]/tot_foi;
            }
  
            {
              double marglik_vec[j];
              // now accumulate likelihood
              for (int k=0;k<(j-1);k++) {
                int m = infectors[k];
                double mask = get_mask(t, n_masks[m], start_mask[m], end_mask[m]);
              
                // Keep adjacency in households
                if (adj_type_tmp[k] == -1) {
                  mask = 1;
                }
  
                // If infector currently not infectious or not connected, marginalize out number of generations
                if (mask * (Pone[m] + Ptwo[m] + Pthree[m] + I[m]) == 0.0) {
                  // Likelihood for one gen is 0
                  marglik_gen_vec[0] = -1e8;
                  for (int kappa = 1; kappa < max_kappa; kappa++) {
                    marglik_gen_vec[kappa] =
                      dbinom(dist_seq[i][ancestors[k]], seq_len_mat[i][ancestors[k]], (kappa+1)*mu, 1) +
                      log(w_mat[kappa][temp_diff_seq[i][ancestors[k]]]);
                      //Rprintf(\"    Maringalizing no infectors seq %i kappa %i ll %f delta %i l %i mu %f prob %f dt %f \\n\", i, kappa, dbinom(dist_seq[i][ancestors[k]], seq_len_mat[i][ancestors[k]], (kappa+1)*mu, 1), dist_seq[i][ancestors[k]],  seq_len_mat[i][ancestors[k]], mu, log(w_mat[kappa][temp_diff_seq[i][ancestors[k]]]), temp_diff_seq[i][ancestors[k]]);
                  }
  
                  marglik_vec[k] = log(rel_rates[k]) + log_sum_exp(marglik_gen_vec, max_kappa);
                  //Rprintf(\"    Element no infectors %i  ll %f prob %f \\n\", k, log_sum_exp(marglik_gen_vec, max_kappa),  log(rel_rates[k]));
                } else {
                  marglik_vec[k] = log(rel_rates[k]) + dbinom(dist_seq[i][ancestors[k]], seq_len_mat[i][ancestors[k]], mu, 1);
                  //Rprintf(\"    Element with infectors %i  ll %f prob %f \\n\", k, dbinom(dist_seq[i][ancestors[k]], seq_len_mat[i][ancestors[k]], mu, 1),  log(rel_rates[k]));
                }
              }
  
              for (int kappa = 0; kappa < max_kappa; kappa++) {
                marglik_gen_vec[kappa] =  dbinom(dist_closest_comm_seq[i], closest_comm_seq_len[i], (kappa+1) * mu, 1) +
                  log(w_mat[kappa][closest_comm_seq_temp_diff[i]]);
                  
                  //Rprintf(\"    Maringalizing seq %i kappa %i ll %f prob %f \\n\", i, kappa, dbinom(dist_closest_comm_seq[i], closest_comm_seq_len[i], (kappa+1) * mu, 1), log(w_mat[kappa][closest_comm_seq_temp_diff[i]]));

              }
  
              // add community by hand
              marglik_vec[j-1] =  log(rel_rates[j-1]) + log_sum_exp(marglik_gen_vec, max_kappa);
              new_genlik = log_sum_exp(marglik_vec, j);
              //Rprintf(\"WITH FOI: infection genlik %i, %f from community %f \\n\",u, new_genlik, log_sum_exp(marglik_gen_vec, max_kappa));
            }
          }
        }
        
        if ((trans[0] + trans[7]) == 1) {
           genlik[u] = new_genlik;
            //Rprintf(\"Infection genlik %i, %f \\n\",u, new_genlik);
        } else {
          if (genlik[u] > new_genlik || genlik[u] == -999) {
            genlik[u] = new_genlik;
            //Rprintf(\"Updating genlik %i, %f \\n\",u, new_genlik);
          }
        }
      } else if (t >= time_seq[u] && gencnt[u] == 0.0)  {
        // Case in which sequence data is available but model does not predict infection
        //genlik[u] = -10;
        gencnt[u] = 1;
        //Rprintf(\"Endline genlik %i, %f \\n\",u, genlik[u]);
      }
    }
  }
  
  // Compute number of infected per block 
  {
    double block_counts[B][2];    // temporary storage of block counts
  
    // Initialize
    for (int b=0; b<B; b++) {
      block_counts[b][0] = 0;
      block_counts[b][1] = 0;
    }
    
    // Compute
    for (u=0; u<U; u++) { 
      block_counts[map_participant_block[u]][map_participant_age_cat[u]] += C[u];
    }
    
    //if (t > (t_endline - 2/365)) {
    //  for (int i=0;i<B;i++) {
    //   for (int j=0;j<2;j++) {
     //   Rprintf(\"Block counts %f b %i a %i: %f\\n\",t,i, j,  block_counts[i][j]);
     //  }
      //}
    //}
    
    // Assign
    for (u=0; u<U; u++) {
      C_block[u] = block_counts[map_participant_block[u]][map_participant_age_cat[u]];
      // if (t > (t_endline - 2/365)) {
       // Rprintf(\"Block counts for %u: %f\\n\",u, C_block[u]);
     // }
    }
  }
    
     //Rprintf(\"Total infected in time %f: %f\\n\",t,cnt_inf);
    ")
  
  if (version == "epi") {
    return(seis_rprocess_epi)
  } else if (version == "gen") {
    return(seis_rprocess_gen)
  }
}

# DMeasure ----------------------------------------------------------------

buildDUnitMeasure <- function(version = NULL) {
  # This is the code for observations.
  # Start with pcr test.
  seis_dunit_measure_gen <- Csnippet(
    "
    double this_time = t;
    double epsilon = 1e-3;
    
    lik = 0;
    if (!ISNA(n_pcr) && use_test_lik == 1) {
      lik += dbinom(nearbyint(pcr_pos), nearbyint(n_pcr), (Pone+Ptwo+Pthree+I) * sens_pcr + (1-(Pone+Ptwo+Pthree+I)) * (1-spec_pcr), 1);
      if (ISNA(dbinom(nearbyint(pcr_pos), nearbyint(n_pcr), (Pone+Ptwo+Pthree+I) * sens_pcr + (1-(Pone+Ptwo+Pthree+I)) * (1-spec_pcr), 1))) {
        Rprintf(\"PCR lik is NA at time %f: %f %f %f %f pcr:%f pcrpos:%f\\n\", t, Pone, Ptwo, Pthree, I, pcr_pos,  n_pcr);
      }
    } 
                 //Rprintf(\"Done pcr\\n\");

    if (!ISNA(sero) && use_sero_lik == 1) {
      double t_log[4];
      double ll[n_samples];    // log-lik to marginalize
      double t_inf;
      
      if (Y == 0) {
        t_inf = 60;
      } else {
        t_inf = (this_time - Y)*365;
      }
      
      t_log[0] = 1;
      t_log[1] = (log(t_inf) - t_log_mean)/t_log_sd;
      t_log[2] = pow(t_log[1], 2);
      t_log[3] = pow(t_log[1], 3);
      
      for (int i=0; i<n_samples; i++) {
        double this_sens;
        this_sens = inv_logit(dot_product_2(t_log, time_sens_coefs[i], 4));
        ll[i] = log_prior + dbinom(nearbyint(sero), 1, C * this_sens + (1-C) * (1-spec_sero), 1);
      }
      
      lik += log_sum_exp(ll, n_samples);
      if (ISNA(log_sum_exp(ll, n_samples))) {
        //Rprintf(\"SERO lik is NA\\n\");
      }
            //Rprintf(\"Done sero\\n\");
      //lik += dbinom(nearbyint(sero), 1, C * sens_sero + (1-C) * (1-spec_sero), 1);
    }
    
    if (!ISNA(seq) && use_gen_lik == 1) {
      if (ISNA(genlik) | isnan(genlik)) {
        lik += -999;
      } else {
        lik += genlik;
      }
      if (ISNA(genlik)) {
        //Rprintf(\"genlik lik is NA\\n\");
      }
    }
    
    if (!ISNA(age_cat)) {
      
      int a = nearbyint(age_cat);
      double lp[n_draws];
      


      if (t > (t_baseline - epsilon) && t < (t_baseline + epsilon)) {
        for (int i=0;i<n_draws;i++) {
          lp[i] = -log(n_draws) + dbinom(nearbyint(C_block), nearbyint(n_block), seroprev_baseline[i][a], 1);
        }
        //Rprintf(\"Done baseline %f\\n\",  log_sum_exp(lp, n_draws));
      } else if (t > (t_middle - epsilon) && t < (t_middle + epsilon)) {
        for (int i=0;i<n_draws;i++) {
          lp[i] = -log(n_draws) + dbinom(nearbyint(C_block), nearbyint(n_block), seroprev_middle[i][a], 1);
        }
       // Rprintf(\"Done middle %f\\n\",  log_sum_exp(lp, n_draws));
    } else if (t > (t_endline - epsilon) && t < (t_endline + epsilon)) {
      for (int i=0;i<n_draws;i++) {
        lp[i] = dbinom(nearbyint(C_block), nearbyint(n_block), seroprev_endline[i][a], 1);
          
        // Endline seroprevalence from SP4
        if (C_block < n_block) {
          double lp2[2];
          lp2[0] = lp[i];
          lp2[1] = pbinom(nearbyint(C_block), nearbyint(n_block), seroprev_endline[i][a], 0, 1);
          lp[i] = log_sum_exp(lp2, 2);
        }
        
        // Add prior
        lp[i] += -log(n_draws);
        
         //if (ISNA(lp[i]) || isnan(lp[i]) || isinf(lp[i])) {
          // Rprintf(\"i : %i, ll: %f\\n\", i, lp[i]);
         //}
      }
      
      // Rprintf(\"Done endline %f \\n\", log_sum_exp(lp, n_draws));
    }
    
   // Rprintf(\"Done age %f\\n\",  (1/n_block) * log_sum_exp(lp, n_draws));
    lik += (1/n_block) * log_sum_exp(lp, n_draws);
  }
  ")
  
  seis_dunit_measure_epi <- Csnippet(
    "
    lik = 0;
    if (!ISNA(n_pcr) && use_test_lik == 1) {
      lik += dbinom(nearbyint(pcr_pos), nearbyint(n_pos), (Pone+Ptwo+Pthree+I) * sens_pcr + (1-(Pone+Ptwo+Pthree+I)) * (1-spec_pcr), 1);
    } 
    
    if (!ISNA(sero) && use_sero_lik == 1) {
      lik += dbinom(nearbyint(sero), 1, C * sens_sero + (1-C) * (1-spec_sero), 1);
    }
  ")
  
  
  if (version == "epi") {
    return(seis_dunit_measure_epi)
  } else if (version == "gen") {
    return(seis_dunit_measure_gen)
  }
  
}

buildDMeasure <- function(version = NULL) {
  
  seis_dmeasure_epi <- Csnippet(
    "
  const double *Pone = &Pone1;
  const double *Ptwo = &Ptwo1;
  const double *Pthree = &Pthree1;
  const double *I = &I1;
  const double *C = &C1;
  const double *n_pcr = &n_pcr1;
  const double *pcr_pos = &pcr_pos1;
  const double *sero = &sero1;

  int u;
  lik = 0;
  
  for (u = 0; u<U; u++) {
  
  if (!ISNA(n_pcr[u])) {
  lik += dbinom(nearbyint(pcr_pos[u]), nearbyint(n_pcr[u]), (Pone[u]+Ptwo[u]+Pthree[u]+I[u]) * sens_pcr + (1-(Pone[u]+Ptwo[u]+Pthree[u]+I[u])) * (1-spec_pcr), 1);
  }
    if (!ISNA(sero[u])) {
      lik += dbinom(nearbyint(sero[u]), 1, C[u] * sens_sero + (1-C[u]) * (1-spec_sero), 1);
    }
  }
  
  ")
  
  
  if (version == "epi") {
    return(seis_dmeasure_epi)
  } else if (version == "gen") {
    return(NULL)
  }
  
}

# RMeasure ----------------------------------------------------------------

buildRMeasure <- function(version = NULL) {
  
  # Add serology
  seis_rmeasure <- Csnippet(
    "
  double *n_pcr = &n_pcr1;
  double *pcr_pos = &pcr_pos1;
  double *sero = &sero1;
  const double *I = &I1;
  const double *Pone = &Pone1;
  const double *Ptwo = &Ptwo1;
  const double *Pthree = &Pthree1;
  const double *C = &C1;
  int Isim, Csim;
  int u;
  

  for (u = 0; u < U; u++) {
         //Rprintf(\"Starting sim, u: %i\\n\",u);

    Isim = nearbyint(I[u]) +  nearbyint(Pone[u]) +  nearbyint(Ptwo[u]) +  nearbyint(Pthree[u]);
    Csim = nearbyint(C[u]);
    if (Isim == 1) {
      pcr_pos[u] = rbinom(1, sens_pcr);
    } else {
      pcr_pos[u] = rbinom(1, (1-spec_pcr));
    }
    if (Csim >= 1) {
      sero[u] = rbinom(1, sens_sero);
    } else {
      sero[u] = rbinom(1, (1-spec_sero));
    }
      //Rprintf(\"Done sim u: %i\\n\",u);
  }

")
  
  return(seis_rmeasure)
}

buildRUnitMeasure <- function(version = NULL) {
  seis_runit_measure <- spatPomp_Csnippet(
    "
  double pcr_pos;
  int Isim;
  Isim = nearbyint(I);
   if (Isim == 1) {
    pcr_pos = rbinom(1, sens_pcr);
  } else {
    pcr_pos = rbinom(1, (1-spec_pcr));
  }
")
  
  return(seis_runit_measure)
}


# Globals -----------------------------------------------------------------

replaceEmpty <- function(x) {
  if (is.null(dim(x))) {
    999
  } else {
    array(999, dim = rep(1, length(dim(x))))
  } 
}

mapReplaceEmpty <- function(x) {
  rapply(x, function(y) {
    if (purrr::is_empty(y)) {
      replaceEmpty(y)
    } else {
      y
    }
  }, 
  how = "replace")
}

checkGlobalTestPerf <- function(test_perf) {
  test_perf %>% 
    mapReplaceEmpty()
}

checkGlobalEpiAdj <- function(adj_epi) {
  adj_epi %>% 
    mapReplaceEmpty()
}

checkGlobalGen <- function(gen) {
  gen %>% 
    mapReplaceEmpty()
}

checkGlobalPars <- function(pars) {
  pars %>% 
    mapReplaceEmpty()
}

checkGlobalLiks <- function(liks) {
  liks %>% 
    mapReplaceEmpty()
}

checkGlobalNatHist <- function(nat_hist) {
  nat_hist %>% 
    mapReplaceEmpty()
}

checkGlobalMasksInfo <- function(masks_info) {
  masks_info %>% 
    mapReplaceEmpty()
}

checkGlobalsVaccInfo <- function(vaccination) {
  vaccination  %>% 
    mapReplaceEmpty() 
}

checkGlobalsSeroSens <- function(sero_sens) {
  sero_sens %>% 
    mapReplaceEmpty()
}

checkGlobalsSeroprev <- function(seroprev) {
  seroprev %>% 
    mapReplaceEmpty()
}

buildGlobalsWrapper <- function(parlist) {
  
  test_perf <- checkGlobalTestPerf(parlist$test_perf)
  adj_epi <- checkGlobalEpiAdj(parlist$adj_epi)
  gen <- checkGlobalGen(parlist$gen)
  pars <- checkGlobalPars(parlist$pars)
  liks <- checkGlobalLiks(parlist$liks)
  nat_hits <- checkGlobalNatHist(parlist$nat_hist)
  masks_info <- checkGlobalMasksInfo(parlist$masks_info)
  vaccination <- checkGlobalsVaccInfo(parlist$vaccination)
  sero_sens <- checkGlobalsSeroSens(parlist$sero_sens)
  seroprev <- checkGlobalsSeroprev(parlist$seroprev)
  
  buildGlobals(
    # Test peformance
    spec_pcr = test_perf$spec_pcr,
    sens_pcr = test_perf$sens_pcr,
    spec_sero = test_perf$spec_sero,
    sens_sero = test_perf$sens_sero,
    # Adjacency for epi
    adj_list = adj_epi$adj_list, 
    adj_type = adj_epi$adj_type, 
    start_list = adj_epi$start_list, 
    end_list = adj_epi$end_list, 
    # Genetic part
    do_gen = gen$do_gen,
    adj_list_ancestors = gen$adj_list_ancestors, 
    adj_list_infectors = gen$adj_list_infectors, 
    adj_type_infectors = gen$adj_type_infectors, 
    dist_closest_comm_seq = gen$dist_closest_comm_seq, 
    closest_comm_seq_len = gen$closest_comm_seq_len, 
    seq_len_mat = gen$seq_len_mat, 
    dist_seq = gen$dist_seq, 
    start_list_seq = gen$start_list_seq, 
    end_list_seq = gen$end_list_seq, 
    n_max_ancestors = gen$n_max_ancestors, 
    temp_diff_seq = gen$temp_diff_seq,
    closest_comm_seq_temp_diff = gen$closest_comm_seq_temp_diff,
    seq_covar = gen$seq_covar,
    w_mat = gen$w_mat,
    mu = gen$mu,
    max_kappa = gen$max_kappa,
    time_seq = gen$time_seq,
    seq_delay = gen$seq_delay,
    # Parameter names for rprocess
    set_expanded = pars$set_expanded,
    set_fixed = pars$set_fixed,
    # Likelihood flags
    use_test_lik = liks$use_test_lik,
    use_sero_lik = liks$use_sero_lik,
    use_gen_lik = liks$use_gen_lik,
    # Natural history parameters
    scov2_gamma = nat_hits$scov2_gamma,
    scov2_eta = nat_hits$scov2_eta,
    scov2_phi = nat_hits$scov2_phi,
    scov2_delta = nat_hits$scov2_delta,
    alpha = nat_hits$alpha,
    # Masks for vacations and quarantine
    n_masks = masks_info$n_masks,
    start_mask = masks_info$start_mask,
    end_mask = masks_info$end_mask,
    t_vacc = vaccination$t_vacc,
    t_vacc_fake = vaccination$t_vacc_fake,
    # Sensitivity of serology
    t_log_mean = sero_sens$t_log_mean,
    t_log_sd = sero_sens$t_log_sd,
    time_sens_coefs = sero_sens$time_sens_coefs,
    # Priors on seroprev at baseline
    seroprev_baseline = seroprev$seroprev_baseline,
    seroprev_middle = seroprev$seroprev_middle,
    seroprev_endline = seroprev$seroprev_endline,
    map_participant_block = seroprev$map_participant_block,
    map_participant_age_cat = seroprev$map_participant_age_cat
  )
}

makeNullGenGlobals <- function() {
  list(
    do_gen = F,
    adj_list_ancestors = NULL, 
    adj_list_infectors = NULL, 
    adj_type_infectors = NULL, 
    dist_closest_comm_seq = NULL, 
    closest_comm_seq_len = NULL, 
    seq_len_mat = NULL, 
    dist_seq = NULL, 
    start_list_seq = NULL, 
    end_list_seq = NULL, 
    n_max_ancestors = NULL, 
    temp_diff_seq = NULL,
    closest_comm_seq_temp_diff = NULL,
    w_mat = NULL,
    mu = NULL,
    max_kappa = NULL
  )
}

buildGlobals <- function(
    spec_pcr,
    sens_pcr,
    spec_sero,
    sens_sero,
    adj_list, 
    adj_type, 
    start_list, 
    end_list, 
    do_gen = F,
    adj_list_ancestors = NULL, 
    adj_list_infectors = NULL, 
    adj_type_infectors = NULL, 
    dist_closest_comm_seq = NULL, 
    closest_comm_seq_len = NULL, 
    seq_len_mat = NULL, 
    dist_seq = NULL, 
    start_list_seq = NULL, 
    end_list_seq = NULL, 
    n_max_ancestors = NULL, 
    temp_diff_seq = NULL,
    closest_comm_seq_temp_diff = NULL,
    seq_covar = NULL,
    w_mat = NULL,
    mu = NULL,
    max_kappa = NULL,
    set_expanded,
    set_fixed,
    use_test_lik,
    use_sero_lik,
    use_gen_lik,
    scov2_gamma,
    scov2_eta,
    scov2_phi,
    scov2_delta,
    alpha,
    n_masks,
    start_mask,
    end_mask,
    t_vacc_fake,
    t_vacc,
    time_seq,
    seq_delay,
    t_log_mean,
    t_log_sd,
    time_sens_coefs,
    age_cat,
    seroprev_baseline,
    seroprev_middle,
    seroprev_endline,
    map_participant_block,
    map_participant_age_cat
) {
  
  # A. Test performance ---
  test_perf <- Csnippet(
    str_glue(
      "
      double spec_pcr = {spec_pcr};
      double sens_pcr = {sens_pcr};
      double spec_sero = {spec_sero};
      double sens_sero = {sens_sero};
    "
    )
  )
  
  # Serology sensitivity
  sens_sero_dat <- str_c(
    matrixToCSnippet(time_sens_coefs, 
                     mat_name = "time_sens_coefs",
                     type = "double"),
    str_glue(
      "
      double t_log_mean = {t_log_mean};
      double t_log_sd = {t_log_sd};
      int n_samples = {nrow(time_sens_coefs)};
      double log_prior = {log(1/nrow(time_sens_coefs))};
      "
    ),
    sep = "\n"
  )
  
  # https://stackoverflow.com/questions/20733590/dot-product-function-in-c-language
  dot_product_C <- "
  double dot_product_2(double v[], const double u[], int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++)
        result += v[i]*u[i];
    return result;
    }
  "
  
  inv_logit_C <- "
  double inv_logit(double x) {
    double result = 1/(1+exp(-x));
    return result;
  }
  "
  
  # B. Vectors and matrices for epi part ---
  vec_mat <- str_c(
    matrixToCSnippet(adj_list, 
                     mat_name = "adj_list",
                     type = "int"),
    vectorToCSnippet(adj_type, 
                     vec_name = "adj_type",
                     type = "int"),
    matrixToCSnippet(start_list, 
                     mat_name = "start_list",
                     type = "int"),
    matrixToCSnippet(end_list, 
                     mat_name = "end_list",
                     type = "int"),
    sep = "\n"
  )
  
  
  # C. Vectors and matrices for gen part ---
  
  if (do_gen) {
    vec_mat_gen <- str_c(
      vectorToCSnippet(adj_list_ancestors, 
                       vec_name = "adj_list_ancestors",
                       type = "int"),
      vectorToCSnippet(adj_list_infectors, 
                       vec_name = "adj_list_infectors",
                       type = "int"),
      vectorToCSnippet(adj_type_infectors, 
                       vec_name = "adj_type_infectors",
                       type = "int"),
      vectorToCSnippet(dist_closest_comm_seq, 
                       vec_name = "dist_closest_comm_seq",
                       type = "int"),
      vectorToCSnippet(closest_comm_seq_len, 
                       vec_name = "closest_comm_seq_len",
                       type = "int"),
      matrixToCSnippet(seq_len_mat, 
                       mat_name = "seq_len_mat",
                       type = "int"),
      matrixToCSnippet(dist_seq, 
                       mat_name = "dist_seq",
                       type = "int"),
      vectorToCSnippet(start_list_seq, 
                       vec_name = "start_list_seq",
                       type = "int"),
      vectorToCSnippet(end_list_seq, 
                       vec_name = "end_list_seq",
                       type = "int"),
      vectorToCSnippet(n_max_ancestors, 
                       vec_name = "n_max_ancestors",
                       type = "int"),
      matrixToCSnippet(temp_diff_seq,
                       mat_name = "temp_diff_seq",
                       type = "int"),
      vectorToCSnippet(closest_comm_seq_temp_diff,
                       vec_name  = "closest_comm_seq_temp_diff",
                       type = "int"),
      matrixToCSnippet(w_mat,
                       mat_name = "w_mat",
                       type = "double"),
      vectorToCSnippet(time_seq,
                       vec_name = "time_seq",
                       type = "double"),
      vectorToCSnippet(seq_covar,
                       vec_name  = "seq_covar",
                       type = "int"),
      sep = "\n"
    )
    
    # D. Genetic parameters
    gen_param <- Csnippet(
      str_glue("double mu = {mu};
               int max_kappa = {max_kappa};
               double seq_delay = {seq_delay};")
    )
    
  } else {
    vec_mat_gen <- ""
    gen_param <- ""
  }
  
  # Function to compute average lambda
  mean_lambda_C <- "
  double mean_lambda(double lambdas[], double probs[]) {
    double res;
    double sum_probs = 0;
    // Normalize probs
    
    for (int i=0;i<3;i++) {
      sum_probs += probs[i];    
    }
    
    for (int i=0;i<3;i++) {
      probs[i] = probs[i]/sum_probs;
    }
    
    res = dot_product_2(lambdas, probs, 3);
    return res;
  }
  "
  
  # Vaccination
  vacc_param <- str_c(
    str_glue("double t_vacc_fake = {t_vacc_fake};"),
    vectorToCSnippet(t_vacc, 
                     vec_name = "t_vacc",
                     type = "double"),
    sep = "\n"
  )
  
  # E. logsumexp function ---
  # https://stackoverflow.com/questions/4169981/logsumexp-implementation-in-c
  logsumexp_C <- "
    double log_sum_exp(double nums[], const int ct) {
    double max_exp = nums[0], sum = 0.0;
    int i;

    for (i = 1 ; i < ct ; i++)
      if (nums[i] > max_exp)
        max_exp = nums[i];

    for (i = 0; i < ct ; i++)
      sum += exp(nums[i] - max_exp);

    return log(sum) + max_exp;
  }"
  
  # Function to determine the variant based on the the time
  # !! This function assumes that time starts on "2021-01-01" 
  get_variant_C <- str_glue("
    int get_variant(const double t) {{
      int ind;
    
      if(t < {dateToTime(dateStartAlpha())}) {{
        ind = 0;    // wildtype
      }} else if(t >= {dateToTime(dateStartAlpha())} && t < {dateToTime(dateStartDelta())}) {{
        ind = 1;    // alpha
      }} else if(t >= {dateToTime(dateStartDelta())} && t < {dateToTime(dateStartOmicron())}) {{
        ind = 2;    // delta
      }} else {{
        ind = 3;    // omicron
      }}
      
    return ind;
   
  }}")
  
  # Function to determine whether the participant was in the class or not at time t
  # !! This function assumes that time starts on "2021-01-01" 
  get_mask_C <- "
    double get_mask(
    double t, 
    const int n_mask, 
    const double start_mask[9],
    const double end_mask[9]) {
    double mask = 1;
    
    if (n_mask > 0) {
      for (int l=0; l<n_mask; l++) {
        if(t >= start_mask[l] && t <= end_mask[l]) {
          mask = 0;
        }
      }
    }
      
    return mask;
   
  }"
  
  # F. Parameter specifications ---
  # Set the same unit param for each lambda type
  set_expanded <- set_expanded %>% map_chr(function(x) {
    if (str_detect(x, "lambda")) {
      str_split(x, "_")[[1]][1]
    } else {
      x
    }
  }) %>% 
    unique()
  
  set_expanded <- paste0("const int ", 
                         set_expanded, 
                         "_unit = 1;\n", 
                         collapse = " ") %>% 
    Csnippet()
  
  set_fixed <- paste0("const int ", 
                      set_fixed, 
                      "_unit = 0;\n", 
                      collapse = " ") %>% 
    Csnippet()
  
  # G. likelihood flags ---
  lik_flags <- Csnippet(
    str_glue(
      "
      int use_test_lik = {use_test_lik};
      int use_sero_lik = {use_sero_lik};
      int use_gen_lik = {use_gen_lik};
    "
    )
  )
  
  # H. SARS-CoV-2 natural history parameters ---
  nat_hist <- str_c(
    vectorToCSnippet(vec = scov2_gamma,
                     vec_name = "scov2_gamma_vec"),
    vectorToCSnippet(vec = scov2_eta,
                     vec_name = "scov2_eta_vec"),
    vectorToCSnippet(vec = scov2_phi,
                     vec_name = "scov2_phi_vec"),
    matrixToCSnippet(mat = scov2_delta,
                     mat_name = "scov2_delta_mat"),
    vectorToCSnippet(vec = alpha,
                     vec_name = "alpha")
  )
  
  # I. Mask information
  mask_info <- str_c(
    vectorToCSnippet(vec = n_masks,
                     vec_name = "n_masks",
                     type = "int"),
    matrixToCSnippet(start_mask,
                     mat_name = "start_mask",
                     type = "double"),
    matrixToCSnippet(end_mask,
                     mat_name = "end_mask",
                     type = "double")
  )
  
  seroprev_info <- str_c(
    matrixToCSnippet(mat = seroprev_baseline,
                     mat_name = "seroprev_baseline",
                     type = "double"),
    matrixToCSnippet(mat = seroprev_middle,
                     mat_name = "seroprev_middle",
                     type = "double"),
    matrixToCSnippet(mat = seroprev_endline,
                     mat_name = "seroprev_endline",
                     type = "double"),
    vectorToCSnippet(vec = map_participant_block,
                     vec_name = "map_participant_block",
                     type = "int"),
    vectorToCSnippet(vec = map_participant_age_cat,
                     vec_name = "map_participant_age_cat",
                     type = "int"),
    str_glue("const int B = {max(map_participant_block) + 1};\n"),
    map_chr(seq_along(getSeroprevalenceDates()), function(x) {
      str_glue("const double t_{names(getSeroprevalenceDates())[x]} = {dateToTime(getSeroprevalenceFakeDates()[x])};")
    }) %>% 
      str_c(collapse = "\n"),
    str_glue("const int n_draws = {nrow(seroprev_endline)};"),
    sep = "\n"
  )
  
  seis_globals <- Csnippet(
    paste(
      test_perf,
      sens_sero_dat,
      dot_product_C,
      inv_logit_C,
      vec_mat,
      vec_mat_gen,
      logsumexp_C,
      # log1p_C,
      # log1m_C,
      mean_lambda_C,
      get_variant_C,
      get_mask_C,
      gen_param,
      set_fixed,
      set_expanded,
      lik_flags,
      nat_hist,
      mask_info,
      vacc_param,
      seroprev_info,
      sep = "\n"
    ))
  
  return(seis_globals)
}


# Parameter transformations -----------------------------------------------

buildParTrans <- function(basic_log_names,
                          basic_logit_names,
                          fixed_pars,
                          U) {
  
  basic_log_names <- setdiff(basic_log_names, fixed_pars)
  basic_logit_names <- setdiff(basic_logit_names, fixed_pars)
  
  log_names <- unlist(lapply(basic_log_names, 
                             function(x, U) paste0(x, 1:U), U))
  logit_names <- unlist(lapply(basic_logit_names, 
                               function(x, U) paste0(x, 1:U), U))
  
  seis_partrans <- pomp::parameter_trans(log = log_names, logit = logit_names)
  
  return(seis_partrans)
}



# Covariates --------------------------------------------------------------

buildCovariates <- function(seis_data) {
  seis_covar <- seis_data %>% 
    select(date, participant, seq_covar = seq) 
  
  seis_covar <- bind_rows(seis_covar,
                          seis_covar %>% distinct(participant) %>% 
                            mutate(date = 0, seq_covar = NA)) %>% 
    arrange(date, participant)
  
  return(seis_covar)
}

filterParNames <- function(x,
                           basic_params,
                           fixed_params) {
  x %>% 
    {.[. %in% basic_params]} %>% 
    {
      if (!is.null(fixed_params)) {
        setdiff(., fixed_params)
      } else {
        .
      }
    }
  
}

makeParLogNames <- function(basic_log_names =  c("lambda_wt", 
                                                 "lambda_a",
                                                 "lambda_d",
                                                 "lambda_o",
                                                 "lambdaA_wt", 
                                                 "lambdaA_a", 
                                                 "lambdaA_d", 
                                                 "lambdaA_o", 
                                                 "lambdaB_wt",
                                                 "lambdaB_a",
                                                 "lambdaB_d",
                                                 "lambdaB_o"),
                            basic_params,
                            fixed_params = NULL) {
  basic_log_names %>% 
    filterParNames(basic_params = basic_params,
                   fixed_params = fixed_params)
  
}

makeParLogitNames <- function(basic_logit_names = c("S_0", "C_0"), 
                              basic_params,
                              fixed_params = NULL) {
  basic_logit_names %>% 
    filterParNames(basic_params = basic_params,
                   fixed_params = fixed_params)
}

# Names -------------------------------------------------------------------

makeStateNames <- function(version = "epi") {
  
  if (version == "epi") {
    seis_unit_statenames <- c("S",         # Susceptible
                              "E",         # Exposed non-infectious
                              "Pone",      # Pre-symptomatic infectious
                              "Ptwo",      # 
                              "Pthree", 
                              "I",         # Symptomatic infectious
                              "Sprime",    # Susceptible with history of exposure
                              "C",         # Any infection/vaccination history
                              "X"          # Time of infection
    )
  } else if (version == "gen") {
    seis_unit_statenames <- c("S", "E", "Pone", "Ptwo", "Pthree", "I", "Sprime",
                              "C", "X", "Y", "genlik", "gencnt", "hazard", "ratio",
                              "C_block")
  }
  
  # seis_statenames <- paste0(rep(seis_unit_statenames, each = U), 1:U)
  
  return(seis_unit_statenames)
}


makeFullParamNames <- function(fixed_params,
                               expanded_params,
                               U) {
  
  expanded_params <- setdiff(expanded_params, fixed_params)
  
  full_paramnames <- c(
    if (length(fixed_params) > 0) {
      paste0(fixed_params, "1")
    }, 
    if (length(expanded_params) > 0) {
      paste0(rep(expanded_params, each = U), 1:U)
    })
  
  return(full_paramnames)
}

# Build pomp object -------------------------------------------------------

buildPompModel <- function(data,
                           covar = NULL,
                           times = "time",
                           units = "participant",
                           obsnames = c("n_pcr", "pcr_pos", "sero"),
                           t0,
                           model_version = "epi",
                           delta_t,
                           verbose = F,
                           init_params = c(
                             "lambda_wt" = .05,        # commnunity FOI
                             "lambda_a" = .05,       
                             "lambda_d" = .05,       
                             "lambda_o" = .05,        
                             "lambdaA_wt" = 15,        # same group FOI
                             "lambdaA_a" = 15,      
                             "lambdaA_d" = 15,        
                             "lambdaA_o" = 15,       
                             "lambdaB_wt" = 3,        # same class FOI
                             "lambdaB_a" = 3,       
                             "lambdaB_d" = 3,        
                             "lambdaB_o" = 3,       
                             "S_0" = .7,             # probability of susceptible at baseline
                             "C_0" = .3             # probability of previously infecged at baseline
                           ) ,
                           fixed_params = NULL,
                           U,
                           globals_parlist
                           
) {
  
  # A. Unit names ---
  unit_statenames <- makeStateNames(version = model_version)
  
  # B. Parameters ---
  # ivp_params <- c("S")
  # ivp_paramnames <- map(ivp_params, ~ paste0(., 1:U) %>% paste0("_0")) %>% unlist()
  
  # Set different parameters for within and across group transmission
  expanded_params <- names(init_params) %>% 
    setdiff(fixed_params) # names for expanded set of parameters to use ibpf
  
  full_paramnames <- makeFullParamNames(fixed_params = fixed_params,
                                        expanded_params = expanded_params,
                                        U = U)
  
  # Make initial parameters
  params <- makeParams(param_names = full_paramnames,
                       basic_params = init_params,
                       fixed = fixed_params,
                       expanded = expanded_params,
                       U = U
  )
  
  
  if (verbose) {
    cat("-- Pomp obj: Done parameters. \n")
  }
  
  # C. Parameter transformations ---
  par_trans <- buildParTrans(fixed_pars = fixed_params,
                             basic_log_names = makeParLogNames(basic_params = names(init_params),
                                                               fixed_params = fixed_params),
                             basic_logit_names = makeParLogitNames(basic_params = names(init_params),
                                                                   fixed_params = fixed_params),
                             U = U)
  
  if (verbose) {
    cat("-- Pomp obj: Done par_trans. \n")
  }
  
  # D. Globals ---
  if (model_version == "epi") {
    globals_parlist$gen <- makeNullGenGlobals()
  }
  
  globals_parlist$pars <- list(
    set_fixed = fixed_params,
    set_expanded = expanded_params
  )
  
  globals <- buildGlobalsWrapper(globals_parlist)
  
  if (verbose) {
    cat("-- Pomp obj: Done globals. \n")
  }
  
  if (!is.null(covar)) {
    seis_data <- spatPomp(
      data = data,
      times = times,
      units = units,
      obsnames = obsnames,
      covar = as.data.frame(covar),
      shared_covarnames = setdiff(colnames(covar), c(times, units)),
      t0 = t0)
  } else {
    seis_data <- spatPomp(
      data = data,
      times = times,
      units = units,
      obsnames = obsnames,
      t0 = t0)
  }
  
  if (verbose) {
    cat("-- Pomp obj: Done pomp data. \n")
  }
  
  # E. Pomp object ---
  seis_pomp <- spatPomp(
    seis_data,
    unit_statenames = unit_statenames,
    paramnames = full_paramnames,
    rinit = buildRinit(version = model_version),
    rprocess = euler(buildRProcess(version = model_version), delta.t = delta_t),
    dunit_measure = buildDUnitMeasure(version = model_version),
    runit_measure = buildRUnitMeasure(version = model_version),
    dmeasure = buildDMeasure(version = model_version),
    rmeasure = buildRMeasure(version = model_version),
    globals = globals,
    partrans = par_trans,
    verbose = verbose
  )
  
  
  if (verbose) {
    cat("-- Pomp obj: Done pomp obj. \n")
  }
  
  # F. Manual settings ---
  coef(seis_pomp) <- params
  seis_pomp@skeleton@skel.fn@statenames <- character()
  seis_pomp@dprior@statenames <- character()
  seis_pomp@rprior@statenames <- character()
  seis_pomp@vmeasure@statenames <- character()
  seis_pomp@emeasure@statenames <- character()
  seis_pomp@dprocess@statenames <- character()
  seis_pomp@accumvars <- character()
  
  # G. Fill covariate table to avoid issues with interpolation ---
  # if (!is.null(covar)) {
  #   
  #   for (i in 1:nrow(seis_pomp@covar@table)) {
  #     ind <- which(!is.na(seis_pomp@covar@table[i, ]))
  #     if (length(ind) > 0) {
  #       ind <- ind[1]
  #       seis_pomp@covar@table[i, 1:ind] <- seis_pomp@covar@table[i, ind]
  #     }
  #   }
  #   seis_pomp@covar@table[is.na(seis_pomp@covar@table)] <- -1
  # }
  
  if (verbose) {
    cat("-- Pomp obj: Done covars. \n")
  }
  
  return(seis_pomp)
}


# Utils -------------------------------------------------------------------


# From spatPomp
to_C_array <- function(v){
  paste0("{", 
         paste0(v, collapse = ","), 
         "}") 
}

matrixToCarray <- function(mat) {
  to_C_array(apply(mat, 1, to_C_array))
}

matrixToCSnippet <- function(mat, 
                             mat_name,
                             type = "double") {
  
  if (is.null(dim(mat))) {
    mat <- matrix(mat, ncol = 1)
  }
  
  if (length(dim(mat)) == 1) {
    mat <- matrix(mat, ncol = 1)
  }
  
  str_glue(
    "const {type} {mat_name}[{nrow(mat)}][{ncol(mat)}] = {matrixToCarray(mat)};"
  )
}

vectorToCSnippet <- function(vec, 
                             vec_name,
                             type = "double") {
  str_glue(
    "const {type} {vec_name}[{length(vec)}] = {to_C_array(vec)};"
  )
}

makeParams <- function(param_names, 
                       basic_params, 
                       fixed, 
                       expanded,
                       U) {
  res <- rep(0, length = length(param_names))
  names(res) <- param_names
  
  # Expand initial values even if in fixed
  if (any(fixed == "S_0") | any(fixed == "C_0")) {
    iv <- fixed[fixed %in% c("S_0", "C_0")]
    fixed <- setdiff(fixed, iv)
    expanded <- c(expanded, iv)
  }
  
  if (!is.null(fixed)) {
    for (p in fixed) res[paste0(p, 1)] <- basic_params[p]
  }
  if (!is.null(expanded)) {
    for (p in expanded) res[paste0(p, 1:U)] <- basic_params[p]
  }
  res
}


# Particle filter ---------------------------------------------------------


getSdRwParam <- function(refine) {
  case_when(
    refine ~ 1.5,
    TRUE ~ 3
  )
}

getSdRwIV <- function(refine) {
  case_when(
    refine ~ 0.25,
    TRUE ~ 0.5
  )
}

#' makeRWSD
#'
#' @param params 
#' @param fixed 
#' @param reg_param_sd 
#' @param ivp_param_sd 
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
makeRWSD <- function(params, 
                     fixed,
                     reg_param_sd = 0.005,
                     ivp_param_sd = 0.2,
                     times,
                     small_step_params = NULL) {
  
  # Make times between observations
  u_times <- unique(times) %>% sort()
  
  # Make time differences
  dt <- c(0, diff(u_times))
  sqrt_dt <- round(sqrt(dt)*1e3)*1e-3  * reg_param_sd
  u_sqrt_dt <- sort(unique(sqrt_dt))
  dt_list <- map_chr(
    1:length(u_sqrt_dt),
    function(x) {
      str_c(
        "c(",
        str_c(
          str_c("'",
                formatC(u_times[sqrt_dt == u_sqrt_dt[x]], digits = 5),
                "'"),
          collapse = ", "),
        ")"
      )
    })
  
  regular_pars <- names(params) %>% 
    str_subset("S_|C_", negate = T)
  
  ivp_pars <- names(params) %>% str_subset("S_|C_")
  
  if (length(fixed) > 0) {
    regular_pars <- regular_pars %>% 
      str_subset(str_c(str_c(fixed, "[0-9]+"), collapse = "|"), negate = T)
    
    ivp_pars <- ivp_pars %>% 
      str_subset(str_c(str_c(fixed, "[0-9]+"), collapse = "|"), negate = T)
  }
  
  if (length(regular_pars) > 0) {
    
    reg_param_txt <- map_chr(
      1:length(dt_list), 
      function(x) {
        if (x == length(dt_list)) {
          str_glue("ifelse(formatC(time, digits = 5) %in% {dt_list[x]}, {u_sqrt_dt[x]}, 0")
        } else {
          str_glue("ifelse(formatC(time, digits = 5) %in% {dt_list[x]}, {u_sqrt_dt[x]}, ")
        }
      }
    ) %>% 
      str_c(collapse = "") %>% 
      str_c(str_c(rep(")", length(u_sqrt_dt)), collapse = ""))
    
    if (!is.null(small_step_params)) {
      reg_small_step_param_txt <- map_chr(
        1:length(dt_list), 
        function(x) {
          if (x == length(dt_list)) {
            str_glue("ifelse(formatC(time, digits = 5) %in% {dt_list[x]}, {u_sqrt_dt[x]/5}, 0")
          } else {
            str_glue("ifelse(formatC(time, digits = 5) %in% {dt_list[x]}, {u_sqrt_dt[x]/5}, ")
          }
        }
      ) %>% 
        str_c(collapse = "") %>% 
        str_c(str_c(rep(")", length(u_sqrt_dt)), collapse = ""))
    }
    
    # Set the SDW of parameters by variant to be non-0 only when the variant
    # is circulating
    reg_pars_vect <- map_chr(
      regular_pars, 
      function(par) {
        # Get variant
        var <- getParamVariant(par)
        var_times <- getVariantDates(var) %>% 
          mutate(across(c("TL", "TR"), ~ dateToTime(.)))
        
        if (is.null(small_step_params)) {
          str_glue("{par} = ifelse(time >= {var_times$TL[1]} & time < {var_times$TR[1]}, {reg_param_txt}, 0)")
        } else {
          if (any(map_lgl(small_step_params, ~ str_detect(par, str_glue("^{.}[0-9]"))))) {
            str_glue("{par} = ifelse(time >= {var_times$TL[1]} & time < {var_times$TR[1]}, {reg_small_step_param_txt}, 0)")
          } else {
            str_glue("{par} = ifelse(time >= {var_times$TL[1]} & time < {var_times$TR[1]}, {reg_param_txt}, 0)")
          }
        }
      })
    
    reg_pars_txt <- str_c(reg_pars_vect, collapse = ", ")
    
  } else {
    reg_pars_txt <- character(0)
  }
  
  if (length(ivp_pars) > 0) {
    ivp_pars_txt <- str_c(
      ivp_pars, " = ivp(", ivp_param_sd, ")", collapse = ", "
    )
  } else {
    ivp_pars_txt <- character(0)
  }
  
  text <- str_c(
    "res <- rw_sd(",
    str_c(reg_pars_txt, ivp_pars_txt, sep = ","),
    ")")
  
  eval(parse(text = text))
  return(res)
}


# Indexing ----------------------------------------------------------------

getParamVariant <- function(par) {
  if (str_detect(par, "_wt")) {
    "wildtype"
  } else if (str_detect(par, "_a")) {
    "alpha"
  } else if (str_detect(par, "_d")) {
    "delta"
  } else if (str_detect(par, "_o")) {
    "omicron"
  } else {
    NA_character_
  }
}

getUnitSID <- function(unit, sugar_ids_mapping) {
  
}

getUnitSchool <- function(unit, 
                          sugar_ids_mapping,
                          date) {
  
}

getUnitGroup <- function(unit, 
                         sugar_ids_mapping,
                         date) {
  
}

getGroupOutbreakDates <- function(groups,
                                  terms,
                                  dt = 7) {
  
  outbreaks_df <- getOutbreakMetadata() %>% 
    filter(group %in% groups,
           term %in% terms)
  
  if (nrow(outbreaks_df) == 0) {
    stop("No data found for", groups, terms)
  }
  
  outbreaks_df %>% 
    select(TL = visite_1 - dt,
           TR = visiter_3 + dt)
  
}
