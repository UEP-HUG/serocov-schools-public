generated quantities {
  vector[Z_int_comb] p_week_community;
  vector[G] p_outbreak;
  matrix[N_sens_pred, max(Z_comb, 1)] sens_pred;
  vector[N] gen_sero;
  vector[E] gen_sero_event_all;
  vector[E] gen_sero_event_obs;
  vector[E] gen_seroprev_event_all;
  vector[E] gen_seroprev_event_obs;
  real default_sens = inv_logit(t_log_default * beta[1:4]);
  vector[N_cohorts] gen_seroprev_cohorts;
  vector[L_cohorts] p_cohort_inf;    // probability of infection
  
  
  // Generate data for cohorts
  for(i in 1:M_cohorts) {
    int cnt = 1;    // counter of intervals per cohort
    
    for (j in starts_cohort[i]:ends_cohort[i]) {
      int v = map_cohort_int_variant[j];   // this variant
      int ind_L = int_cohort_starts[j];
      int ind_R = int_cohort_ends[j];
      int o = map_cohort_int_outbreak_group[j];
      real hz;    // the hazard for computation of prob of infection in interval
      int k = map_cohort_int_int_comb[j];
      
      // Initialize probabilities
      if (cnt == 1) {
        array[0] real x_r;
        array[0] int x_i;
        array[2] real par;
        par[1] = p_0_mu[1];
        par[2] = p_0_sd[1];
        gen_seroprev_cohorts[ind_L] = integrate_1d(to_integrate, 0, 1, par, x_r, x_i, 1e-6);
      }
      
      if (o == 0) {
        array[1] real x_r;
        array[0] int x_i;
        array[2] real par;
        par[1] = X_cohort_int[j, ] * beta_int;
        par[2] = alpha_i_int_sd[k];
        x_r[1] = dt_cohort[j]/dt_scale;
        p_cohort_inf[j] = integrate_1d(to_integrate, 0, 1, par, x_r, x_i, 1e-6);
      } else {
        array[1] real x_r_c;
        array[0] real x_r_o;
        array[0] int x_i;
        array[2] real par_c;
        array[2] real par_o;
        par_c[1] = X_cohort_int[j, ] * beta_int;
        par_c[2] = alpha_i_int_sd[k];
        par_o[1] = log_alpha[o];
        par_o[2] = alpha_i_sd[v];
        x_r_c[1] = inf_window/dt_scale;    // assuming outbreaks last inf_window
        p_cohort_inf[j] = 1 - integrate_1d(to_integrate2, 0, 1, par_c, x_r_c, x_i, 1e-6) *  integrate_1d(to_integrate2, 0, 1, par_o, x_r_o, x_i, 1e-6);
      }
      
      
      // Propagate in time qith forward equation
      gen_seroprev_cohorts[ind_R] = gen_seroprev_cohorts[ind_L] + (1-gen_seroprev_cohorts[ind_L]) * p_cohort_inf[j];
      
      cnt += 1;
    }
  }
  
  
  // -- Gen.A: Weekly community infection probablity ----
  for (i in 1:Z_int_comb) {
    array[1] real x_r;
    array[0] int x_i;
    array[2] real par;
    par[1] = X_int_comb[i, ] * beta_int;
    par[2] = alpha_i_int_sd[i];
    x_r[1] = 7/dt_scale;    // 7 days for weekly estimates
    p_week_community[i] = integrate_1d(to_integrate, 0, 1, par, x_r, x_i, 1e-6);
  }
  
  // -- Gen.B: Outbreak infection probability ----
  for (i in 1:G) {
    array[1] real x_r_c;
    array[0] real x_r_o;
    array[0] int x_i;
    array[2] real par_c;
    array[2] real par_o;
    int j = map_outbreak_group_int_comb[i];
    int v = map_outbreak_variant[i];
    par_c[1] = X_int_outbreak_group[i, ] * beta_int;
    par_c[2] = alpha_i_int_sd[j];
    par_o[1] = log_alpha[i];
    par_o[2] = alpha_i_sd[v];
    x_r_c[1] = inf_window/dt_scale;    // 7 days for weekly estimates
    p_outbreak[i] = 1 - integrate_1d(to_integrate2, 0, 1, par_c, x_r_c, x_i, 1e-6) *  integrate_1d(to_integrate2, 0, 1, par_o, x_r_o, x_i, 1e-6);
  }
  
  
  // -- Gen.C: Predicted time-varying serology senstivity ----
  if (Z_comb == 0) {
    for (i in 1:N_sens_pred) {
      sens_pred[i, 1] = inv_logit(t_log_pred_mat[i, ] * beta[1:4]);
    }
  } else {
    for (z in 1:Z_comb) {
      for (i in 1:N_sens_pred) {
        sens_pred[i, z] = inv_logit(append_col(t_log_pred_mat[i, ], X_sens_comb[z, ]) * beta);
      }
    }
  }
  
  
  // -- Gen.D: Posterior retrodictive checks ---- 
  // ---- Gen.D.1: Serology results for all time points ---- 
  for (i in 1:N) {
    // Generate with control sensitivity
    gen_sero[i] = bernoulli_rng(sens * p[i] + (1-spec) * (1-p[i]));
  }
  
  // ---- Gen.D.2: Serology results for serology-only time points ---- 
  // For individuals without test data compute average
  for (i in ind_only_sero) {
    int m = map_obs_participant[i];    // this participant
    int l = map_obs_interval[i];       // this interval
    int n_int = max(0, l-starts[m]+1);           // number of intervals before this one to compute probabilities for
    int cnt = 0;                // counter of prob intervals
    vector[n_int] p_int;        // probability of infection within each interval
    row_vector[n_int] cum_p_int;    // relatie probability of infection within each interval
    row_vector[n_int] this_gen_sero;           // log likelihoods to marginalize
    
    
    if (l == 0) {
      // This is for first observation per participant, assume constant delay of dt_default
      gen_sero[i] = bernoulli_rng(default_sens * p[i] + (1-spec) * (1-p[i]));
      
    } else {
      // Compute mean of generated quantity
      // Compute interval probabilities of new infection
      for (j in starts[m]:l) {
        cnt += 1;    // update counter of intervals
        // Add tolerance of 1e-8 to avoid numerical underflows
        p_int[cnt] = p_inf[j];
      }
      
      // Compute relative probabilities
      for (j in 1:n_int) {
        if (j == 1) {
          cum_p_int[j] = p_int[j] + 1e-7;
        } else {
          cum_p_int[j] = p_int[j] * prod(1-p_int[1:(j-1)]) + 1e-7;
        }
      }
      
      // Now marginalize
      for (j in 1:n_int) {
        real n_int2 = n_int;
        real this_sens;
        row_vector[Z] this_sens_covar;
        
        // Thes are the time covariates for cubic regression
        this_sens_covar[1:4] = make_t_log_covar(dt_int_inf[m][l, j], t_log_mean, t_log_sd);
        
        if (n_sens_covar > 0) {
          this_sens_covar[5:Z] = X_sens_full[i, 5:Z];
        }
        
        // Compute sensitivity
        this_sens = inv_logit(this_sens_covar * beta);
        
        // Marginalize infection over previous intervals
        this_gen_sero[j] = bernoulli_rng(this_sens * p[i] + (1-spec) * (1-p[i]));
      }
      
      // Compute weighted average of generated outcome
      gen_sero[i] = dot_product(cum_p_int/sum(cum_p_int), this_gen_sero);
    }
  }
  
  // ---- Gen.D.3: Serology results for serology+test time points ---- 
  for (i in 1:N_tests) {
    int j = ind_prev_tests[i];
    real this_sens = inv_logit(X_sens_full[j, ] * beta);
    gen_sero[j] =  bernoulli_rng(this_sens * p[j] + (1-spec) * (1-p[j]));
  }
  
  // ---- Gen.D.4: Aggregate by event and compute seroprevalence ---- 
  // We need to differentiate between the total number of results for all participants
  // and for those that were actually observed.
  {
    // Initialize counters of number of participants to compute seroprevalences
    vector[E] cnt_all = rep_vector(0, E);
    vector[E] cnt_obs = rep_vector(0, E);
    
    for (i in 1:E) {
      gen_sero_event_all[i] = 0;
      gen_sero_event_obs[i] = 0;
    }
    
    // Generat by event for all participants with any data (not only those observed in event)
    for (i in 1:N) {
      int l = map_obs_event[i];
      if (l > 0) {
        gen_sero_event_all[l] += gen_sero[i];
        cnt_all[l] += 1;
      }
    }
    
    // Generat by event for observed values
    for (i in 1:N_obs) {
      int j = ind_sero_obs[i];
      int l = map_obs_event[j];
      if (l > 0) {
        gen_sero_event_obs[l] += gen_sero[j];
        cnt_obs[l] += 1;
      }
    }
    
    // Now compute seroprevalences
    gen_seroprev_event_all = gen_sero_event_all ./ cnt_all;
    gen_seroprev_event_obs = gen_sero_event_obs ./ cnt_obs;
  }
}