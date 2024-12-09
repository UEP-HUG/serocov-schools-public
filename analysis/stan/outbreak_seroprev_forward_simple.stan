// This is a model for individual-level seroprevalence using the forward equation
// at the individual level.

// For quadratic form based on Azman, A.S., et al.s, 2020. Vibrio cholerae O1 
// transmission in Bangladesh: insights from a nationally representative serosurvey. 
// The Lancet Microbe, 1(8), pp.e336-e343.
// https://github.com/HopkinsIDD/Bangladesh-Cholera-Serosurvey/blob/master/source/overall-analysis.stan


functions {
   // To integrate out variability
  real to_integrate(
    real x,          // Function argument
    real xc,         // Complement of function argument
    //  on the domain (defined later)
    array[] real theta,    // parameters
    array[] real x_r,      // data (real)
    array[] int x_i) {     // data (integer)
    
    real mu = theta[1];
    real sigma = theta[2];
    real res;
    
    if (size(x_r) == 1) {
      res = 1-exp(-exp(mu + sigma * inv_Phi(x)) * x_r[1]);
    } else {
      res = 1-exp(-exp(mu + sigma * inv_Phi(x)));
    }
    return res;
    }
  // Function to make covariates for cubic regression of sensitivity
  row_vector make_t_log_covar(real t, real t_log_mean, real t_log_sd) {
    row_vector[4] t_log;
    t_log[1] = 1;
    t_log[2] = (log(t) - t_log_mean)/t_log_sd;
    t_log[3] = t_log[2]^2;
    t_log[4] = t_log[2]^3;
    
    return t_log;
  }
}

data {
  // Integers
  int<lower=1> N;    // number of timepoints (can have no observations at baseline)
  int<lower=1> N_obs;    // number of serological observations
  int<lower=0> N_tests;    // number of observations of previously positives
  int<lower=1> M;    // number of participants
  int<lower=1> J;    // number of outbreaks
  // int<lower=1> K;    // number of terms
  // int<lower=1> S;    // number of schools
  int<lower=0> G;    // numger of school groups x outbreak combinations
  // int<lower=0> Z;    // number of covariates
  int<lower=0> L;    // total number of participants x intervals
  int<lower=0> V;    // total number of variant periods
  int<lower=0> max_n_intervals;
  int<lower=0> N_sens_pred;
  int<lower=0> E;    // total number of events (outbreaks, baselines)
  int<lower=0> n_sens_covar;
  int<lower=0> Z_comb;
  int<lower=0> n_covar_int;
  int<lower=0> L_outbreaks;    // number of intervals that correspond to outbreaks
  int<lower=0> N_all_tests;
  int<lower=0> Z_int_comb;    // number of unique combinations of interval covariates
  int<lower=0> L_cohorts;
  int<lower=1> M_cohorts;    // number of cohorts
  int<lower=1> N_cohorts;
  int<lower=1> N_age_cat;
  
  // Data
  array[N] real<lower=0> t_pos_test;      // number of positive tests in interval
  array[N] real<lower=1> times;   // map of interval to variant
  array[N] int<lower=0, upper=1> y;
  array[N_tests] int<lower=0, upper=1> prev_pos;
  array[N_all_tests] int<lower=0, upper=1> y_tests;
  array[N_cohorts] real<lower=1> cohort_times;
  
  // Maps
  array[M] int<lower=1, upper=L> starts;    // starts for each interval
  array[M] int<lower=1, upper=L> ends;      // ends for each interval
  array[M] int<lower=1,upper=N_age_cat> map_part_age_cat;    // map from participants to age categories
  array[L] int<lower=0> int_starts;    // starts for each interval
  array[L] int<lower=0> int_ends;      // ends for each interval
  array[L] int<lower=1, upper=V> map_int_variant;   // map of interval to variant
  array[L] int<lower=0, upper=G> map_int_outbreak_group;    // map from interval to outbreak x group
  array[L] int<lower=1,upper=Z_int_comb> map_int_int_comb; //// map from interval to interval combination
  array[L] int<lower=0, upper=L_outbreaks> map_int_obint;         // map from intervals to set of intervals with outbreaks
  array[L] int<lower=0, upper=1> map_int_vaccination;    // map from intervals to vaccination status
  array[N] int<lower=1, upper=M> map_obs_participant;
  array[N] int<lower=0, upper=L> map_obs_interval;    // map to intervals. If value is 0 means baseline observation with no prior interval
  // array[N] int<lower=1, upper=L> map_obs_time;
  // array[N_tests] int<lower=1> map_tests_int;
  array[N] int<lower=0, upper=E> map_obs_event;
  array[N_all_tests] int<lower=1, upper=V> map_test_variant;
  array[N_all_tests] int<lower=0, upper=G> map_test_outbreak;
  array[N_all_tests] int<lower=0, upper=L>  map_test_interval;
  array[N_all_tests] int<lower=0, upper=L_outbreaks> map_test_obint;         // map from intervals to set of intervals with outbreaks
  array[G] int<lower=0, upper=V> map_outbreak_variant;
  array[G] int<lower=1, upper=Z_int_comb> map_outbreak_group_int_comb;
  array[L_cohorts] int<lower=1, upper=V> map_cohort_int_variant;
  array[L_cohorts] int<lower=0, upper=G> map_cohort_int_outbreak_group;
  array[L_cohorts] int<lower=0> int_cohort_starts;    // starts for each interval
  array[L_cohorts] int<lower=0> int_cohort_ends;      // ends for each interval
  array[M_cohorts] int<lower=1, upper=L> starts_cohort;    // starts for each interval
  array[M_cohorts] int<lower=1, upper=L> ends_cohort;      // ends for each interval
  array[L_cohorts] int<lower=1,upper=Z_int_comb> map_cohort_int_int_comb; //// map from interval to interval combination
  
  // Indices
  array[N_obs] int<lower=1, upper=N> ind_sero_obs; 
  array[N_obs-N_tests] int<lower=1, upper=N> ind_only_sero; 
  array[N_tests] int<lower=1, upper=N> ind_prev_tests;
  
  // Controls
  int<lower=0> control_tp;
  int<lower=0> control_fp;
  int<lower=0> N_pos_control;
  int<lower=0> N_neg_control;
  real<lower=0> sens_pcr;
  real<lower=0> spec_pcr;
  
  // Generated data
  vector[N_sens_pred] dt_sens_pred; 
  int<lower=0> T_max;    // maximum delay after reported case
  
  // Priors
  real mu_log_lambda;
  real<lower=0> sd_log_lambda;
  real mu_log_alpha;
  real<lower=0> sd_log_alpha;
  
  // Covariates
  matrix[N, n_sens_covar] X_sens;    // covariate matrix for sensitivity
  matrix[Z_comb, n_sens_covar] X_sens_comb;    // unique combinations of sentitivity covariates
  matrix[L, n_covar_int]  X_int;    // Covariate matrix for intervals or variants/vacations
  matrix[Z_int_comb, n_covar_int] X_int_comb;    // matrix of unique covariate combinations
  matrix[G, n_covar_int] X_int_outbreak_group;    // covariate matrix to compute outbreaks
  matrix[L_cohorts, n_covar_int] X_cohort_int;
  
  // Other
  real<lower=0> dt_control;   // time delay from infection in controls (~ 11 months in Michielin et al. 2023)
  real<lower=0> inf_window;
  real<lower=0> inf_window_outbreak;    // infection window for outbreaks
  real<lower=0> dt_exposure;
  vector<lower=0,upper=1>[N] vaccinated;    // whether participant was vaccinated prior to the observation date
  
  // Baseline Dec 2020
  array[N_age_cat] int<lower=0> baseline_pos;
  array[N_age_cat] int<lower=0> N_baseline;
  real<lower=0,upper=1> sens_roche;    // Sensitivity of Roche-S test done in Dec 2022
  real<lower=0,upper=1> spec_roche;    // Sensitivity of Roche-S test done in Dec 2022
  real tol;    // tolerance for probability computation
  
}

transformed data {
  array[L] real<lower=0> dt;
  array[N] real<lower=0> dt_inf;    // time to infection
  array[L_cohorts] real<lower=0> dt_cohort;
  
  real dt_scale = 10;
  int Z = 4 + n_sens_covar;    // total number of covariates for sensitivity
  matrix[N, Z] X_sens_full;    // covariate matrix for sensitivity
  real dt_default = 30;
  vector[T_max] t_log_raw;
  real t_log_mean;
  real t_log_sd;
  matrix[N, 4] t_log_mat;    // values of log of delay from infection to serology
  row_vector[4] t_log_control;
  row_vector[4] t_log_default;
  row_vector[4] t_log_1;
  row_vector[4] t_log_500;
  
  matrix[N_sens_pred, 4] t_log_pred_mat;    // values of log of delay from infection to serology for prediction
  array[M] matrix<lower=0>[L, max_n_intervals] dt_int_inf;    // array of matrices of times to infection
  
  array[0] real x_r_mean_p_0;
  array[0] int x_i_mean_p_0;
  
  // -- Transdat.A: times to infection ----
  
  // ---- Transdat.A.1: unknown infection date ----
  // Define possible times to infection when only serology is available
  for(m in 1:M) {
    for (i in 1:L) {
      for (j in 1:max_n_intervals) {
        dt_int_inf[m][i, j] = 0;
      }
    }
  }
  
  for (i in ind_only_sero) {
    int m = map_obs_participant[i];    // this participant
    int l = map_obs_interval[i];       // this interval
    int n_int = max(0, l-starts[m]+1);           // number of intervals before this one to compute probabilities for
    int cnt = 0;                // counter of prob intervals
    
    // Set dt for all observations that are not the first
    if (l > 0) {
      // Compute interval probabilities of new infection
      for (j in starts[m]:l) {
        cnt += 1;    // update counter of intervals
        
        // Delay from this serology to the midpoint of the interval
        // dt_int_inf[m][l, cnt] = times[int_ends[l]] - (times[int_ends[j]] + times[int_starts[j]])/2;
        // Try to assign start of interval
        dt_int_inf[m][l, cnt] = times[int_ends[l]] - times[int_starts[j]];
        
        if (dt_int_inf[m][l, cnt] == 0) {
          dt_int_inf[m][l, cnt] = 1;
        }
        
      }
    }
  }
  
  // ---- Transdat.A.2: known infection date ----
  // Compute delays from pcr/antigen positive tests to serologies
  for (i in 1:N) {
    if (t_pos_test[i] > 0) {
      dt_inf[i] = (times[i] - t_pos_test[i]) + dt_exposure;
    } else {
      dt_inf[i] = 0;
    }
  }
  
  
  // -- Transdat.B: time-varying serology auxiliaries ----
  for(i in 1:T_max){
    t_log_raw[i] = log(i);
  }
  
  t_log_mean = sum(t_log_raw)/T_max;
  t_log_sd = sd(t_log_raw);
  
  // Initialize time deltas
  for (i in 1:L) {
    dt[i] = 0;
  }
  
  // Compute time deltas for each interval for computation of infection probability
  for (i in 1:L) {
    dt[i] = times[int_ends[i]] - times[int_starts[i]];
  }
  
  // Compute time deltas for each interval for computation of infection probability
  for (i in 1:L_cohorts) {
    dt_cohort[i] = cohort_times[int_cohort_ends[i]] - cohort_times[int_cohort_starts[i]];
  }
  
  
  for(i in 1:N){
    // int j = ind_prev_tests[i];
    t_log_mat[i, ] = make_t_log_covar(dt_inf[i], t_log_mean, t_log_sd);
  }
  
  // Concatenate
  if (n_sens_covar > 0) {
    X_sens_full = append_col(t_log_mat, X_sens);
  } else {
    X_sens_full = t_log_mat;
  }
  
  for(i in 1:N_sens_pred){
    t_log_pred_mat[i, ] = make_t_log_covar(dt_sens_pred[i], t_log_mean, t_log_sd);
  }
  
  // Covariates for control set
  t_log_control = make_t_log_covar(dt_control, t_log_mean, t_log_sd);
  t_log_default = make_t_log_covar(dt_default, t_log_mean, t_log_sd);
  t_log_1 = make_t_log_covar(1, t_log_mean, t_log_sd);
  t_log_500 = make_t_log_covar(500, t_log_mean, t_log_sd);
  
}

parameters {
  
  // Initial values at baseline
  vector[M] p_0_i_raw;    // individual level random effects
  vector[N_age_cat] p_0_mu;            // mean of probability of baseline prior infection by age category (logit-scale)
  vector<lower=0>[N_age_cat] p_0_sd;   // sd of probability of baseline prior infection by age category(logit-scale)
  
  // Outbreak exposure
  vector[G] log_alpha;
  // vector[L_outbreaks] alpha_i_raw;    // individual level outbreak random efects
  // vector<lower=0>[G] alpha_i_sd;      // sd of individual level random effects (REs)
  
  // Community exposures
  vector[n_covar_int] beta_int;      // parameters for interval covariates (variants and vacations)
  // vector[L] alpha_i_int_raw;         // individual level outbreak random efects for intervals
  // vector<lower=0>[Z_int_comb] alpha_i_int_sd;      // sd of individual level random effects (REs)
  
  // Serology tests
  real<lower=0, upper=1> spec;      // specificity
  vector[Z] beta;                   // regression coefficients of time-varying sensitivity, one for intercept and 3 for cubic regression
  
}

transformed parameters {
  vector[N] p;        // probability of historical infection at each time point. There are n_max_intervals + 1 available time points
  vector[M] p_0;      // initial probability of prior exposure
  vector[L] p_inf;    // probability of infection
  vector[L] eta_int;    // linear term of interval rate effects 
  real sens = inv_logit(t_log_control * beta[1:4]);    // serology sensivitiy for controls
  vector[N_age_cat] mean_p_0;    // expectation of baseline probability of prior infection
  
  // -- Transpar.A: Rate of community infection ----
  
  // Interval rate of infection from community 
  eta_int = exp(X_int * beta_int);
  
  
  // -- Transpar.B: Forward equation for probability of previous exposure ----
  
  // Compute the initial probability of exposure for all participants
  p_0 = inv_logit(p_0_mu[map_part_age_cat] + p_0_sd[map_part_age_cat] .* p_0_i_raw);
  
  // Approximation of expected value based on Deaunizeau (2017)
  // https://doi.org/10.48550/arXiv.1703.00091
  for (i in 1:N_age_cat) {
    mean_p_0[i] = inv_logit(p_0_mu[i]/sqrt(1 + 0.368 * p_0_sd[i]^2));
  }
  
  
  for(i in 1:M) {
    int cnt = 1;
    
    for (j in starts[i]:ends[i]) {
      int v = map_int_variant[j];   // this variant
      int ind_L = int_starts[j];
      int ind_R = int_ends[j];
      int o = map_int_outbreak_group[j];
      real hz;    // the hazard for computation of prob of infection in interval
      
      // Initialize probabilities
      if (cnt == 1) {
        p[ind_L] = p_0[i];   
      }
      
      if (o == 0) {
        hz = eta_int[j] * dt[j]/dt_scale;
      } else {
        int q = map_int_obint[j];    // index of the oubtreak interval (one per participant per outbreak)
        hz = eta_int[j] * dt[j]/dt_scale + exp(log_alpha[o]);
      }
      
      p_inf[j] = (1-exp(-hz));
      
      // Propagate in time. Here the infection probaility cannot go down
      p[ind_R] = p[ind_L] + (1-p[ind_L]) * p_inf[j];
      
      cnt += 1;
    }
  }
}
model {
  
  // For case where we only have serology
  // Need to marginalize out unknown dates of infection because here
  // were are assuming the sens at dt_control days.
  // y[ind_only_sero] ~ bernoulli(sens * p[ind_only_sero] + (1-spec) * (1-p[ind_only_sero]));
  real default_sens = inv_logit(t_log_default * beta[1:4]);
  vector[N_all_tests] p_inf_tests;    // the probability of infection prior to the test
  vector[N] p_sero;    // the vector of probability of true sero-positivity accounting for vaccination
  
  // First set the probabilities of seropositivity accounting for vaccination
  p_sero = p .* (1-vaccinated) + vaccinated;
  
  
  // -- Mod.A: All PCR/antigen test results ----
  // The test part of the likelihood 
  for (i in 1:N_all_tests) {
    int v = map_test_variant[i];
    int l = map_test_interval[i];
    int g = map_test_outbreak[i];
    int q = map_test_obint[i];
    real hz;
    
    if (g == 0) {
      // Only community infection
      hz = eta_int[l] * inf_window/dt_scale;
    } else {
      // Also account for outbreak
      hz = eta_int[l] * inf_window/dt_scale + exp(log_alpha[g]);
    }
    
    p_inf_tests[i] = 1-exp(-hz);
  }
  
  // Likelihood of test result
  y_tests ~ bernoulli(sens_pcr * p_inf_tests + (1-spec_pcr) * (1-p_inf_tests));
  
  
  // -- Mod.B: Serology  ----
  
  // ---- Mod.B.1: Serology only ---- 
  
  // For these participants we need to integrate out the unknown dates of 
  // possible infection to account for time-varying sensitivity.
  
  for (i in ind_only_sero) {
    int m = map_obs_participant[i];    // this participant
    int l = map_obs_interval[i];       // this interval
    int n_int = max(0, l-starts[m]+1);           // number of intervals before this one to compute probabilities for
    int cnt = 0;                // counter of prob intervals
    vector[n_int] p_int;        // probability of infection within each interval
    vector[n_int] cum_p_int;    // relatie probability of infection within each interval
    vector[n_int] dt_int;       // time delay to serology
    vector[n_int] ll;           // log-likelihoods of set of infection dates to marginalize
    real ll_inf;    // total log-likelihood of infection
    
    if (l == 0) {
      // This is for first observation per participant, assume constant delay of dt_default
      // TODO: change this to account for history
      target += bernoulli_lpmf(y[i]| default_sens * p_sero[i] + (1-spec) * (1-p_sero[i]));
      
    } else {
      // Marginalize over intervals
      // Compute interval probabilities of new infection
      for (j in starts[m]:l) {
        cnt += 1;    // update counter of intervals
        
        if (map_int_vaccination[j] == 1) {
          // Set to 1 in the event of vaccination
          p_int[cnt] = 1 - tol;
        } else {
          p_int[cnt] = p_inf[j];
        }
      }
      
      // Compute relative probabilities
      for (j in 1:n_int) {
        if (j == 1) {
          cum_p_int[j] = p_int[j] + tol;
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
        this_sens = this_sens_covar * beta;
        
        // Marginalize infection over previous intervals
        ll[j] = log(cum_p_int[j]/sum(cum_p_int)) + bernoulli_logit_lpmf(y[i]| this_sens);
      }
      
      // Total log-likelihood of infection at different dates 
      ll_inf = log_sum_exp(ll);
      
      // print("Index: ", i, " ll_inf: ", ll_inf, " log(p): ", log(p[i]), 
      // " log(1-p): ", log1m(p[i]) , " spec ", spec,
      // " ll ", bernoulli_lpmf(y[i]| (1-spec)));
      // 
      // Add to target the combination of infection and non-infection
      target += log_mix(p_sero[i], ll_inf,  bernoulli_lpmf(y[i]| (1-spec)));
    }
  }
  
  
  // ---- Mod.B.2: Previous PCR/antigen and serology ---- 
  // For these participants we can take a multinomial likelihood of the 
  // joint probability of serology and test result.
  {
    vector[2] p_comb;
    // 1. sero+ PCR+
    // 2. sero- PCR+
    for (i in 1:N_tests) {
      array[2] int res;
      int j = ind_prev_tests[i];
      real x = p_sero[j];
      // real this_sens = inv_logit(logit_sens + dt_inf[i] * beta);
      real this_sens = inv_logit(X_sens_full[j, ] * beta);
      real tot_pcomb  = 0;
      
      p_comb[1] = this_sens * sens_pcr * x + (1-spec) * (1-spec_pcr) * (1-x);
      p_comb[2] = (1-this_sens) * sens_pcr * x + spec * (1-spec_pcr) * (1-x);
      tot_pcomb = sum(p_comb);
      for (k in 1:2) {
        p_comb[k] = p_comb[k]/tot_pcomb;
      }
      // Put this in transformed data
      if (y[j] == 1) {
        res[1] = 1;
        res[2] = 0;
      } else {
        res[2] = 1;
        res[1] = 0;
      }
      
      target += multinomial_lpmf(res| p_comb);
    }
  }
  
  // -- Mod.C: Priors ----
  
  // ---- Mod.C.1: Community infections during intervals ----
  beta_int[1] ~ normal(mu_log_lambda, sd_log_lambda);
  beta_int[2:n_covar_int] ~ normal(0, 1);
  
  // ---- Mod.C.2: Outbreak infections ----
  log_alpha ~ normal(mu_log_alpha, sd_log_alpha);
  
  // ---- Mod.C.3: Serolgy sens/sepc controls ----
  control_tp ~ binomial(N_pos_control, sens);
  control_fp ~ binomial(N_neg_control, 1-spec);
  spec ~ beta(4, 1);
  
  // ---- Mod.C.4: Time-varying sensitivity ----
  // Time coefficients
  beta[1] ~ normal(2, .5);
  beta[3] ~ normal(-1, .5);
  
  // Prior on sensitivity at 1 day post exposure
  (t_log_1 * beta[1:4]) ~ normal(-3, .75);  // logit-scale, inv_logit(-3) ~= 0.075
  // Prior on sensitivity at 500 day post exposure
  (t_log_500 * beta[1:4]) ~ normal(0, .5);  // logit-scale, inv_logit(1) ~= 0.75
  
  // Other covariates
  if (Z > 4) {
    beta[5:Z] ~ normal(0, .5);
  }
  
  // ---- Mod.C.5: Initial exposure status ----
  p_0_mu ~ normal(-1.35, .5);    // logit-scale, inv_logit(-1.35) ~= 0.2
  p_0_sd ~ normal(0, 2.5);
  p_0_i_raw ~ std_normal();
  // Data from seroprevalence survey in December 2020
  baseline_pos ~ binomial(N_baseline, mean_p_0 * sens_roche + (1-spec_roche) * (1-mean_p_0));
}

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
        gen_seroprev_cohorts[ind_L] = integrate_1d(to_integrate, 0, 1, par, x_r, x_i, tol);
      }
      
      if (o == 0) {
        array[1] real x_r;
        array[0] int x_i;
        array[2] real par;
        par[1] = X_cohort_int[j, ] * beta_int;
        // par[2] = alpha_i_int_sd[k];
        x_r[1] = dt_cohort[j]/dt_scale;
        p_cohort_inf[j] = 1-exp(-exp(par[1]) * x_r[1]);
        // p_cohort_inf[j] = integrate_1d(to_integrate, 0, 1, par, x_r, x_i, tol);
      } else {
        array[1] real x_r_c;
        array[0] real x_r_o;
        array[0] int x_i;
        array[2] real par_c;
        array[2] real par_o;
        par_c[1] = X_cohort_int[j, ] * beta_int;
        // par_c[2] = alpha_i_int_sd[k];
        par_o[1] = log_alpha[o];
        // par_o[2] = alpha_i_sd[o];
        x_r_c[1] = inf_window/dt_scale;    // assuming outbreaks last inf_window
        p_cohort_inf[j] =  1-exp(-exp(par_c[1]) * x_r_c[1] - exp(par_o[1]));
        // p_cohort_inf[j] = 1 - integrate_1d(to_integrate2, 0, 1, par_c, x_r_c, x_i, tol) *  integrate_1d(to_integrate2, 0, 1, par_o, x_r_o, x_i, tol);
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
    // par[2] = alpha_i_int_sd[i];
    x_r[1] = 7/dt_scale;    // 7 days for weekly estimates
    // p_week_community[i] = integrate_1d(to_integrate, 0, 1, par, x_r, x_i, tol);
    p_week_community[i] = 1-exp(-exp(par[1]) * x_r[1]);
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
    // par_c[2] = alpha_i_int_sd[j];
    par_o[1] = log_alpha[i];
    // par_o[2] = alpha_i_sd[i];
    x_r_c[1] = inf_window_outbreak/dt_scale;    // 7 days for weekly estimates
    p_outbreak[i] = 1 - exp(-exp(par_c[1])*x_r_c[1] - exp(par_o[1]));
    // p_outbreak[i] = 1 - integrate_1d(to_integrate2, 0, 1, par_c, x_r_c, x_i, tol) *  integrate_1d(to_integrate2, 0, 1, par_o, x_r_o, x_i, tol);
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
