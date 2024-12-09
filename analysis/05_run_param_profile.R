# This script performs a global parameter search

# Preamble ----------------------------------------------------------------

library(stringr)
library(magrittr)
library(spatPomp)
library(optparse)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(tictoc)

source("analysis/utils.R")
source("analysis/pomp_utils.R")

option_list <- list(
  make_option(c("-f", "--run_file"), 
              default = "generated_data/run_data_4e535e_sb_sa_sc_sd_20_21_21_22_gen_genlik_true.rds", 
              action ="store", type = "character", help = "Processed run file"),
  make_option(c("-n", "--n_prof"), 
              default = 3, 
              action ="store", type = "numeric", help = "Number of profile points"),
  make_option(c("-m", "--n_rep"), 
              default = 2, 
              action ="store", type = "numeric", help = "Number of repetitions per profile point"),
  make_option(c("-i", "--init_ind"), 
              default = 1, 
              action ="store", type = "numeric", help = "Parameter set to use"),
  make_option(c("-p", "--profile_param"), 
              default = "lambda_a", 
              action ="store", type = "character", help = "Parameter to profile over"),
  make_option(c("-l", "--lower_bound"), 
              default = .1, 
              action ="store", type = "numeric", help = "Lower Bound"),
  make_option(c("-u", "--upper_bound"), 
              default = .2, 
              action ="store", type = "numeric", help = "Upper Bound"),
  make_option(c("-r", "--redo"), 
              default = FALSE, 
              action ="store", type = "logical", help = "Redo computations"),
  make_option(c("-x", "--maxlik_start"), 
              default = FALSE, 
              action ="store", type = "logical", help = "Start from max lik estimates")
)

# Deparse options
opt <- parse_args(OptionParser(option_list = option_list))

print(opt)

set.seed(opt$n_start)

# Make output file and directory names
# The pomp object file
pomp_file <- makePompObjFile(run_file = opt$run_file) 

# The ibpf profile object basename
pomp_profile_file <- makePompProfileFile(run_file = opt$run_file,
                                         param = opt$profile_param,
                                         lower = opt$lower_bound,
                                         upper = opt$upper_bound,
                                         n_prof = opt$n_prof,
                                         n_rep = opt$n_rep,
                                         maxlik_start = opt$maxlik_start) 

# The directory where to save results
pomp_profile_dir <- pomp_profile_file %>% 
  str_remove("\\.rds")
dir.create(pomp_profile_dir)

# The file name for this run id (used for multiple initial parameter runs)
this_pomp_profile_file <- str_c(
  pomp_profile_dir, "/",
  pomp_profile_file %>% 
    str_remove("generated_data/") %>% 
    str_replace(".rds", str_glue("_{opt$init_ind}.rds"))
)


# A. Load data ---------------------------------------------------------------

tic("Data setup")

# Get the run data built with 04_prepare_data_for_model.R
run_data <- readRDS(opt$run_file)

# Unpack some of the elements for further use
fixed_params <- run_data$fixed_param    # The fixed parameter that do not require optimizatoin
init_params <- run_data$init_params     # The default initial parameter
U <- run_data$sugar_ids_mapping %>% nrow()    # The total number of participants ('units' in the spatPomp jargon)

# The timestep for simulation, here set to 1 day. All times and rate parameters
# in the pomp object are specified in fractions of years from Jan 1st 2021.
delta_t <- 1/365  

# Compile pomp object (moved here for cluster compiling)
pomp_obj <- buildPompModel(
  data = as.data.frame(run_data$data),
  covar = as.data.frame(run_data$covar),
  times = "time",
  units = "participant",
  obsnames = c("n_pcr", "pcr_pos", "sero", "seq", "age_cat", "n_block"),
  t0 = min(run_data$data$time) - delta_t,
  model_version = run_data$opt$model_version,
  delta_t = delta_t,
  verbose = T,
  init_params = run_data$init_params,
  fixed_params = run_data$fixed_params,
  U = run_data$U,
  globals_parlist = run_data$globals_parlist
)

# names for expanded set of parameters to use ibpf
expanded_params <- setdiff(names(init_params), fixed_params)

# Full parameter names that can be subject to optimization (each unit can
# have a distinct parameter value in the bpf)
full_paramnames <- makeFullParamNames(fixed_params = fixed_params,
                                      expanded_params = expanded_params,
                                      U = U)

# These are parameters that are shared across units
shared_pars <- makeSharedParams(init_params = init_params,
                                fixed_params = fixed_params)

# Define the random walk sd for each parameter. We only perturb variant-specific
# parameters within the time window of variant circulation.
profile_bpf_rw <- makeRWSD(params = pomp_obj@params, 
                           fixed = c(fixed_params, opt$profile_param),
                           reg_param_sd = getSdRwParam(refine = TRUE), 
                           ivp_param_sd = getSdRwIV(refine = TRUE),
                           times = run_data$times,
                           small_step_params = c("lambdaA_a", "lambdaA_d", "lambdaA_o"))

toc(quiet = FALSE, func.toc = printTicTocTime)


# B. Make initial parameters -------------------------------------------------

if (opt$maxlik_start) {
  # Load search results to get relevant starting points for profiling
  max_lik_starts <- getMaxLikStarts(hash = extractHash(opt$run_file),
                                    lik_thresh = 15,
                                    redo = FALSE)
  
  init_params_df <- map_df(
    1:(opt$n_prof * opt$n_rep), 
    function(x) {
      
      # First choose what run to base sampling by sampling according to 
      # likelihood weight (exp(-ll_i/2)/sum(exp(-ll_j/2)))
      this_fit_id <- sample(max_lik_starts$fit_id, 1, prob = max_lik_starts$weight)
      
      res <- max_lik_starts %>% 
        filter(fit_id == this_fit_id) %>% 
        select(one_of(names(pomp_obj@params))) %>% 
        # apply perturbations
        mutate(across(everything(), ~ perturbInitParam(., w = .5))) %>% 
        # enforce constraints as initial values for C and S cannot exceed 1.
        mutate(across(contains("C"), ~ pmin(1, .)),
               across(contains("S"), ~ pmin(1, .)))
      
      res
    })
} else {
  
  # Base initial parameters on global search results
  init_params_df <- makeInitParams(n_start = opt$n_prof * opt$n_rep,
                                   refined = TRUE,
                                   fixed_params = fixed_params,
                                   init_params = init_params)
}


# Set the profiled parameter to a constant value across parameter starts
profile_values <- seq(opt$lower_bound, 
                      opt$upper_bound,
                      length.out = opt$n_prof) %>% 
  rep(times = opt$n_rep)

init_params_df <- init_params_df %>% 
  mutate(across(contains(opt$profile_param), function(x) {profile_values}))

# D. Run profiling --------------------------------------------------------

set.seed(getSeed(x = opt$init_ind))

# Select the initial parameter vector for this job array id
if (opt$maxlik_start) {
  
  # Select the initial parameter vector for this job array id
  params_start <- init_params_df %>% 
    slice(opt$init_ind) %>% 
    unlist()
  
} else {
  this_init_params <- init_params_df %>% 
    slice(opt$init_ind) %>% 
    pivot_longer(cols = everything()) %>% 
    deframe()
  
  # Expand parameters to all units
  params_start <- makeParams(param_names = full_paramnames,
                             basic_params = this_init_params,
                             fixed = fixed_params,
                             expanded = expanded_params,
                             U = U
  )
}


tic("IBPF")

if (!file.exists(this_pomp_profile_file) | opt$redo) {
  
  # Run the block particle filter
  res_ibpf <- ibpf(
    pomp_obj,
    params = params_start,
    Nbpf = 100,
    Np = 2.5e3,
    rw.sd = profile_bpf_rw,
    unitParNames = c("C_0"),
    sharedParNames = shared_pars,
    block_list = run_data$block_list,
    cooling.fraction.50 = .5,
    spat_regression = .5,
    verbose = T
  )
  
  
  # Save in case next fails
  saveRDS(
    list(
      res_ibpf = res_ibpf
    ),
    file = this_pomp_profile_file
  )
  
} else {
  res_ibpf <- readRDS(this_pomp_profile_file)$res_ibpf
  
  if (is.null(res_ibpf)) {
    stop("Computation of ll already done.")
  }
  
  cat("Loaded res_ibpf from file. \n")
}

toc(quiet = FALSE, func.toc = printTicTocTime)


tic("Log-lik estimation")
# Compute estimates of the log-likelihood. We need multiple estimates to compute
# a mean and sd.
bpf_list <- replicate(
  5,
  bpfilter(res_ibpf,
           Np = 7.5e3,
           block_list = run_data$block_list,
           verbose = F)
) 

toc(quiet = FALSE, func.toc = printTicTocTime)

# E. Unpack ------------------------------------------------------------------

# Extract the log-likelihoods
bpf_ll <- purrr::map(bpf_list, ~spatPomp::logLik(.))

# Compute mean and sd
ll_stats <- pomp::logmeanexp(bpf_ll, se = T)

# Compute the block mean and variance
block_stats <- computeBlockLikStats(bpf_list)

# Pack results for further exploration
res_df <- tibble(run = opt$init_ind,
                 ll = ll_stats["est"],
                 ll_se = ll_stats["se"]) %>% 
  bind_cols(res_ibpf@params %>% 
              enframe() %>% 
              pivot_wider())

# F. Save results ------------------------------------------------------------

saveRDS(
  list(
    # res_ibpf = res_ibpf,
    res_df = res_df,
    bpf_ll = bpf_ll,
    # bpf_list = bpf_list,
    init_params = params_start,
    # profile_bpf_rw = profile_bpf_rw,
    block_stats = block_stats
  ),
  file = this_pomp_profile_file
)
