# This script performs a global parameter search.

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

# utility functions
source("analysis/utils.R")
source("analysis/pomp_utils.R")

# options for Rscript 
option_list <- list(
  make_option(c("-f", "--run_file"), 
              default = "generated_data/run_data_4e535e_sb_sa_sc_sd_20_21_21_22_gen_genlik_true.rds", 
              action ="store", type = "character", help = "Processed run file"),
  make_option(c("-n", "--n_start"), 
              default = 1, 
              action ="store", type = "numeric", help = "Number of parameter sets to explore"),
  make_option(c("-i", "--init_ind"), 
              default = 1, 
              action ="store", type = "numeric", help = "Parameter set to use"),
  make_option(c("-r", "--refine"), 
              default = TRUE, 
              action ="store", type = "logical", help = "Whether to try a refinement run based preliminary global run and profiling."),
  make_option(c("-x", "--restart"), 
              default = FALSE, 
              action ="store", type = "logical", help = "Whether to restart run based on previous results."),
  make_option(c("-d", "--redo"), 
              default = FALSE, 
              action ="store", type = "logical", help = "Whether to restart run based on previous results.")
)

# Deparse options
opt <- parse_args(OptionParser(option_list = option_list))

print(opt)

set.seed(opt$n_start)

# Make output file and directory names
# The pomp object file
pomp_file <- makePompObjFile(run_file = opt$run_file) 

# The ibpf object basename
pomp_fit_file <- makePompFitFile(run_file = opt$run_file, refine = opt$refine)

# The directory where to save results
pomp_fit_dir <- pomp_fit_file %>% str_remove("\\.rds")
dir.create(pomp_fit_dir)

# The file name for this run id (used for multiple initial parameter runs)
this_pomp_fit_file <- str_c(
  pomp_fit_dir, "/",
  pomp_fit_file %>% 
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
delta_t <- getDeltaT()  

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
search_bpf_rw <- makeRWSD(params = pomp_obj@params, 
                          fixed = fixed_params,
                          reg_param_sd = getSdRwParam(refine = opt$restart | opt$refine), 
                          ivp_param_sd = getSdRwIV(refine = opt$restart | opt$refine),
                          times = run_data$times,
                          small_step_params = c("lambdaA_a", "lambdaA_d", "lambdaA_o"))

toc(quiet = FALSE, func.toc = printTicTocTime)

# B. Make initial parameters -------------------------------------------------

if (opt$restart)  {
  
  # If refining previous run load previous results
  init_params_df <- getMaxLikStarts(hash = extractHash(opt$run_file),
                                    lik_thresh = NULL,
                                    redo = FALSE)
} else {
  
  # Base initial parameters on global search results
  init_params_df <- makeInitParams(n_start = opt$n_start,
                                   refined = opt$refine,
                                   fixed_params = fixed_params,
                                   init_params = init_params)
}


# C. Run iterated filtering -----------------------------------------------

set.seed(getSeed(x = opt$init_ind))

# Select the initial parameter vector for this job array id
if (opt$restart) {
  # If refining no need to expand 
  params_start <- init_params_df %>% 
    slice(opt$init_ind) %>% 
    select(-run, -contains("ll")) %>% 
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

if (!file.exists(this_pomp_fit_file) | opt$redo) {
  tic("IBPF")
  # Run the block particle filter
  res_ibpf <- ibpf(pomp_obj,
                   params = params_start,
                   Nbpf = 100,
                   Np = 2.5e3,
                   rw.sd = search_bpf_rw,
                   unitParNames = c("C_0"),
                   sharedParNames = shared_pars,
                   block_list = run_data$block_list,
                   cooling.fraction.50 = .5,
                   spat_regression = .5,
                   verbose = T
  )
  
  toc(quiet = FALSE, func.toc = printTicTocTime)
  
  # Save in case next fails
  saveRDS(
    list(
      res_ibpf = res_ibpf,
      params_start = params_start,
      search_bpf_rw = search_bpf_rw
    ),
    file = this_pomp_fit_file
  )
} else {
  bpf_ll <- readRDS(this_pomp_fit_file)$bpf_ll
  
  if (!is.null(bpf_ll)) {
    stop("Computation of ll already done.")
  }
  
  res_ibpf <- readRDS(this_pomp_fit_file)$res_ibpf
  
  cat("Loaded res_ibpf from file. \n")
}


tic("Log-lik estimation")
# Compute estimates of the log-likelihood. We need multiple estimates to compute
# a mean and sd.
bpf_list <- replicate(
  5,
  bpfilter(res_ibpf,
           Np = 5e3,
           block_list = run_data$block_list,
           verbose = F)
) 

toc(quiet = FALSE, func.toc = printTicTocTime)

# D. Unpack results -------------------------------------------------------

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

# E. Save results -----------------------------------------------------------

saveRDS(
  list(
    res_ibpf = res_ibpf,
    res_df = res_df,
    bpf_ll = bpf_ll,
    bpf_list = bpf_list,
    params_start = params_start,
    block_stats = block_stats,
    search_bpf_rw = search_bpf_rw
  ),
  file = this_pomp_fit_file
)
