# This script makes the simulations for the intervention comparison


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
  make_option(c("-c", "--config"), 
              default = "./analysis/configs/preschool_omicron_like_w1_b1.yml", 
              action ="store", type = "character", help = "config for this simulation")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load scenario config
config <- yaml::read_yaml(opt$config)

# File names --------------------------------------------------------------

# Make simulation directory
dir.create("generated_data/simulations")
sim_dir <- str_extract(opt$config, "(.)*(?=_w)") %>% 
  str_replace("analysis/configs", "generated_data/simulations")
dir.create(sim_dir)

# This file name
this_config <- str_split(opt$config, "configs/")[[1]][2]
this_sim_file <- str_c(
  sim_dir, "/",
  this_config %>% 
    str_replace("yml", "rds")
)

# Make synthetic school data ----------------------------------------------
run_data <- makeSimRunData(config = config)


# Params ------------------------------------------------------------------

# Apply reductions
this_params <- applyInterventions(params = config$param_specs,
                                  interventions = config$intervention_specs)

init_params = c(
  replicateParams(params = this_params),
  "S_0" = .9,        # probability of susceptible at baseline
  "C_0" = .1         # probability of susceptible at baseline
)  

base_pars <- names(init_params)
fixed_params <- "S_0"


# Prepare pomp object -----------------------------------------------------

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
  model_version = "gen",
  delta_t = delta_t,
  verbose = T,
  init_params = init_params,
  fixed_params = fixed_params,
  U = run_data$U,
  globals_parlist = run_data$globals_parlist
)



# Run simulations ---------------------------------------------------------

tic("Run scenarios")

sims <- runSimulations(pomp_obj,
                       nsim = config$sim_specs$n_sim,
                       step = config$sim_specs$step_size,
                       tmin = dateToTime(config$sim_specs$date_start),
                       tmax = dateToTime(config$sim_specs$date_end))

toc(quiet = FALSE, func.toc = printTicTocTime)

agg_sim <- sims %>% 
  aggregateSim(seis_unit_statenames = makeStateNames(version = "gen"),
               by_school = F,
               by_group = F,
               simplify_states = T)

saveRDS(
  list(
    agg_sim = agg_sim,
    config = config
  ),
  file = this_sim_file
)
