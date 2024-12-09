# This script draws from the smoothing distribution at given parameter values

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
  make_option(c("-s", "--hash"), 
              default = "4e53", 
              action ="store", type = "character", help = "Processed run file"),
  make_option(c("-n", "--n_traj"), 
              default = 1, 
              action ="store", type = "numeric", help = "Number of parameter sets to explore"),
  make_option(c("-i", "--this_draw"), 
              default = 1, 
              action ="store", type = "numeric", help = "This sample"),
  make_option(c("-d", "--redo"), 
              default = TRUE, 
              action ="store", type = "logical", help = "Whether to restart run based on previous results.")
)


# Deparse options
opt <- parse_args(OptionParser(option_list = option_list))

set.seed(as.numeric(str_c(opt$this_draw, opt$n_traj)))

# Make output file and directory names

# Get fit directory
fit_dir <- getResDir(opt$hash, what = "fit")

# Make the trajectory directory
pomp_traj_dir <- str_replace(fit_dir, "pomp_fit", "pomp_traj")
dir.create(pomp_traj_dir)

# The trajectory object basename
pomp_traj_file <- makePompTrajFile(run_file = getRunDataFile(opt$hash))

# The file name for this run id (used for multiple initial parameter runs)
this_pomp_traj_file <- str_c(
  pomp_traj_dir, "/",
  pomp_traj_file %>% 
    str_remove("generated_data/") %>% 
    str_replace(".rds", str_glue("_{opt$this_draw}.rds"))
)

# Load fit results --------------------------------------------------------

# Load run data
run_data <- readRDS(getRunDataFile(opt$hash))

# Get fit results to draw from
maxlik_res <- getMaxLikStarts(hash = opt$hash,
                              lik_thresh = 10,
                              redo = TRUE)

# Select one random draw
this_fit_id <- sample(maxlik_res$fit_id, 1, prob = maxlik_res$weight)

res_ibpf <- readRDS(maxlik_res$path[this_fit_id])$res_ibpf

this_params <- maxlik_res %>% 
  filter(fit_id == this_fit_id) %>% 
  select(one_of(names(res_ibpf@params))) %>%
  unlist()

# Run filtering -----------------------------------------------------------

set.seed(getSeed(x = opt$this_draw))

tic("Filtering")

filter_trajs <- purrr::map(
  1:opt$n_traj, 
  function(x) {
    
    res <- bpfilter(res_ibpf,
                    block_list = run_data$block_list,
                    filter_traj = T,
                    Np = 25e3,
                    verbose = F)
    
    res@filter.traj
  })

toc(quiet = FALSE, func.toc = printTicTocTime)


# Save results ------------------------------------------------------------

saveRDS(
  filter_trajs,
  file = this_pomp_traj_file
)

