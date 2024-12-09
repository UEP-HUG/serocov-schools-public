# This script computes statistics of interest from the smoothing distribution

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
library(forcats)

# utility functions
source("analysis/utils.R")
source("analysis/pomp_utils.R")

# options for Rscript 
option_list <- list(
  make_option(c("-s", "--hash"), 
              default = "4e53", 
              action ="store", type = "character", help = "Processed run file"),
  make_option(c("-r", "--redo"), 
              default = FALSE, 
              action ="store", type = "logical", help = "Whether to restart run based on previous results."),
  make_option(c("-o", "--outbreak_id"), 
              default = NULL, 
              action ="store", type = "numeric", help = "Whether to extract detailed data for one outbreak")
)


# Deparse options
opt <- parse_args(OptionParser(option_list = option_list))

# Make output file and directory names
# Get fit directory
fit_dir <- getResDir(opt$hash, what = "fit")

# Make the trajectory directory
pomp_traj_dir <- str_replace(fit_dir, "pomp_fit", "pomp_traj")

# The trajectory object basename
pomp_traj_file <- makePompTrajFile(run_file = getRunDataFile(opt$hash),
                                   outbreak_id = opt$outbreak_id)

# The file name for this run id (used for multiple initial parameter runs)
pomp_traj_stats_file <- str_c(
  pomp_traj_dir, "/",
  pomp_traj_file %>% 
    str_remove("generated_data/") %>% 
    str_replace(".rds", "_stats.rds")
)

# Get all the files excluding the eventual stats file already computed
traj_files <- dir(pomp_traj_dir, full.names = TRUE) %>% 
  str_subset("stats|compiled", negate = T)

pomp_traj_compiled_file <- pomp_traj_stats_file %>% 
  str_replace("stats", "compiled") 

model_fidelity_file <- pomp_traj_stats_file %>% 
  str_replace("stats", "fidelity") 

if (is.null(opt$outbreak_id)) {
  cat("-- Computing traj stats for all times and schools \n")
} else {
  cat(str_glue("-- Computing traj stats for outbreak {opt$outbreak_id} \n"))
}

# Load data ---------------------------------------------------------------

# Run data
run_data <- readRDS(getRunDataFile(opt$hash))

# School groups
this_school_groups <- makeSchoolGroupMap(run_data = run_data)

# Filtering times
times <- run_data$data$time %>% unique() %>% sort()


# Compute statistics ------------------------------------------------------

# Variables for which to extract trajectories for
if (is.null(opt$outbreak_id)) {
  
  # Variables for figures
  traj_variables <- c("C", "hazard", "ratio")
  
  # Extract data for all schools and times 
  this_school <- NULL
  time_min <- NULL
  time_max <- NULL
  
} else {
  
  # Only save state variables
  traj_variables <- c("C", "Pone", "Ptwo", "Pthree", "I")
  
  # Get outbreak info
  outbreak_metadata <- getOutbreakMetadata() %>% 
    filter(outbreak_id == opt$outbreak_id)
  
  this_school <- outbreak_metadata$school[1]
  time_min <- dateToTime(min(outbreak_metadata$visite_1) - 7)
  time_max <- dateToTime(max(outbreak_metadata$visite_3) + 7)
  
}


tic("Trajectory extraction")

if (!file.exists(pomp_traj_compiled_file) | opt$redo) {
  
  # Loop over trajectory files and extract trajectories for selected variables
  traj_df <- map_df(
    seq_along(traj_files), 
    function(i) {
      
      # Load trajectory file
      res_list <- readRDS(traj_files[i])
      
      # Extract trajectories
      trajs <- extractTraj(res_list = res_list,
                           to_keep = traj_variables,
                           school_groups = this_school_groups,
                           times = times,
                           time_min = time_min,
                           time_max = time_max,
                           school = this_school) %>% 
        mutate(traj_file = i)
      
      cat("-- Done traj file", i, "/", length(traj_files), ":\t", traj_files[i], "\n")
      
      trajs
    }) %>% 
    # Update draw to be unique to each file
    mutate(draw = str_c(traj_file, draw, sep = "-"))
  
  # Add age categories
  traj_df <- traj_df %>% 
    inner_join(run_data$full_sugar_ids %>% 
                 distinct(sugar_id, age_cat))
  
  # Save for further use
  saveRDS(traj_df, file = pomp_traj_compiled_file)
  
  
  ## Extract model fidelity -----
  if (!is.null(opt$outbreak_id)) {
    # Compute stats for each date for which data is available
    data_by_date <- getDataOutbreak(ob_ids = outbreak_id) %>% 
      filter(sugar_id %in% traj_df$sugar_id) %>% 
      group_by(what, date, school, group) %>% 
      summarise(n = n(),
                who = list(unique(sugar_id)),
                n_pos = sum(result == 1))
    
    # Filter dates from trajectories
    sim_stats <- map_df(unique(data_by_date$what),
                        function(x) {
                          dat <- data_by_date %>% filter(what == x)
                          computeOutbreakStatsByDate(
                            traj_df = traj_df,
                            dates = dat$date,
                            ids = dat$who
                          ) %>% 
                            mutate(what = ifelse(what == "sim_infected", "pcr/antigen", "sero")) %>% 
                            filter(what == x)
                        })
    
    
    model_fidelity <- data_by_date %>% 
      inner_join(sim_stats %>% 
                   group_by(what, date, group) %>% 
                   slice(1) %>% 
                   ungroup()) %>% 
      distinct() %>% 
      mutate(outbreak_id = outbreak_id)
    
    # Save
    saveRDS(model_fidelity, file = model_fidelity_file)
  }
} else {
  traj_df <- readRDS(pomp_traj_compiled_file)
}

toc(quiet = FALSE, func.toc = printTicTocTime)


tic("Stat computation")

# Compute statistics
traj_stats <- map_df(
  traj_variables, 
  function(x) {
    computeTrajStats(traj_df = traj_df,
                     this_compartment = x) %>% 
      mutate(what = x)
  })

toc(quiet = FALSE, func.toc = printTicTocTime)

# Save --------------------------------------------------------------------

saveRDS(
  traj_stats,
  file = pomp_traj_stats_file
)

