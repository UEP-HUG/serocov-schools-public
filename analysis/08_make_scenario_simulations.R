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

Sys.setenv("REDO_DATA" = FALSE)


# Make simulations ---------------------------------------------------------
# Read in default config
config_dir <- "analysis/configs/"
dconfig <- yaml::read_yaml(str_c(config_dir, "default_config.yml")) 

# Simulation length in number of days
sim_length <- 60
 
# Define scenarios
param_scenarios <- list(
  alpha_like = list(
    lambda = .045,
    lambdaA = 1.1,
    lambdaB = .09,
    date_start = dateStartAlpha(),
    date_end = dateStartAlpha() + sim_length
  ),
  delta_like = list(
    lambda = .28,
    lambdaA = 1,
    lambdaB = .06,
    date_start = dateStartDelta(),
    date_end = dateStartDelta() + sim_length
  ),
  omicron_like = list(
    lambda = 9.5,
    lambdaA = 1.3,
    lambdaB = .4,
    date_start = dateStartOmicron(),
    date_end = dateStartOmicron() + sim_length
  )
)

# Class sizes
group_size_stats <- getSchoolCounts() %>% 
  mutate(institution = case_when(school %in% c("SA", "SB") ~ "school",
                                 TRUE ~ "preschool")) %>% 
  # Remove admin
  filter(group != "G1") %>% 
  group_by(school, age_cat) %>% 
  mutate(n_group = n()) %>% 
  ungroup() %>% 
  group_by(institution, age_cat) %>% 
  summarise(mean_size = mean(n_tot),
            sd_size = sd(n_tot),
            mean_groups = mean(n_group),
            sd_groups = sd(n_group)) %>% 
  # To list for scenarios
  group_by(institution) %>% 
  group_map(function(x, y) {
    res <- list(
      institution = y$institution[1],
      n_groups = ceiling(x$mean_groups[x$age_cat == "child"]),
      pupils_per_group = ceiling(x$mean_size[x$age_cat == "child"]),
      teachers_per_group = ceiling(x$mean_size[x$age_cat == "adult"])
    )
    res
  }) %>% 
  magrittr::set_names(map_chr(., ~ .$institution))


# Define within and between group interventions
within_intervention <- seq(0, 1, by = .1)
between_intervention <- seq(0, 1, by = .1)

# Write configs
walk(names(param_scenarios), function(x) {
  walk(names(group_size_stats), function(y) {
    walk(within_intervention, function(w) {
      walk(between_intervention, function(b) {
        
        # Initialize to default
        config <- dconfig
        
        # Scenario name for parsing
        scenario_name <- str_glue("{y}_{x}_w{w}_b{b}")
        config$name <- scenario_name
        
        # Set parameters
        config$param_specs <- param_scenarios[[x]]
        
        # Set group sizes
        config$group_specs <- group_size_stats[[y]]
        config$group_specs$household_prob <- dconfig$group_specs$household_prob
        
        # Set interventions
        config$intervention_specs$within_group <- w
        config$intervention_specs$between_group <- b
        
        yaml::write_yaml(config, 
                         file = str_c(config_dir, scenario_name, ".yml"))
      })
    })
    
  })
})
