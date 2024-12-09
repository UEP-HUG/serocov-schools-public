# This script is to determine the natural history parameters for variants
# The analysis consists in two steps:
#   - Fit Erlang distributions to incubation periods
#   - Find optimal value of k_P, k_E, k_I, and mu_I to match generation times

# Preamble ----------------------------------------------------------------
library(tidyverse)

source("analysis/utils.R")

do_plots <- F

# Load review -------------------------------------------------------------

review_data <- read_csv("data/02_review_natural_history/review_sarscov2_params.csv") %>%
  mutate(row = row_number()) %>% 
  group_by(row) %>% 
  group_modify(function(x, y) {
    if (is.na(x$estimate_scale[1])) {
      x$estimate_scale <- (x$estimate_sd[1])^2/x$estimate_mean[1]
      x$estimate_shape <- x$estimate_mean[1]/x$estimate_scale
    }
    
    if (is.na(x$estimate_mean[1])) {
      x$estimate_mean <- x$estimate_shape * x$estimate_scale
      x$estimate_sd <- sqrt(x$estimate_shape) * x$estimate_scale
    }
    x
  }) %>% 
  ungroup() %>% 
  select(-row)

saveRDS(review_data, "generated_data/nathist_review_param.rds")

if (do_plots) {
  # Plot gamma densities
  param_samples <- review_data %>% 
    filter(!is.na(estimate_shape)) %>% 
    group_by(reference, variant, parameter) %>% 
    group_modify(function(x,y) {
      tibble(rnd = rgamma(1e4,
                          shape = x$estimate_shape[1], 
                          scale = x$estimate_scale[1]))
    }) 
  
  param_samples %>% 
    ggplot(aes(x = rnd, color = variant)) +
    geom_density(aes(lty = reference, color = variant)) +
    facet_grid(reference ~ parameter) +
    theme_bw() +
    coord_cartesian(xlim = c(0, 15)) +
    scale_linetype_manual(values = rep(1, 6))
}


# Step 1: Erlang distributions of incubation period -----------------------

# Function to optimize scale parameter for given integer value of shape
optimizeErlangFit <- function(sample) {
  map_df(
    1:10, 
    function(x) {
      res <- optim(.1,
                   # Objective function to minimize: log-likelihood
                   function(param, shape, sample) {
                     scale <- param[1]
                     # Compute gamma density
                     ll <- sum(dgamma(sample, shape = shape, scale = scale, log = T))
                     -ll
                   }, 
                   shape = x, 
                   sample = sample, 
                   method = "Brent",
                   lower = 0, upper = 100)
      tibble(
        ll = res$value,
        shape = x,
        scale = res$par
      )
    })
}

# Fit Erlangs to all natural history parameters (only incubation periods used afterwards)
erlang_fits <- review_data %>% 
  filter(!is.na(estimate_shape)) %>% 
  group_by(reference, variant, parameter, 
           estimate_mean, estimate_sd, 
           estimate_shape, estimate_scale) %>% 
  group_modify(function(x, y) {
    
    sample <- rgamma(1e4,
                     shape = y$estimate_shape[1], 
                     scale = y$estimate_scale[1])
    
    erlang_res <- optimizeErlangFit(sample) %>%
      slice_min(ll)
    
    tibble(erlang_shape = erlang_res$shape[1],
           erlang_scale = erlang_res$scale[1])
  }) 


# Save
erlang_fits %>% 
  select(reference, variant, parameter, contains("estimate"), contains("erlang"), -contains("rnd")) %>% 
  distinct() %>% 
  ungroup() %>% 
  filter(parameter != "serial interval") %>% 
  saveRDS("generated_data/review_nathist_erlang.rds")


# Find optimal E/P/I parameters -------------------------------------------
# Then find the optimal values of k_E, k_P, k_I and mu_I to match the mean and variance of inferred generation times
# We use the following data:
# Incubation period: Galmiche 2023
# Generation time:
#   - Alpha and Delta: Hart2022
#   - Omicron: Modification of Delta of Hart2022
# 
variants <- c("alpha", "omicron", "delta")

variant_target_params <- map_df(variants, function(v) {
  # Extract the target incubation period and generation time parameters from literature
  inc <- erlang_fits %>% 
    filter(variant == v, 
           parameter == "incubation period")
  
  # Optimal Erlang shape parameter
  k_inc <- inc$erlang_shape[1]
  # Inverse optimal Erlang mean
  gamma <- 1/inc$estimate_mean[1]
  
  if (v != "omicron") {
    gen <- erlang_fits %>% 
      filter(variant == v, 
             parameter == "generation time",
             str_detect(reference, "Hart"))
  } else {
    gen <- erlang_fits %>% 
      filter(variant == "delta", 
             parameter == "generation time",
             str_detect(reference, "Hart")) %>% 
      mutate(estimate_mean = estimate_mean * .9)
  }
  
  e_gen <- gen$estimate_mean[1]
  v_gen <- gen$estimate_sd[1]^2
  
  # Add information about relative strength of transmission
  alpha_P_I <- case_when(
    v == "alpha" ~ 2.8,  # Hart 2022
    v == "delta" ~ 3.9,  # Hart 2022
    v == "omicron" ~ 3.9  # Assumed to be the same as Delta
  )
  
  tibble(
    variant = v,
    k_inc = k_inc, 
    gamma = gamma,
    e_gen = e_gen,
    v_gen = v_gen,
    alpha_P_I = alpha_P_I
  )
})

all_params <- map_df(
  variants, 
  function(v) {
    
    dat <- variant_target_params %>% 
      filter(variant == v)
    
    all_res <- optimizeErlangs(k_inc = dat$k_inc,
                               gamma = dat$gamma,
                               alpha_P_I = dat$alpha_P_I,
                               e_gen = dat$e_gen,
                               v_gen = dat$v_gen) %>% 
      mutate(e_gen = dat$e_gen,
             v_gen = dat$v_gen,
             variant = v) %>% 
      slice_min(err)
  })

saveRDS(all_params, "generated_data/variant_nathis_erlang_fits.rds")

# Save parameters for instantiating model
nathist_params <- all_params %>% 
  select(variant, k_P, k_I, mu) %>% 
  inner_join(variant_target_params %>% 
               select(variant, k_inc, gamma, alpha_P_I)) %>% 
  mutate(scale_inc = 1/(k_inc*gamma),
         scale_I = 1/mu)

saveRDS(nathist_params, "generated_data/variant_nathis_erlang_param.rds")
