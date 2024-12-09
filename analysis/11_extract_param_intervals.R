# This script extracts the parameter intervals for each run subset

# Preamble ----------------------------------------------------------------

library(tidyverse)
library(future)

source("analysis/utils.R")
source("analysis/pomp_utils.R")

plan(multisession, workers = 5)


# Functions ---------------------------------------------------------------


parseProfileParam <- function(x) {
  x %>% 
    str_extract(
      str_c(rep(c("lambda","lambdaA", "lambdaB"), times = 3),
            rep(c("_a", "_d", "_o"), each = 3),
            collapse = "|")
    )  %>% 
    str_c("1")
}

getProfiles <- function(hash,
                        redo = TRUE) {
  
  out_file <- str_glue("generated_data/profiles_{hash}.rds")
  
  if (!file.exists(out_file) | redo) {
    
    # Directories with profiling results
    profile_dirs <- dir("generated_data/", pattern = "pomp_profile", full.names = T) %>% 
      str_subset(hash) %>% 
      str_subset("9_by_5") 
    
    # Loop over files and load likelihood values
    profiles <- furrr::future_map_dfr(
      profile_dirs, 
      function(x) {
        
        profile_param <-parseProfileParam(x)
        
        map_df(dir(x, full.names = T),
               function(y) {
                 out <- readRDS(y)
                 
                 if ("res_df" %in% names(out)) {
                   out$res_df %>% 
                     select(ll, ll_se, profile_param) %>%
                     rename(c("value" = profile_param)) %>% 
                     mutate(id = str_extract(y, "[0-9]+(?=\\.rds)") %>% as.numeric(),
                            profile_param = profile_param) 
                 } else {
                   tibble(ll = NA,
                          ll_se = NA,
                          id = str_extract(y, "[0-9]+(?=\\.rds)") %>% as.numeric(),
                          profile_param = profile_param)
                 }
               })
        
      })
    
    
    saveRDS(profiles, file = out_file)
    
  } else {
    profiles <- readRDS(out_file)
  }
  
  profiles
  
}


parseRunFileInfo <- function(run_file) {
  
  hash <- str_extract(run_file, "(?<=run_data_)(.)*") %>% 
    str_sub(1, 4)
  
  # What setting is it
  setting <- case_when(!str_detect(run_file, "sa") ~ "preschools",
                       !str_detect(run_file, "sc") ~ "schools",
                       TRUE ~ "all")
  
  # Use of genetic likelihood
  gen_like <- str_extract(run_file, "true|false") %>% 
    as.logical() %>% 
    ifelse("with_gen", "no_gen")
  
  tibble(hash = hash,
         setting = setting,
         gen_like = gen_like)
}

plotRawProfiles <- function(profiles,
                            sub_profiles) {
  profiles %>%
    left_join(sub_profiles %>% mutate(in_set = TRUE)) %>%
    replace_na(list(in_set = FALSE)) %>%
    ggplot(aes(x = value, y = ll, color = in_set)) +
    geom_point(alpha = .5) +
    geom_errorbar(aes(ymin = ll - ll_se * 1.96,
                      ymax = ll + ll_se * 1.96),
                  width = 0, alpha = .3) +
    facet_wrap(~profile_param, scales = "free_x") +
    theme_bw()
}


applyProfileFilter <- function(profiles,
                               hash) {
  
  # Get of profile filters
  filters <- getProfileFilter()
  
  if (hash %in% names(filters)) {
    cat("-- Applying profile filter for hash", hash, "\n")
    res <- filters[[hash]](profiles)
  } else {
    cat("-- No filters found for hash", hash, ". Returning raw profiles. \n")
    res <- profiles
  }
  
  res %>% 
    filter(!is.na(ll))
}





computeMCAP <- function(hashes,
                        mcap_span = .75,
                        redo_profiles = FALSE,
                        do_filter = TRUE) {
  
  purrr::map(hashes, function(hash) {
    # A. Get the profiles
    sub_profiles <- getSubProfiles(hash = hash,
                                   redo = redo_profiles,
                                   do_filter = do_filter)
    
    # B. MCAP intervals
    mcap_fits <- sub_profiles %>% 
      group_by(profile_param) %>% 
      group_map(function(x, y) {
        mcap2(logLik = x$ll, 
              parameter = x$value,
              span = mcap_span
        ) %>% 
          append(list(param = y$profile_param[1]))
      })
    
    # C. Extract MCAP smooth
    mcap_smooth <- map_df(mcap_fits, function(x) {
      x$fit %>% 
        as_tibble() %>% 
        mutate(profile_param = x$param %>% str_remove("1")) %>% 
        mutate(max_ll = x$max_ll,
               delta = x$delta) %>% 
        mutate(thresh_ll = max_ll - delta)
    }) 
    
    list(mcap_fits = mcap_fits,
         mcap_smooth = mcap_smooth)
  }) %>% 
    magrittr::set_names(hashes)
}

extractParams <- function(hashes,
                          mcap_span = .75,
                          redo_profiles = FALSE,
                          do_filter = TRUE) {
  
  map_df(hashes, function(hash) {
    
    mcap_dat <- computeMCAP(hashes = hash,
                            mcap_span = mcap_span,
                            redo_profiles = redo_profiles,
                            do_filter = do_filter) %>% 
      .[[1]]
    
    # D. Extract param estimates
    param_estimates <- map_df(mcap_dat$mcap_fits, 
                              function(x) {
                                tibble(
                                  param = x$param,
                                  mle = x$mle,
                                  max_ll = x$max_ll,
                                  lo = x$ci[1],
                                  hi = x$ci[2]
                                )
                              }) %>% 
      mutate(variant = getVarDict()[str_extract(param, "(?<=_)[a,d,o]")],
             param = str_extract(param, "lambdaA|lambdaB|lambda")) %>% 
      bind_cols(parseRunFileInfo(getRunDataFile(hash)))
    
    
    saveRDS(param_estimates, str_glue("generated_data/param_estimates_{hash}.rds"))
    
    param_estimates
  })
}



getSubProfiles <- function(hash, 
                           redo_profiles = FALSE,
                           n_min_profiles = 3,
                           do_filter = TRUE) {
  # A. Get the profiles
  profiles <- getProfiles(hash = hash,
                          redo = redo_profiles)
  
  if(do_filter) {
    sub_profiles <- applyProfileFilter(profiles = profiles,
                                       hash = hash) 
  } else {
    cat("-- Skipping filter.\n")
    sub_profiles <- profiles
  }
  
  sub_profiles <- sub_profiles %>% 
    group_by(profile_param, value) %>%
    arrange(profile_param, value, desc(ll)) %>% 
    slice(1:4) %>% 
    ungroup()
  
  
  # Keep only parameters for which we have at least 4 datapoints
  keep_params <- sub_profiles %>% 
    count(profile_param, value) %>% 
    count(profile_param) %>% 
    filter(n > n_min_profiles)
  
  sub_profiles <- inner_join(sub_profiles,  
                             keep_params)
  
  sub_profiles
}


# Helper to debug profile
plotProfiles <- function(hash,
                         redo_profiles = FALSE,
                         mcap_span = .75,
                         do_filter = TRUE) {
  
  # A. Get the profiles
  sub_profiles <- getSubProfiles(hash = hash,
                                 redo = redo_profiles,
                                 do_filter = do_filter)
  # B. MCAP intervals
  mcap_dat <- computeMCAP(hashes = hash,
                          mcap_span = mcap_span,
                          redo_profiles = redo_profiles,
                          do_filter = do_filter) %>% 
    .[[1]]
  
  y_lim <- sub_profiles %>% 
    group_by(profile_param) %>% 
    summarise(min_ll = min(ll - 2.5 * ll_se),
              max_ll = max(ll + 2.5 * ll_se)) %>%
    ungroup() %>% 
    summarise(low = mean(min_ll),
              hi = max(max_ll)) %>% 
    unlist()
  
  
  sub_profiles %>% 
    mutate(profile_param = str_remove(profile_param, "1"),
           variant = getVarDict()[str_extract(profile_param, "(?<=_)[a,d,o]")]) %>% 
    arrange(profile_param, id) %>% 
    mutate(full_id = str_c(profile_param, id)) %>% 
    # Compute custom dodge
    group_by(profile_param) %>% 
    group_modify(function(x, y) {
      range <- diff(range(x$value))
      spacing <- range/50
      n_elem <- x %>% group_by(value) %>% slice(1) %>% nrow()
      n_values <- x %>% count(value) %>% pull(n) %>% max()
      values <- seq(1, n_values, by = 1)
      std_values <- (values-mean(values))/sd(values) * spacing
      x %>% 
        arrange(value, desc(ll)) %>% 
        mutate(spaced_value = value + std_values[(row_number()%%n_values + 1)])
    }) %>% 
    mutate(id= str_c(profile_param, id, sep = "-")) %>% 
    ggplot(aes(x = value, y = ll, color = variant, group = id)) +
    geom_errorbar(aes(ymin = ll - ll_se*1.96, ymax = ll + ll_se*1.96), 
                  width = 0, lwd = .4, alpha = .5) +
    geom_point(alpha = .5)  +
    geom_line(data = mcap_dat$mcap_smooth,
              inherit.aes = F,
              aes(x = parameter, y = smoothed, group = profile_param),
              color = "black",
              lty = 2) +
    geom_line(data = mcap_dat$mcap_smooth,
              inherit.aes = F,
              aes(x = parameter, y = quadratic, group = profile_param),
              color = "red") +
    geom_hline(data = mcap_dat$mcap_smooth %>% 
                 group_by(profile_param) %>% 
                 slice(1),
               aes(yintercept = thresh_ll), lty = 3, lwd = .6) +
    geom_vline(data = mcap_dat$mcap_smooth %>% 
                 filter(smoothed == max_ll),
               aes(xintercept = parameter), lty = 3, lwd = .6) +
    facet_wrap(~ profile_param, scales = "free_x") +
    theme_bw() +
    scale_color_manual(values = variantColors()) +
    coord_cartesian(ylim = y_lim) +
    ggtitle(parseRunFileInfo(getRunDataFile(hash)) %>% unlist() %>% str_c(collapse = " "))
  
}



getProfileFilter <- function() {
  
  filters <- list(
    "4e53" = function(df) {
      df %>%
        filter(!(str_detect(profile_param, "lambda_a") & (value > .078))) %>%
        filter(!(str_detect(profile_param, "lambda_d") & (value < .25))) %>%
        filter(!(str_detect(profile_param, "lambda_o") & (value > 19))) %>%
        filter(!(str_detect(profile_param, "lambdaB_o") & (value < .2))) %>%
        filter(!(str_detect(profile_param, "lambdaA_a") & (value < .7 | value > 2))) %>%
        filter(!(str_detect(profile_param, "lambdaA_d") & (value < .49))) %>%
        filter(!(str_detect(profile_param, "lambdaA_o") & (value < .5)))
    },
    "75da" = function(df) {
      df %>%
        filter(!(str_detect(profile_param, "lambda_a") & (value < .05))) %>%
        filter(!(str_detect(profile_param, "lambda_d") & (value < .5))) %>%
        filter(!(str_detect(profile_param, "lambdaA_o") & (value < .5 | (value > 1.25 & value < 1.5) | value > 2.25))) %>%
        filter(!(str_detect(profile_param, "lambdaB_o") & (value < .25)))
    },
    "11fe" = function(df) {
      df  %>%
        filter(!(str_detect(profile_param, "lambda_o") & (value > 15))) %>%
        filter(!(str_detect(profile_param, "lambda_d") & (value < .25))) %>%
        filter(!(str_detect(profile_param, "lambdaA_o") & (value < .8 | value > 2))) %>%
        filter(!(str_detect(profile_param, "lambdaB_o") & value < .11))
    } ,
    "5ea2" = function(df) {
      df  %>%
        filter(!(str_detect(profile_param, "lambda_a") & (value > .06))) %>%
        filter(!(str_detect(profile_param, "lambda_d") & (value < .05 | value > .8))) %>%
        filter(!(str_detect(profile_param, "lambdaA_o") & (value < .75 | value > 2))) %>%
        filter(!(str_detect(profile_param, "lambdaB_o") & value > .4))
    } ,
    "65cc" = function(df) {
      df  %>%
        filter(!(str_detect(profile_param, "lambda_a") & (value < .055))) %>%
        filter(!(str_detect(profile_param, "lambda_d") & (value > .55))) %>%
        filter(!(str_detect(profile_param, "lambdaA_a") & (value > 1.7 | value < .5)))  %>%
        filter(!(str_detect(profile_param, "lambdaA_d") & (value < 0.47 | value > 1.5)))  %>%
        filter(!(str_detect(profile_param, "lambdaA_o") & (value < .5)))  %>%
        filter(!(str_detect(profile_param, "lambdaB_d") & (value < .1))) %>%
        filter(!(str_detect(profile_param, "lambdaB_o") & (value < .2 | value > .5)))
    } ,
    "6823" = function(df) {
      df %>%
        filter(!(str_detect(profile_param, "lambda_a") & (value < .04))) %>%
        filter(!(str_detect(profile_param, "lambda_d") & (value < .1))) %>%
        filter(!(str_detect(profile_param, "lambdaA_a") & value > 1.5)) %>%
        filter(!(str_detect(profile_param, "lambdaA_o") & value < .5)) %>%
        filter(!(str_detect(profile_param, "lambdaB_a") & (value < 0.025 | value > .2))) %>%
        filter(!(str_detect(profile_param, "lambdaB_o") & value < .1))
    }
  )
  filters
}

# Run hashes
# [1] "4e53"  "75da"  "11fe" "5ea2"  "65cc"  "6823" 

plotProfiles(hash = "65cc",
             mcap_span = .85,
             redo_profiles = FALSE,
             do_filter = TRUE)


makeProfilePlots <- function(hash,
                             redo_profiles = TRUE,
                             do_filter = TRUE,
                             mcap_span = .75) {
  
  # A. Get the profiles
  profiles <- getProfiles(hash = hash,
                          redo = redo_profiles)
  
  # sub_profiles <- profiles
  sub_profiles <- getSubProfiles(hash = hash,
                                 redo = redo_profiles,
                                 do_filter = do_filter)
  
  
  # A. Raw profile ----
  cat("--- Raw profile figure", hash, " \n")
  
  p_raw <- plotRawProfiles(profiles, sub_profiles) +
    ggtitle(parseRunFileInfo(getRunDataFile(hash)) %>% unlist() %>% str_c(collapse = " "))
  
  ggsave(p_raw, 
         filename = str_glue("figures/raw_profiles_{hash}.png"),
         width = 10,
         height = 7)
  
  # B. Full profiles ----
  
  # Get MCAP data
  mcap_dat <- computeMCAP(hashes = hash,
                          mcap_span = mcap_span,
                          redo_profiles = redo_profiles,
                          do_filter = do_filter) %>% 
    .[[1]]
  
  
  # Extract param estimates
  param_estimates <- extractParams(hashes = run_hashes,
                                   mcap_span = mcap_span,
                                   redo_profiles = redo_profiles,
                                   do_filter = do_filter)
  
  # Compute y axis limits
  y_lim <- sub_profiles %>% 
    group_by(profile_param) %>% 
    summarise(min_ll = min(ll - 2.5 * ll_se),
              max_ll = max(ll + 2.5 * ll_se)) %>%
    ungroup() %>% 
    summarise(low = mean(min_ll),
              hi = max(max_ll)) %>% 
    unlist()
  
  
  
  # E. Profile figure
  cat("--- Profile figure", hash, "\n")
  
  p_profiles <- sub_profiles %>% 
    mutate(profile_param = str_remove(profile_param, "1"),
           variant = getVarDict()[str_extract(profile_param, "(?<=_)[a,d,o]")]) %>% 
    arrange(profile_param, id) %>% 
    mutate(full_id = str_c(profile_param, id)) %>% 
    # Compute custom dodge
    group_by(profile_param) %>% 
    group_modify(function(x, y) {
      range <- diff(range(x$value))
      spacing <- range/50
      n_elem <- x %>% group_by(value) %>% slice(1) %>% nrow()
      n_values <- x %>% count(value) %>% pull(n) %>% max()
      values <- seq(1, n_values, by = 1)
      std_values <- (values-mean(values))/sd(values) * spacing
      x %>% 
        arrange(value, desc(ll)) %>% 
        mutate(spaced_value = value + std_values[(row_number()%%n_values + 1)])
    }) %>% 
    mutate(id= str_c(profile_param, id, sep = "-")) %>% 
    ggplot(aes(x = value, y = ll, color = variant, group = id)) +
    geom_errorbar(aes(ymin = ll - ll_se*1.96, ymax = ll + ll_se*1.96), 
                  width = 0, lwd = .4, alpha = .5) +
    geom_point(alpha = .5)  +
    geom_line(data = mcap_dat$mcap_smooth,
              inherit.aes = F,
              aes(x = parameter, y = smoothed, group = profile_param),
              color = "black",
              lty = 2) +
    geom_line(data = mcap_dat$mcap_smooth,
              inherit.aes = F,
              aes(x = parameter, y = quadratic, group = profile_param),
              color = "red") +
    geom_hline(data = mcap_dat$mcap_smooth %>% 
                 group_by(profile_param) %>% 
                 slice(1),
               aes(yintercept = thresh_ll), lty = 3, lwd = .6) +
    facet_wrap(~ profile_param, scales = "free_x") +
    theme_bw() +
    scale_color_manual(values = variantColors()) +
    coord_cartesian(ylim = y_lim) +
    ggtitle(parseRunFileInfo(getRunDataFile(hash)) %>% unlist() %>% str_c(collapse = " "))
  
  
  ggsave(p_profiles, 
         filename = str_glue("figures/profiles_{hash}.png"),
         width = 10,
         height = 7)
  
}


# Define run sets ---------------------------------------------------------

run_hashes <- c("4e53", "75da", "11fe", "5ea2", "65cc", "6823")

# Parse settings
run_info <- map_df(run_hashes, function(h) {
  # Get the run data file
  run_file <- getRunDataFile(hash = h)
  # Parse info
  parseRunFileInfo(run_file)
})


# Load the files ----------------------------------------------------------


param_estimates <- extractParams(hashes = run_hashes,
                                 redo_profiles = FALSE,
                                 do_filter = TRUE,
                                 mcap_span = .85)


# walk(run_hashes, ~ makeProfilePlots(.))

saveRDS(param_estimates, file = "generated_data/all_param_estimates.rds")

# Compile param estimates -------------------------------------------------

# Scraps ------------------------------------------------------------------
# 
# 
fits <- dir("generated_data", pattern = "profile", full.names = TRUE) %>%
  str_subset("4e53") %>%
  str_subset("lambdaA_d") %>%
  # str_subset("h1.25") %>%
  dir(full.names = TRUE) %>%
  map_df(function(x) {
    res <- readRDS(x)
    res$res_df %>%
      mutate(file = x)
  })  %>%
  mutate(run = str_extract(file, "[0-9]+(?=\\.rds)") %>% as.numeric())

fits %>%
  select(run, ll, ll_se, contains("lambda")) %>%
  pivot_longer(cols = contains("lambda")) %>%
  mutate(param = str_remove(name, "[0-9]+")) %>%
  group_by(run, param) %>%
  summarise(value = mean(value)) %>%
  filter(param != "lambdaA_d") %>%
  bind_rows(
    fits %>%
      select(run, ll, ll_se, contains("C_")) %>%
      pivot_longer(cols = contains("C_0")) %>%
      group_by(run) %>%
      summarise(value = sum(value>.5)) %>%
      mutate(param = "n_C")
  ) %>%
  inner_join(
    fits %>%
      distinct(run, ll, ll_se, lambdaA_d = lambdaA_d1)
  ) %>%
  filter(!(str_detect(param, "wt"))) %>%
  ggplot(aes(x = lambdaA_d, y = value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~param, scales = "free_y") +
  theme_bw()

