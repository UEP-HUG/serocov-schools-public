# This script aims at modeling the seroprevalence by outbreak to prodived
# a probe for testing the pomp model outputs.


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(posterior)
library(bayesplot)
library(cmdstanr)
library(here)

source("analysis/utils.R")
source("analysis/stan_utils.R")

# Whether to redo the study data processing (school groups and household info)
Sys.setenv("REDO_DATA" = "FALSE")

# Parameters
do_plots <- F

# Glossary of variables and column names
# - term: school term
# - group: school group
# - test: PCR/antigen anti-SARS-CoV-2 tests
# - obs/observation: serology result
# - int/observation interval: period between two (potentially fake) sampling points 


# A. Load contextual dat -----------------------------------------------------

# Composition of school groups
school_groups <- getDataGroups()

# Composition of households
households <- getDataHouseholds() %>% 
  mutate(age = case_when(age == 2021 ~ 40,
                         TRUE ~ age))

# Vaccination data
vaccination_data <- getDataVaccination() %>% 
  mutate(earliest_vacc_date = coalesce(d1_date, d2_date, d3_date)) %>% 
  select(sugar_id, earliest_vacc_date)

# Long format of school outbreak metadata. For definitions of "visites" see 
# data description notebook.
outbreak_metadata <- getOutbreakMetadata() %>% 
  select(outbreak_id, school, group, visite_1, visite_3) %>% 
  pivot_longer(cols = contains("visite"),
               names_to = "event",
               values_to = "date") %>% 
  mutate(event = str_c("v", str_extract(event, "[1-3]"))) %>% 
  addTerms(date_col = "date")

# Get the dates of other study sampling dates besides outbreaks (already in long format)
timeline_metadata <- getTimelineMetadata()

# This is the available serology and test data at baseline and endline
timeline_data <- getDataStartEnd()

# Combine the metadata. We need to expand the timeline metadata which is only by school.
all_metadata <- bind_rows(
  outbreak_metadata,
  timeline_metadata %>% 
    group_by(school, event, term) %>% 
    summarise(date = mean(c(TL, TR))) %>% 
    # Expand to include all groups
    left_join(
      school_groups %>% 
        select(school, group, term)
    ) %>% 
    # Add fake outbreak id to identify these events
    mutate(outbreak_id = str_c(event, group, term, sep = "-")) %>% 
    filter(!is.na(date))
)

# B. Load outbreak data ------------------------------------------------------

# Merged PCR/antigen data from study and ARGOS [origin: 01_parse_questionnaire.R]
merged_test_data <- readRDS("generated_data/merged_test_argos_data.rds")

# Pull all the serology/test data by outbreak
outbreak_data <- getDataOutbreak(ob_ids = getAllOutbreaks()) %>% 
  mutate(
    # There are some outbreaks for which visite_2 is NA
    visite_2 = coalesce(visite_2, visite_1 + 7),
    # New columns by visite for whether the data is in visite or not
    across(
      contains("visite"), 
      function(x) {
        date >= (x  - 2) & date <= (x + 2)
      }, 
      .names = "v{str_extract(.col, '[1-3]')}")
  )

# Unique sugar ids in data
u_sugar_ids <- outbreak_data %>%
  distinct(school, sugar_id, term)

# C. Raw longitudinal serology data ---------------------------------------

# Combine outbreak and baseline/endline data
raw_sero_long <- bind_rows(
  outbreak_data %>% 
    select(outbreak_id, sugar_id, age_cat, school, 
           group, term, date, v1, v3, result, what) %>% 
    mutate(
      event = case_when(v1 ~ "v1", 
                        v3 ~ "v3", 
                        T ~ NA_character_)
    ) %>% 
    select(-v1, -v3),
  timeline_data %>% 
    select(-TL, TR) %>% 
    mutate(outbreak_id = str_c(event, group, term, sep = "-"))
) %>%     
  # !! Keep only serologies here
  filter(what == "sero") %>%
  # !! Remove data outside of either outbreak visits or baseline/endlines
  filter(!is.na(event)) %>%
  # !! Remove admin
  filter(group != "G1") %>% 
  # !! Remove Magnolias
  filter(school != "SE")


# D. Longitudinal data processing --------------------------------------------
# This first step joins the test data to determine whether participants had
# prior positive test results.

sero_long_s1 <- raw_sero_long %>% 
  left_join(
    merged_test_data %>% 
      select(sugar_id, test_date, test_result)) %>% 
  # Keep only test data prior to serology
  filter(test_date < date | is.na(test_date)) %>% 
  group_by(outbreak_id, sugar_id, age_cat, school, group, term, date, result, event) %>% 
  # Determine whether past test data is available and when. Differentiate between
  # the latest and earliest positive test data for modeling time-varying 
  # serology sensitivity.
  summarise(
    prev_pos = any(test_result == "positive", na.rm = T),
    latest_pos_date = max(test_date[test_result == "positive"], na.rm = T),
    earliest_pos_date = min(test_date[test_result == "positive"], na.rm = T),
    has_prev_data = any(!is.na(test_result))
  ) %>% 
  ungroup() %>% 
  left_join(raw_sero_long, .) %>% 
  arrange(sugar_id, date) %>% 
  # !! overwrite prev_pos when date not known
  mutate(prev_pos = case_when(
    is.infinite(latest_pos_date) ~ FALSE,
    T ~ prev_pos
  )) 


# E. Make fake sampling dates ------------------------------------------------

# Add missing outbreak v1 or v3 dates for participants for which we have serologies 
# on either date but not both.
missing_outbreak_dates <- sero_long_s1 %>%
  # Focus on outbreaks as other events have hyphens in their name
  filter(!str_detect(outbreak_id, "-")) %>%
  group_by(sugar_id, outbreak_id, school, group) %>% 
  summarise(
    n = n(),
    event = case_when(
      "v1" %in% event ~ "v3",
      T ~ "v1")
  ) %>% 
  # Each participant should have 2 data points (visite 1 and 3)
  filter(n < 2) %>% 
  inner_join(outbreak_metadata) %>% 
  ungroup()

# Missing participants from outbreaks for which we have data at other time
# points but none in specific outbreak.
missing_outbreak_participants <- school_groups %>%
  distinct(sugar_id, school, group, term) %>% 
  # Make full list of participant/group/outbreaks
  right_join(outbreak_metadata) %>% 
  # Add flag whether we alread have data for them
  left_join(
    sero_long_s1 %>%
      distinct(sugar_id, school, outbreak_id, group, term) %>% 
      mutate(in_data = T)
  ) %>% 
  filter(is.na(in_data)) %>% 
  select(-in_data)


# Missing baseline/endline data
missing_data_other_events <- timeline_metadata %>% 
  filter(!is.na(TL)) %>% 
  group_by(school, event, term) %>% 
  group_modify(function(x, y) {
    
    # Available data for this event
    avail_dat <- sero_long_s1 %>% 
      inner_join(y)
    
    # Expand to all groups and match on this school/event
    inner_join(
      school_groups %>% 
        select(sugar_id, school, group, term),
      y
    ) %>%
      bind_cols(x) %>%
      rename(date = TL) %>%
      select(-TR) %>% 
      # Keep only fake data for participants that are in school but for which
      # we do not have available event data
      filter(
        sugar_id %in% sero_long_s1$sugar_id,
        !(sugar_id %in% avail_dat$sugar_id)
      ) %>% 
      select(-all_of(colnames(y)))
  }) %>% 
  ungroup() %>% 
  # Create unique id following same pattern as timeline data
  mutate(outbreak_id = str_c(event, group, term, sep = "-"))


# We need to add fake baselines in 20/21 for all participants to model infection during 20/21.
# First determine which participant do not have a 20/21 baseline
missing_baseline_20_21 <- sero_long_s1 %>% 
  bind_rows(missing_data_other_events) %>%
  distinct(sugar_id, school, outbreak_id) %>%
  group_by(sugar_id) %>%
  mutate(has_baseline_20_21 = any(str_detect(outbreak_id, "baseline") & str_detect(outbreak_id, "20/21"))) %>%
  filter(!has_baseline_20_21) %>% 
  distinct(sugar_id, school)

# Get the metadata for the 20/21 baseline sampling dates. We need to add fake data
# for La Decouverte because it was not sampled in 20/21.
baseline_metadata <- timeline_metadata %>% 
  filter(event == "baseline", term == "20/21") %>% 
  # Make fake baseline date for La Decouverte in 2021
  bind_rows(
    tribble(
      ~school, ~event, ~TL, ~TR, ~term,
      "SB", "baseline", as.Date("2021-03-17"), as.Date("2021-03-17"), "20/21"
    )
  )

# Make explicit missing data for baselines 20/21
missing_data_baseline_20_21 <- missing_baseline_20_21 %>% 
  inner_join(baseline_metadata) %>% 
  rename(date = TL) %>%
  select(-TR) %>% 
  mutate(outbreak_id = str_c(event, school, term, sep = "-"))


# Fake data to model changes of SARS-CoV-2 variants dates for everyone
variant_changes <- map_df(
  c("delta", "omicron"), 
  function(x) {
    expand.grid(
      sugar_id = unique(sero_long_s1$sugar_id),
      term = getAllTerms()
    ) %>% 
      as_tibble() %>% 
      inner_join(
        getVariantDates(x) %>% 
          select(variant, date = TL) %>% 
          addTerms()
      ) %>% 
      mutate(
        event = "variant_change",
        outbreak_id = str_c(event, term, variant, sep = "-")
      ) %>% 
      select(-variant)
  })


# Fake data for vacations for everyone
vacations <- sero_long_s1 %>% 
  distinct(sugar_id, school) %>% 
  inner_join(
    # This is vacation dates for all participants
    getVacations() %>% 
      addTerms(date_col = "TL") %>% 
      inner_join(
        expand.grid(name = getVacations()$name,
                    school = getAllSchools())
      )
  ) %>%
  mutate(name = str_c("vacations_", name)) %>% 
  rename(event = name) %>% 
  pivot_longer(
    cols = c("TL", "TR"),
    names_to = "what",
    values_to = "date"
  ) %>% 
  mutate(
    event = case_when(
      what == "TL" ~ str_c(event, "_start"),
      T ~ str_c(event, "_end")
    )
  ) %>% 
  mutate(outbreak_id = str_c(event, school, term, sep = "-")) %>% 
  distinct() %>% 
  select(-what)


# Vaccination
vacc_events <- vaccination_data %>%
  filter(sugar_id %in% raw_sero_long$sugar_id) %>% 
  rename(date = earliest_vacc_date) %>% 
  addTerms() %>% 
  mutate(event = "vaccination",
         outbreak_id = str_c("vaccination", term, str_sub(sugar_id, 1, 4), sep = "-"))

# F. Combine serology and fake data --------------------------------------

# Combine all data
sero_long_s2 <- sero_long_s1 %>% 
  bind_rows(missing_outbreak_dates,
            missing_outbreak_participants,
            missing_data_other_events,
            missing_data_baseline_20_21,
            variant_changes,
            vacations,
            vacc_events) %>%
  addTerms() %>% 
  select(-age_cat, -n) %>% 
  # Add age category information
  inner_join(households %>% select(sugar_id, age_cat)) %>% 
  arrange(sugar_id, date)  %>% 
  replace_na(list(prev_pos = FALSE,
                  latest_pos_date = NA,
                  has_prev_data = FALSE))

# Fill the schools and groups for missing groups
sero_long_s2 <- sero_long_s2 %>% 
  select(-school, -group) %>% 
  left_join(school_groups %>% distinct(sugar_id, school, term))  %>% 
  left_join(school_groups %>% distinct(sugar_id, school, group, term))

# First schools
for (i in which(is.na(sero_long_s2$school))) {
  sero_long_s2$school[i] <- school_groups$school[school_groups$sugar_id == sero_long_s2$sugar_id[i]] %>% first()
}

# Then groups
for (i in which(is.na(sero_long_s2$group))) {
  sero_long_s2$group[i] <- school_groups$group[school_groups$sugar_id == sero_long_s2$sugar_id[i]] %>% first()
}


# Plot of sample collection timelines by participant to check if all participants
# have data for all events.
if (do_plots) {
  
  sero_long_s2 %>% 
    mutate(
      event_type = map_chr(outbreak_id, 
                           ~ makeEventType(., simplify = F, add_terms = T)),
      event_type = factor(event_type, levels = names(getEventTypeColors(with_terms = TRUE)))
    ) %>% 
    ggplot(aes(x = date, y = sugar_id)) +
    geom_line(aes(group = sugar_id), lwd = .2, alpha = .1) +
    geom_point(aes(color = event_type, alpha = is.na(result))) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("Longitudinal serological samples") +
    scale_alpha_manual(values = c(1, .2)) +
    scale_color_manual(values = getEventTypeColors(with_terms = TRUE))
}


# G. Fraction seropositive by outbreak -----------------------------------

# Compute of fraction seropositive by school and outbreak for display. 
all_event_data <- bind_rows(
  sero_long_s2 %>% 
    group_by(outbreak_id, school, group, event) %>% 
    computeProp(res_col = "result")  %>% 
    left_join(all_metadata) %>% 
    distinct(),
  # Add school-level data
  sero_long_s2 %>% 
    mutate(outbreak_id = case_when(
      str_detect(outbreak_id, "-") ~ str_c(event, school, term, sep = "-"),
      T ~ outbreak_id)
    ) %>% 
    group_by(outbreak_id, school, event) %>% 
    computeProp(res_col = "result") %>% 
    mutate(group = school)  %>% 
    left_join(
      all_metadata %>% 
        mutate(
          outbreak_id = case_when(
            str_detect(outbreak_id, "-") ~ str_c(event, school, term, sep = "-"),
            T ~ outbreak_id)
        ) %>% 
        distinct() %>% 
        group_by(outbreak_id, event, school) %>% 
        summarise(date = mean(date, na.rm = T))
    )
) %>% 
  addTerms() %>% 
  ungroup() %>% 
  mutate(full_id = str_c(outbreak_id, group, sep = "-"),
         is_school = group == school) %>% 
  filter(group != "G1") %>% 
  mutate(
    event_type = map_chr(outbreak_id, 
                         ~ makeEventType(., simplify = F, add_terms = F)),
    event_type = factor(event_type, levels = names(getEventTypeColors(with_terms = FALSE)))
  )

saveRDS(all_event_data, "generated_data/all_event_data.rds")

# Plot the fraction of seropositives. This is interesting to see if there
# are time-trends.
if (do_plots) {
  
  pd <- position_dodge(2)
  
  all_event_data %>% 
    ggplot(aes(x = date, y = prop,
               color = event_type, 
               alpha = is_school,
               lwd = is_school)) +
    geom_line(data = all_event_data %>% 
                filter(!is_school),
              inherit.aes = F,
              aes(x = date, y = prop,
                  alpha = is_school,
                  group = str_c(group, term, sep = "-")),
              lwd = .1,
              position = pd) +
    geom_point(data = all_event_data %>% 
                 filter(!is_school),
               aes(size = n, group = full_id), position = pd) +
    geom_point(data = all_event_data %>% 
                 filter(is_school),
               aes(size = n)) +
    geom_errorbar(data = all_event_data %>% 
                    filter(!is_school),
                  aes(ymin = lo, ymax = hi, group = full_id), width = 0, lwd = .35, position = pd) +
    geom_errorbar(data = all_event_data %>% 
                    filter(is_school),
                  aes(ymin = lo, ymax = hi), width = 0, lwd = .35) +
    theme_bw() +
    scale_alpha_manual(values = c(.3, 1)) +
    scale_linewidth_manual(values = c(.2, .5)) +
    scale_color_manual(values = getEventTypeColors(with_terms = F)) +
    facet_grid(school~.) +
    labs(y = "Proportion seropositive")
  
}


# H. Empirical proportion false-negative serologies ----------------------

# Compute proportion of false negatives
prop_false_negs <- sero_long_s2 %>% 
  filter(prev_pos) %>% 
  select(sugar_id, age_cat, result, date, earliest_pos_date, latest_pos_date) %>% 
  pivot_longer(
    cols = contains("pos_date"), 
    names_to = "what",
    values_to = "ref_date") %>% 
  mutate(variant = map_chr(ref_date, ~ getVariant(.))) %>% 
  mutate(delay = as.numeric(difftime(date, ref_date, units = "days")),
         delay_cat = cut(delay, c(0, 10, 30, 90, 360, Inf), include.lowest = T),
         false_neg = result == 0,
         emp_sens = result == 1) %>%
  group_by(age_cat, what, delay_cat) %>% 
  mutate(mean_delay = mean(delay)) %>% 
  group_by(what, age_cat, delay_cat, mean_delay) %>% 
  computeProp(res_col = "emp_sens") 

saveRDS(prop_false_negs, "generated_data/prop_false_negs.rds")

if (do_plots) {
  pd <- position_dodge(width = 1)
  
  prop_false_negs  %>% 
    ggplot(aes(x = mean_delay, y = prop)) +
    geom_point(aes(size = n), position = pd) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, position = pd) +
    theme_bw() +
    facet_grid(what ~ age_cat) +
    labs(x = "time from test+ to serology", 
         y = "proportion false negatives",
         title = "Proportion false positives by time to positive PCR/antigen test")
  
}

# if (do_plots) {
#   # Most participants are in 1 or 2 outbreak
#   sero_long_s2 %>% 
#     group_by(sugar_id) %>% 
#     summarise(n = length(unique(outbreak_id))) %>% 
#     ggplot(aes(x= n)) +
#     geom_histogram() +
#     theme_bw()
# }


# I. Analysis serology dataset -------------------------------------------
# We need to filter out uncesseray data.

# Make unique terms in which participants were sampled to filter out unecessary
# fake data.
u_terms <- sero_long_s2 %>% 
  group_by(sugar_id, term) %>% 
  summarise(any_sero_data = any(!is.na(result))) %>% 
  ungroup() %>% 
  filter(term == "20/21" | (term == "21/22" & any_sero_data)) %>% 
  select(-any_sero_data)

# Maximum date of serology or PCR/antigen test to remove fake data that come
# later.
max_obs_date <- raw_sero_long %>% 
  group_by(sugar_id) %>% 
  filter(date == max(date[!is.na(result)], na.rm = T)) %>% 
  select(sugar_id, max_obs_date = date) %>% 
  bind_rows(
    merged_test_data %>% 
      filter(sugar_id %in% raw_sero_long$sugar_id) %>% 
      group_by(sugar_id) %>% 
      filter(test_date == max(test_date[!is.na(test_result)], na.rm = T)) %>% 
      mutate(test_date = test_date + 60) %>% 
      select(sugar_id, max_obs_date = test_date)
  ) %>% 
  group_by(sugar_id) %>% 
  slice_max(max_obs_date, n = 1, with_ties = F) %>% 
  ungroup()

# Make the analysis dataset
sero_long <- sero_long_s2 %>% 
  select(sugar_id, outbreak_id, school, group, term, event, date, result, 
         prev_pos,earliest_pos_date, age_cat) %>% 
  ungroup() %>% 
  select(-any_of(c("n"))) %>% 
  # !! take distinct
  distinct() %>% 
  # !! take one group per participant
  group_by(sugar_id, date, event) %>% 
  slice(1) %>% 
  ungroup()  %>%
  # # !! remove adults in administrative group not involved in outbreaks
  filter(group != "G1") %>%
  # !! remove participants that are not in the households file
  filter(!is.na(age_cat)) %>%  
  # !! remove all fake data after last serology or pcr test
  left_join(max_obs_date) %>% 
  filter(date <= max_obs_date | is.na(max_obs_date)) %>%
  # !! remove all fake data for which we don't have any information in school term
  inner_join(u_terms) %>% 
  # # !! remove participants with only 1 time point
  add_count(sugar_id) %>%
  filter(n > 1) %>%
  arrange(sugar_id, date, outbreak_id) %>% 
  mutate(time_id = row_number()) %>% 
  # Factor reordering for furhter analysis
  mutate(age_cat = factor(age_cat) %>% forcats::fct_relevel("adult")) %>% 
  arrange(time_id)


# Plot of analysis dataset to see if data was dropped correctly
if (do_plots) {
  
  sero_long %>% 
    mutate(
      event_type = map_chr(outbreak_id, ~ makeEventType(., simplify = F, add_terms = TRUE)),
      event_type = factor(event_type, levels = names(getEventTypeColors(with_terms = T)))
    ) %>% 
    ggplot(aes(x = date, y = sugar_id)) +
    geom_line(aes(group = sugar_id), lwd = .2) +
    geom_point(aes(color = event_type, pch = prev_pos)) +
    theme_bw() +
    scale_color_manual(values = getEventTypeColors(with_terms = T)) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("Analysis longitudinal dataset")
}


# Add checks
# TODO

# Save analysis dataset
saveRDS(sero_long, "generated_data/sero_long.rds")

# Variables for stan model
M <- sero_long %>% distinct(sugar_id) %>% nrow()    # Number of participants
times <- as.numeric(sero_long$date)    # time in numeric days of dates of samples

N_obs <- sum(!is.na(sero_long$result))    # Number of observations
ind_sero_obs <- which(!is.na(sero_long$result))    # Indices of observations (not all rows have data as there is fake data)

# Observations vector, filled with 0 when NA (these won't be used)
y <- sero_long$result
y[is.na(y)] <- 0

# Vaccination status
vaccinated <- sero_long %>% 
  left_join(
    vaccination_data
  ) %>% 
  mutate(earliest_vacc_date = ifelse(is.na(earliest_vacc_date), 
                                     as.Date("2100-01-01"), 
                                     earliest_vacc_date),
         vaccinated = date >= (earliest_vacc_date + 14)) %>% 
  pull(vaccinated)

# J. Define observation intervals --------------------------------------------------------

intervals <- makeIntervals(df = sero_long, 
                           group_var = "sugar_id")

# Plot to check if intervals make sens
if (do_plots) {
  intervals %>% 
    ggplot(aes(y = sugar_id)) +
    geom_linerange(aes(xmin = TL, xmax = TR), lwd = .1) +
    geom_point(data = sero_long, aes(x = date)) +
    facet_grid(school ~ ., scales = "free") +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("Observation intervals by particiapant")
}


# Variables for stan model
L <- nrow(intervals)    # Number of observation intervals

# Start/end timepoints of intervals
int_starts <- intervals$tid_L    # timepoint id (row of sero_long) of left interval bound
int_ends <- intervals$tid_R      # timepoint id of right interval bound

# Map from serology to interval
map_obs_interval <- map_dbl(
  sero_long$time_id, 
  function(x) {
    ind <- which(intervals$tid_L == x)
    # If no match this is for the first observation of each participant
    if (length(ind) == 0) {
      return(0)
    }
    ind
  })

# Map to variants
map_int_variant_num <- map_chr(intervals$TL, ~ getVariant(.))
u_variants <- map_int_variant_num %>% sort() %>% unique()
map_int_variant <- map_dbl(map_int_variant_num, function(x){which(u_variants == x)})
V <- length(u_variants)    # total number of variants in dataset

# Sarts end ends of intervals for each participant. This assumes that the intervals
# dataframe is sorted by sugar_id and date
start_ends_df <- intervals %>% 
  group_by(sugar_id) %>% 
  summarise(start = min(int_id),
            end = max(int_id))

# Extract vectors
starts <- start_ends_df$start
ends <- start_ends_df$end

# Maximum number of intervals per participant
max_n_intervals <- intervals %>% 
  count(sugar_id) %>% 
  pull(n) %>% 
  max()


# Map from intervals to vaccination
map_int_vaccination <- str_detect(intervals$outbreak_id, "vaccination") %>% 
  as.numeric()

# K. Group mappings ------------------------------------------------------

# Combinations of outbreak/groups
outbreak_group_comb <- sero_long %>% 
  ungroup() %>% 
  distinct(outbreak_id, school, group) %>% 
  filter(!str_detect(outbreak_id, "-")) %>% 
  # Add visite_1 date for variant information
  inner_join(
    getOutbreakMetadata() %>% 
      select(outbreak_id, school, group, visite_1)
  ) %>%  
  arrange(as.numeric(outbreak_id), school, group) %>% 
  mutate(og_id = row_number())

# Map from outbreak to variant
map_outbreak_variant <- map_chr(outbreak_group_comb$visite_1, ~ getVariant(.)) %>% 
  map_dbl(function(x){which(u_variants == x)})

outbreak_group_all <- intervals %>% 
  left_join(outbreak_group_comb) %>% 
  replace_na(list(og_id = 0))


# L. Times to previous tests -------------------------------------------------

t_pos_test <- map_dbl(1:nrow(sero_long), function(x) {
  
  tp <- sero_long$earliest_pos_date[x]
  if (is.na(tp) | is.infinite(tp)) {
    return(0)
  } else {
    tp
  }
})

ind_prev_tests <- which(t_pos_test > 0)
N_tests <- sum(t_pos_test > 0)

dt_sens_pred <- seq(1, 500, by = 5)
saveRDS(dt_sens_pred, "generated_data/dt_sens_pred.rds")

prev_pos_df <-  sero_long %>% 
  ungroup() %>% 
  filter(prev_pos) %>% 
  mutate(dt = difftime(date, earliest_pos_date) %>% as.numeric()) 


# M. Covariates for sensitivity ------------------------------------------

X_sens <- model.matrix(~ age_cat, 
                       data = sero_long) %>% 
  # Remove first column which is intercept
  .[, -1, drop = F]

# Get unique combinations of covariates for generated quantities
X_sens_comb <- X_sens %>% unique()

# Expand to dataframe for plotting
X_sens_comb_df <- X_sens_comb %>% 
  as_tibble() %>% 
  rename(age_cat = age_catchild) %>% 
  mutate(across(everything(), 
                function(x) {
                  l <- levels(sero_long[[cur_column()]])
                  print(l)
                  l[x + 1]
                })) %>% 
  mutate(comb_row = row_number())

saveRDS(X_sens_comb_df, "generated_data/X_sens_comb_df.rds")


# Map of observations to baseline/outbreak visite dates
# Combinations of outbreak/groups including baseline
event_comb <- sero_long %>% 
  filter(!str_detect(event, "vaccination")) %>% 
  group_by(outbreak_id, event, school, group) %>% 
  summarise(date = mean(date),
            n_obs = sum(!is.na(result)),
            n_pos = sum(result, na.rm = T),
            n_tot = n()) %>% 
  arrange(outbreak_id, school, group) %>% 
  ungroup() %>% 
  mutate(event_id = row_number(),
         event_type = case_when(
           !str_detect(outbreak_id, "-") ~ str_c("outbreak-", outbreak_id),
           T ~ str_c(str_extract(outbreak_id, "baseline|last_visit|debut_school_year"),
                     str_extract(outbreak_id, "20/21|21/22"),
                     sep = "-")
         )) %>% 
  mutate(what = str_extract(event_type, "(.)*(?=-)")) %>% 
  addTerms()

saveRDS(event_comb, "generated_data/event_comb.rds")

# Map from observation to event
map_obs_event <- sero_long %>% 
  left_join(
    event_comb %>% 
      select(outbreak_id, event, group, event_id)
  ) %>% 
  pull(event_id)

map_obs_event[is.na(map_obs_event)] <- 0

# Covariates for intervals
X_int_df <- intervals %>% 
  select(variant, vacations) 

# Make model matrix for regression
X_int <- model.matrix(~ variant + vacations, data = X_int_df)
n_covar_int <- ncol(X_int)

# matrix of unique covariate combinations to define changes in community transmission (variants and vacations)
X_int_comb <- unique(X_int)   
Z_int_comb <- nrow(X_int_comb)    

# Define which elements of X_int_comb are periods in school (only variant effects)
map_outbreak_int_variant <- apply(X_int_comb, 1, function(x) {
  if (sum(x) == 1) {
    TRUE
  } else if (sum(x) == 2 && (x["variantdelta"] == 1 || x["variantomicron"] == 1)) {
    TRUE
  } else {
    FALSE
  }
}) %>% 
  which() %>% 
  as.integer()

# Map from interval to unique interval combination
map_int_int_comb <- apply(X_int, 1, function(x){
  which(
    map_lgl(1:nrow(X_int_comb), ~ identical(x, X_int_comb[.,]))
  )
})

# Setup for outbreaks as well
X_int_outbreak_group_df <- intervals %>% 
  distinct(variant, vacations, outbreak_id, group, school) %>% 
  filter(!str_detect(outbreak_id, "-"))  %>%  
  arrange(as.numeric(outbreak_id), school, group) %>% 
  left_join(outbreak_group_comb) %>% 
  mutate(variant = map_chr(visite_1, ~ getVariant(.))) %>% 
  distinct()

X_int_outbreak_group <- model.matrix(~ variant + vacations, data = X_int_outbreak_group_df)

map_outbreak_group_int_comb <- apply(X_int_outbreak_group, 1, function(x){
  which(
    map_lgl(1:nrow(X_int_comb), ~ identical(x, X_int_comb[.,]))
  )
})

# Checks
testthat::expect_identical(colnames(X_int), colnames(X_int_outbreak_group))


# N. PCR/antigen data -----------------------------------------------------
# We use PCR/antigen data to inform the community infections

# Get all test data from included participants
proc_test_data <- merged_test_data %>% 
  mutate(test_id = row_number()) %>% 
  inner_join(sero_long %>% 
               distinct(sugar_id, term, school, group)) %>%
  arrange(sugar_id, test_date) %>% 
  # !! keep only test data within min and max dates
  filter(test_date >= min(sero_long$date),
         test_date <= max(sero_long$date)) %>% 
  add_count(sugar_id) %>% 
  # add outbreak information
  left_join(getOutbreakMetadata() %>% 
              mutate(visite_2 = coalesce(visite_2, visite_1 + 5)) %>% 
              select(-variant)) %>%
  # determine if test was within an outbreak
  mutate(in_outbreak = case_when(
    is.na(visite_1) ~ FALSE,
    test_date >= (visite_1 - 2) & test_date <= (visite_2 + 5) ~ TRUE,
    T ~ FALSE)
  ) %>% 
  group_by(sugar_id, test_date) %>% 
  mutate(any_outbreak = any(in_outbreak)) %>%
  ungroup() %>% 
  # !! Keep only data in outbreak or outside of any outbreaks
  filter((any_outbreak & in_outbreak) | !any_outbreak) %>% 
  mutate(outbreak_id = case_when(in_outbreak ~ outbreak_id, 
                                 T ~ NA_character_),
         visite_1 = case_when(in_outbreak ~ visite_1, 
                              T ~ NA)) %>% 
  distinct(sugar_id, school, group, term, test_date, test_result, test_type,
           source, in_outbreak, outbreak_id, visite_1) %>% 
  ungroup() %>% 
  # add interval information
  left_join(intervals %>% select(-outbreak_id)) %>% 
  filter(test_date >= TL,
         test_date < TR) %>% 
  arrange(sugar_id, test_date) 

proc_test_data <- bind_rows(
  proc_test_data %>%
    filter(!is.na(outbreak_id)) %>% 
    # !! keep only positive tests
    group_by(outbreak_id, sugar_id) %>% 
    arrange(outbreak_id, sugar_id, desc(test_result)) %>% 
    slice(1) %>% 
    ungroup(),
  proc_test_data %>%
    filter(is.na(outbreak_id))
) %>% 
  arrange(test_date, sugar_id)


saveRDS(proc_test_data, "generated_data/test_data_long.rds")


if (do_plots) {
  proc_test_data %>% 
    arrange(term, group) %>% 
    ggplot(aes(x = test_date, y = sugar_id)) +
    geom_line(lwd = .1) +
    geom_point(aes(color = in_outbreak)) +
    theme_bw() +
    facet_grid(school ~ ., scales = "free") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("PCR/antigen data")
}


# Variables for stan model

# Map from intervals to outbreak id
map_int_outbreak_group <- outbreak_group_all$og_id
# Map from intervals with outbreak to all intervals
u_int_obint <- which(map_int_outbreak_group>0)
map_int_obint <- map_dbl(
  intervals$int_id, 
  function(x) {
    res <- which(u_int_obint == x)
    if (length(res) > 0) {
      res
    } else {
      0
    }
  })

# Make maps from tests to outbreak info
test_outbreak_data <- proc_test_data %>%
  select(-visite_1) %>% 
  left_join(outbreak_group_comb) %>% 
  # Set to outside outbreak transmission if visite 1
  mutate(og_id = case_when(
    is.na(og_id) ~ 0,
    test_date < (visite_1 - 2) ~ 0,
    in_outbreak ~ og_id,
    T ~ 0
  )) 

# Make map to outbreak id (0: out of outbreak)
map_test_outbreak <- test_outbreak_data %>% 
  pull(og_id)

# Make map to oubtreak interval
map_test_interval <- array(NA, dim = nrow(proc_test_data))

for (i in 1:nrow(proc_test_data)) {
  res <- intervals %>% 
    filter(
      sugar_id == proc_test_data$sugar_id[i],
      (TL - 2) <= proc_test_data$test_date[i],
      (TR + 7) >= proc_test_data$test_date[i]
    ) %>% 
    pull(int_id)
  
  # From interval to outbreak interval
  if (length(res) == 1) {
    map_test_interval[i] <- res
  } else if (length(res) > 1) {
    map_test_interval[i] <- min(res)
  } else {
    map_test_interval[i] <- 0
  }
}

# Check if assignment went well
testthat::expect_false(any(map_test_interval == 0))

# Make map to outbreak interval
map_test_obint <- array(NA, dim = nrow(proc_test_data))

for (i in 1:nrow(proc_test_data)) {
  if (!proc_test_data$in_outbreak[i]) {
    map_test_obint[i] <- 0
  } else {
    res <- intervals %>% 
      filter(
        sugar_id == proc_test_data$sugar_id[i],
        (TL - 2) <= proc_test_data$test_date[i],
        (TR + 7) >= proc_test_data$test_date[i],
        !str_detect(outbreak_id, "-")
      ) %>% 
      pull(int_id)
    
    # From interval to outbreak interval
    if (length(res) == 1) {
      map_test_obint[i] <- which(u_int_obint == res)
    } else if (length(res) > 1) {
      map_test_obint[i] <- which(u_int_obint ==  min(res))
    } else {
      map_test_obint[i] <- 0
    }
  }
}

# Check if assignment went well
testthat::expect_false(any(map_test_obint > 0 & map_test_outbreak == 0))
testthat::expect_false(any(map_test_obint == 0 & map_test_outbreak > 0))

# Get the variant period for each test. u_variants already defined above
# for assinging variant periods to serology intervals
map_test_variant <- map_chr(proc_test_data$test_date, ~ getVariant(.)) %>% 
  map_dbl(function(x){which(u_variants == x)})


# O. Control data ------------------------------------------------------------

# Michielin, et al., 2023. Clinical sensitivity and specificity of a 
# high-throughput microfluidic nano-immunoassay combined with capillary blood
# microsampling for the identification of anti-SARS-CoV-2 Spike IgG serostatus. 
# Plos one, 18(3), p.e0283149. https://doi.org/10.1371/journal.pone.0283149
dt_control <- 11 * 30

control_data <- readxl::read_xlsx("data/09_other_data/journal.pone.0283149.s004.xlsx") %>% 
  janitor::clean_names()

# Positives
pos_dat <- control_data %>% 
  filter(!is.na(nia_mitra_result)) %>% 
  filter(roche_result == "pos") %>% 
  count(nia_mitra_result)

pos_dat_stats <- tibble(
  prop = 1 - min(pos_dat$n)/sum(pos_dat$n),
  lo = Hmisc::binconf(x = max(pos_dat$n), n = sum(pos_dat$n))[2],
  hi = Hmisc::binconf(x = max(pos_dat$n), n = sum(pos_dat$n))[3],
  age_cat = "adult",
  group = "Michielin et al.",
  dt = dt_control
)

saveRDS(pos_dat_stats, "generated_data/pos_dat_stats.rds")

neg_dat <- control_data %>% 
  filter(!is.na(nia_mitra_result)) %>% 
  filter(roche_result == "neg") %>% 
  count(nia_mitra_result)


# P. Baseline data SP2 -------------------------------------------------------

# Get seropositivity data in the general population in Dec 2020 (Stringhini et al., 2021)
baseline_data <- getBaselineSeroPopData()

# Compute corresponding age categories for study participants
age_data <- sero_long %>% 
  distinct(sugar_id) %>% 
  inner_join(
    households %>% 
      select(sugar_id, age, age_cat_simple = age_cat)
  ) %>% 
  mutate(
    age_cat = map_chr(age,  function(x) {
      with(baseline_data, age_cat[age_L <= x & x <= age_R])
    })
  ) %>% 
  arrange(sugar_id)

# Check that ordering of praticipants is consistent
testthat::expect_identical(age_data$sugar_id[1], intervals$sugar_id[1])

# Data for stan
u_age_cats <- unique(age_data$age_cat)
N_age_cat <- length(u_age_cats)
map_part_age_cat <- map_dbl(age_data$age_cat, ~ which(u_age_cats == .))

# Select data for age categories in data
baseline_pos <- map_dbl(u_age_cats, function(x) with(baseline_data, n_pos[age_cat == x]))
N_baseline <- map_dbl(u_age_cats, function(x) with(baseline_data, n_tot[age_cat == x]))

# Save for later use

list(
  u_age_cats = u_age_cats,
  baseline_data = baseline_data,
  obs_counts = tibble(age_cat = u_age_cats[map_part_age_cat]) %>% 
    count(age_cat, name = "n_tot")
) %>% 
  saveRDS("generated_data/baseline_seroprev_data.rds")

# Q. Define intervals by cohort for generation --------------------------------

group_events <- sero_long %>% 
  filter(str_detect(outbreak_id, "vacations|variant", negate = T)) %>% 
  group_by(school, group, term, outbreak_id, event) %>% 
  summarise(date = mean(date)) 

# Define cohort
cohort_events <- group_events %>% 
  ungroup() %>% 
  mutate(
    prev_group = map_chr(as.character(group), ~ getPerviousGroup(.)),
    cohort = case_when(
      term == "21/22" & prev_group != "first" ~ str_c(prev_group, "20/21", sep = "-"),
      T ~ str_c(group, term, sep = "-")
    )
  )


# Add 20/21 basline data for all as before
missing_cohort_baseline_20_21 <- cohort_events %>% 
  distinct(cohort, school, outbreak_id) %>%
  group_by(cohort) %>%
  mutate(has_baseline_20_21 = any(str_detect(outbreak_id, "baseline") & str_detect(outbreak_id, "20/21"))) %>%
  filter(!has_baseline_20_21) %>% 
  distinct(cohort, school)

# Make explicit missing data for baselines 20/21
missing_data_cohort_baseline_20_21 <- missing_cohort_baseline_20_21 %>% 
  inner_join(baseline_metadata) %>% 
  rename(date = TL) %>%
  select(-TR) %>% 
  mutate(outbreak_id = str_c(event, school, term, sep = "-"))

# Fake data to model changes of SARS-CoV-2 variants dates for everyone
cohort_variant_changes <- map_df(
  c("delta", "omicron"), 
  function(x) {
    expand.grid(
      cohort = unique(cohort_events$cohort),
      term = getAllTerms()
    ) %>% 
      as_tibble() %>% 
      inner_join(
        getVariantDates(x) %>% 
          select(variant, date = TL) %>% 
          addTerms()
      ) %>% 
      mutate(
        event = "variant_change",
        outbreak_id = str_c(event, term, variant, sep = "-")
      ) %>% 
      select(-variant)
  })


# Fake data for vacations for everyone
cohort_vacations <- cohort_events %>% 
  distinct(cohort, school) %>% 
  inner_join(
    # This is vacation dates for all participants
    getVacations() %>% 
      addTerms(date_col = "TL") %>% 
      inner_join(
        expand.grid(name = getVacations()$name,
                    school = getAllSchools())
      )
  ) %>%
  mutate(name = str_c("vacations_", name)) %>% 
  rename(event = name) %>% 
  pivot_longer(
    cols = c("TL", "TR"),
    names_to = "what",
    values_to = "date"
  ) %>% 
  mutate(
    event = case_when(
      what == "TL" ~ str_c(event, "_start"),
      T ~ str_c(event, "_end")
    )
  ) %>% 
  mutate(outbreak_id = str_c(event, school, term, sep = "-")) %>% 
  distinct() %>% 
  select(-what)

# Add final date
max_date <- max(cohort_events$date)
cohort_final_date <- expand_grid(
  cohort = unique(cohort_events$cohort),
  date = max_date,
  event = "last_date"
) %>% 
  addTerms() %>% 
  mutate(outbreak_id = str_c(event, cohort, term, sep = "-")) %>% 
  filter(!(cohort %in% (cohort_events %>% filter(date == max_date) %>% pull(cohort))))

full_cohort_data <- cohort_events %>% 
  bind_rows(missing_data_baseline_20_21,
            cohort_vacations,
            cohort_variant_changes,
            cohort_final_date) %>% 
  mutate(
    event_type = map_chr(outbreak_id, ~ makeEventType(., simplify = F, add_terms = TRUE)),
    event_type = factor(event_type, levels = names(getEventTypeColors(with_terms = T)))
  ) %>% 
  ungroup() %>% 
  arrange(cohort, date) %>% 
  mutate(time_id = row_number())

cohort_groups <- distinct(cohort_events, cohort, school, group, term)

# First schools
for (i in which(is.na(full_cohort_data$school))) {
  full_cohort_data$school[i] <- cohort_groups$school[cohort_groups$cohort == full_cohort_data$cohort[i]] %>% first()
}

# Then groups
for (i in which(is.na(full_cohort_data$group))) {
  full_cohort_data$group[i] <- cohort_groups$group[cohort_groups$cohort == full_cohort_data$cohort[i]] %>% first()
}


if (do_plots) {
  full_cohort_data  %>% 
    ggplot(aes(x = date, y = cohort)) +
    geom_line(aes(group = cohort), lwd = .2) +
    geom_point(aes(color = event_type)) +
    theme_bw() +
    scale_color_manual(values = getEventTypeColors(with_terms = T)) +
    facet_grid(school ~ ., scales = "free") +
    ggtitle("Analysis longitudinal dataset")
}

cohort_intervals <- makeIntervals(df = full_cohort_data, 
                                  group_var = "cohort")

# Make necessary mappings for stan
N_cohorts <- nrow(full_cohort_data)
L_cohorts <- nrow(cohort_intervals)
X_cohort_int <- model.matrix(~ variant + vacations, 
                             data = cohort_intervals %>% 
                               select(variant, vacations))

# number of unique cohorts
u_cohorts <- full_cohort_data %>% 
  distinct(cohort) %>%
  pull(cohort)

M_cohorts <- length(u_cohorts)

# Map from cohort interval to variant
map_cohort_int_variant <- map_chr(cohort_intervals$TL, ~ getVariant(.)) %>% 
  map_dbl(function(x){which(u_variants == x)})

# Sarts end ends of intervals for each cohort
start_ends_chort_df <- cohort_intervals %>% 
  group_by(cohort) %>% 
  summarise(start = min(int_id),
            end = max(int_id))

# Extract vectors
starts_cohort <- start_ends_chort_df$start
ends_cohort <- start_ends_chort_df$end

cohort_outbreak_group_all <- cohort_intervals %>% 
  left_join(outbreak_group_comb) %>% 
  replace_na(list(og_id = 0))

map_cohort_int_outbreak_group <- cohort_outbreak_group_all$og_id

int_cohort_starts <- cohort_intervals$tid_L
int_cohort_ends <- cohort_intervals$tid_R
cohort_times <- as.numeric(full_cohort_data$date)

# Map from interval to unique interval combination
map_cohort_int_int_comb <- apply(X_cohort_int, 1, function(x){
  which(
    map_lgl(1:nrow(X_int_comb), ~ identical(x, X_int_comb[.,]))
  )
})

u_cohort_int_obint <- which(map_cohort_int_outbreak_group>0)
map_cohort_int_obint <- map_dbl(
  cohort_intervals$int_id, 
  function(x) {
    res <- which(u_cohort_int_obint == x)
    if (length(res) > 0) {
      res
    } else {
      0
    }
  })

# Map to age categories (for children only)
# TODO

# Save
saveRDS(full_cohort_data, here("generated_data/full_cohort_data.rds"))

# R. Data for stan -----------------------------------------------------------

data <- list(
  # Integers
  N = sero_long %>% nrow(),
  M = M,
  J = sero_long %>% distinct(outbreak_id) %>% nrow(),
  S = sero_long %>% distinct(school) %>% nrow(),
  G = outbreak_group_comb %>% nrow(),
  K = sero_long %>% distinct(term) %>% nrow(),
  L_outbreaks = sum(map_int_outbreak_group>0),
  N_tests = N_tests,
  N_obs = N_obs,
  V = V,
  L = L,
  E = max(map_obs_event),
  N_all_tests = nrow(proc_test_data),
  max_n_intervals = max_n_intervals,
  Z_comb = 0,
  n_sens_covar = 0,
  n_covar_int = n_covar_int,
  Z_int_comb = Z_int_comb,
  M_cohorts = M_cohorts,
  L_cohorts = L_cohorts,
  N_cohorts = N_cohorts,
  N_age_cat = N_age_cat,
  
  # Maps
  map_obs_participant = makeMap(sero_long, "sugar_id"),
  map_obs_interval = map_obs_interval,
  starts = starts,
  ends = ends,
  int_starts = int_starts,
  int_ends = int_ends,
  map_int_variant = map_int_variant,
  map_int_outbreak_group = map_int_outbreak_group,
  map_test_variant = map_test_variant,
  map_test_outbreak = map_test_outbreak,
  map_outbreak_variant = map_outbreak_variant,
  map_obs_event = map_obs_event,
  map_int_obint = map_int_obint,
  map_int_vaccination = map_int_vaccination,
  map_test_obint = map_test_obint,
  map_test_interval = map_test_interval,
  map_int_int_comb = map_int_int_comb,
  map_outbreak_group_int_comb = map_outbreak_group_int_comb,
  map_cohort_int_variant = map_cohort_int_variant,
  starts_cohort = starts_cohort,
  ends_cohort = ends_cohort,
  map_cohort_int_outbreak_group = map_cohort_int_outbreak_group,
  int_cohort_starts = int_cohort_starts,
  int_cohort_ends = int_cohort_ends,
  map_cohort_int_int_comb = map_cohort_int_int_comb,
  map_part_age_cat = map_part_age_cat,
  map_outbreak_int_variant = map_outbreak_int_variant,
  
  # Data
  y = y,
  times = times,
  t_pos_test = t_pos_test,
  prev_pos = sero_long$prev_pos[ind_prev_tests],
  y_tests = proc_test_data$test_result == "positive",
  cohort_times = cohort_times,
  
  # Indices
  ind_sero_obs = ind_sero_obs,
  ind_only_sero = setdiff(ind_sero_obs, ind_prev_tests),
  ind_prev_tests = ind_prev_tests,
  
  # Controls for serology sens/spec
  control_tp = pos_dat$n[pos_dat$nia_mitra_result == "pos"],
  control_fp = neg_dat$n[neg_dat$nia_mitra_result == "pos"],
  N_pos_control = sum(pos_dat$n),
  N_neg_control = sum(neg_dat$n),
  sens_pcr = .9,
  spec_pcr = .99,
  
  # Generated quantities
  dt_sens_pred = dt_sens_pred,
  N_sens_pred = length(dt_sens_pred),
  T_max = 500,
  
  # Priors
  mu_log_lambda = -1,
  sd_log_lambda = 2,
  mu_log_alpha = -2,
  sd_log_alpha = 2.5,
  mu_mu_log_alpha = -2,
  sd_mu_log_alpha = 2.5,
  # Covariate matrices
  # X_sens = X_sens,
  # n_sens_covar = ncol(X_sens),
  dt_control = dt_control,
  
  # Other
  inf_window = 7,
  inf_window_outbreak = 30,
  dt_exposure = 3,
  vaccinated = vaccinated,
  
  # Covariates
  X_sens_comb =  array(NA, dim = c(nrow(X_sens_comb), 0)),
  # X_sens_comb = X_sens_comb,
  # Z_comb = nrow(X_sens_comb),
  X_sens = array(NA, dim = c(nrow(X_sens), 0)),
  X_int = X_int,
  X_int_comb = X_int_comb,
  X_int_outbreak_group = X_int_outbreak_group,
  X_cohort_int = X_cohort_int,
  
  # Baseline in dec 2020 https://doi.org/10.1016/S1473-3099(21)00054-2
  # Stringhini et al. 2021s
  baseline_pos = baseline_pos,
  N_baseline = N_baseline,
  sens_roche = .98,    # Sensitivity of Roche-S test done in Dec 2022
  spec_roche = .99,     # Specificity of Roche-S test done in Dec 2022
  
  
  # Integral tolerance
  tol = 1e-5
)

# Save data for later use
saveRDS(outbreak_group_comb, "generated_data/outbreak_group_comb.rds")
saveRDS(outbreak_metadata, "generated_data/outbreak_metadata_stat_model.rds")

# S. Run stan ----------------------------------------------------------------

# Compile the stan model
model <- cmdstan_model("analysis/stan/outbreak_seroprev_forward_simple.stan")

# Draw samples from the posterior
cmdstan_fit <- model$sample(data = data,
                            chains = 4,
                            parallel_chains = 4,
                            iter_warmup = 250,
                            iter_sampling = 1000,
                            max_treedepth = 12,
                            init = .1, 
                            refresh = 10
                            # save_warmup = TRUE,
                            # output_dir = "generated_data",
                            # output_basename = "test_pooled"
)

cmdstan_fit$save_object("generated_data/stan_model_fit_non_pooled_simple.rds")

# Compute generated quantities
genquant <- model$generate_quantities(
  fitted_params = cmdstan_fit,
  data = data,
  parallel_chains = 4
)

genquant$save_object("generated_data/stan_model_genquant_non_pooled_simple.rds")


# T. Save objects for pomp model run -----------------------------------------

list(
  T_max = data$T_max,
  beta_draws = cmdstan_fit$draws("beta") %>% 
    as_draws() %>% 
    as_draws_df() %>% 
    as_tibble()
) %>% 
  saveRDS("generated_data/time_sens_beta_draws_v2.rds")

pred_prob <- cmdstan_fit$summary("p", .cores = 4)
saveRDS(pred_prob, here("generated_data/pred_prob.rds"))
