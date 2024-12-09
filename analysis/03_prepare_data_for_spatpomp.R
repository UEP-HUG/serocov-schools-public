# This script prepares the data for the pomp model

# Preamble ----------------------------------------------------------------
library(stringr)
library(magrittr)
library(spatPomp)
library(optparse)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(lubridate)

source("analysis/utils.R")
source("analysis/pomp_utils.R")

Sys.setenv("REDO_DATA" = FALSE)

option_list <- list(
  make_option(c("-s", "--schools"), 
              default = "SA", 
              action ="store", type = "character", help = "Schools filter"),
  make_option(c("-t", "--terms"), 
              default = "20/21, 21/22", 
              action ="store", type = "character", help = "Terms filter"),
  make_option(c("-m", "--model_version"), 
              default = "gen", 
              action ="store", type = "character", help = "Model type"),
  make_option(c("-g", "--use_genlik"), 
              default = TRUE, 
              action ="store", type = "logical", help = "Make data rport"),
  make_option(c("-d", "--make_data_report"), 
              default = FALSE, 
              action ="store", type = "logical", help = "Make data rport")
)

opt <- parse_args(OptionParser(option_list = option_list))

filter_schools <- parseFilter(opt$schools) 
filter_terms <- parseFilter(opt$terms)

run_file <- makeRunDataFile(filter_schools = filter_schools, 
                            filter_terms = filter_terms,
                            model_version = opt$model_version,
                            use_gen_lik = opt$use_genlik)

pomp_file <- makePompObjFile(run_file = run_file) 
pomp_fit_file <- makePompFitFile(run_file = run_file)
traceplot_file <- makeTraceplotFile(run_file = run_file)
simplot_file <- makeSimplotFile(run_file = run_file)
markdown_file <- makeMarkdownFile(run_file = run_file)

if (filter_schools[1] == "all") {
  filter_schools <- NULL
}

# A. School group composition ------------------------------------------------

# We here use all epi data for all participants with consent. We start by only
# modeling participants within schools (not households)

households <- getDataHouseholds()
school_groups <- getDataGroups()

# Filter participants with consent and in schools (exclude households)
included_participants <- households %>% 
  # !! this implicetely eliminates household members
  inner_join(school_groups) %>% 
  # !! Remove admin which are not in contact with children
  filter(group != "G1") %>% 
  subsetData(schools = filter_schools,
             terms = filter_terms) %>% 
  distinct(sugar_id, hug_codbar, age_cat, term, school, group)

# Get school group compositions
group_composistions <- getSchoolCounts() %>% 
  subsetData(schools = filter_schools,
             terms = filter_terms) %>% 
  # !! Remove admin
  filter(group != "G1")

# Expand participant ids to account for unsampled participants
full_sugar_ids <- expandIds(obs = included_participants, 
                            all_counts = group_composistions) %>% 
  mutate(added_id = str_detect(sugar_id, "^M"),
         hug_codbar = coalesce(hug_codbar, sugar_id)) %>% 
  left_join(getDataHouseholds() %>% 
              select(sugar_id, hh_id)) %>% 
  replace_na(list(hh_id = "missing")) %>% 
  arrange(term, school, group, sugar_id)

# Unique set of sugar ids to use in pomp Model
u_sugar_ids <- full_sugar_ids %>% distinct(sugar_id) %>% pull(sugar_id)
U <- length(u_sugar_ids)    # number of distinct participants
u_sugar_ids_num <- makeUIdsNum(u_sugar_ids)

sugar_ids_mapping <- tibble(
  sugar_id = u_sugar_ids,
  id = u_sugar_ids_num
)

# B. Epi data ----------------------------------------------------------------

# Test data - - -
test_data <- readRDS("generated_data/combined_test_data.rds") %>% 
  filter(sugar_id %in% included_participants$sugar_id)

test_data_proc <- test_data %>% 
  mutate(day = dateToDay(test_date)) %>% 
  select(sugar_id, date = test_date, day, n_test, n_pos)

# Serological data - - -
sero_data <- readRDS("generated_data/all_serologies.rds")  %>% 
  filter(sugar_id %in% included_participants$sugar_id)

sero_data_proc <- sero_data %>% 
  mutate(day = dateToDay(date_rdv)) %>% 
  select(sugar_id, date = date_rdv, day, sero_result = nia_epfl_interp)  %>% 
  mutate(sero_result = ifelse(sero_result == "positive", 1, 0))

# C. Adjacency ---------------------------------------------------------------
# Make adjacency between participants
# Adjacency types:
#   - 0: in same group
#   - 1: in same school

adj_list_df <- buildAdjacencyDF(full_sugar_ids = full_sugar_ids,
                                u_sugar_ids = u_sugar_ids) %>% 
  rename(year = term) %>% 
  arrange(year, from, to, type) %>% 
  distinct() %>% 
  # !! Keep only closest connection
  group_by(from, to) %>% 
  slice_min(type) %>% 
  ungroup() %>% 
  filter(from != to) %>% 
  arrange(year, from, to) %>% 
  mutate(from_chr = from,
         to_chr = to,
         from = map_dbl(from, ~ which(u_sugar_ids_num == .) - 1),
         to = map_dbl(to, ~ which(u_sugar_ids_num == .) - 1))

adj_epi <- list(
  adj_list = select(adj_list_df, from, to),
  adj_type = adj_list_df$type,
  start_list =  getStartEndAdj(adj_list_df, "start") - 1,
  end_list = getStartEndAdj(adj_list_df, "end") - 1
)

# Check if everyone is connected
testthat::expect_equal(sum(!(u_sugar_ids_num %in% adj_list_df$from_chr)), 0)
testthat::expect_equal(sum(!(u_sugar_ids_num %in% adj_list_df$to_chr)), 0)

# D. Masks ----------------------------------------------------------------

# Make the information on whether participants are in school or not for the 
# model to consider.
# mask = 0 indicates that participant is excluded from FOI computations
masks_df <- makeMaskDF(full_sugar_ids = full_sugar_ids,
                       u_sugar_ids = u_sugar_ids,
                       test_data = test_data_proc,
                       group_closures = getGroupColsures(),
                       vacations = getVacations(),
                       quarantine_duration = getQuarantinePolicy()) 

# Make to matrix for pomp
masks_info <- makeMaskInfo(masks_df = masks_df)

# E. Genetic data ------------------------------------------------------------

# Phylo metadata
phylo_metadata <- getPhyloMetadata()

# Plot sample collection dates
phylo_data <- getPhyloData(phylo_metadata = phylo_metadata,
                           sugar_ids_mapping = sugar_ids_mapping, 
                           terms = filter_terms,
                           adj_list_df = adj_list_df)

# Make matrix of probability of time difference of infections
max_kappa <- 5   # maximum number of generation times
max_delay <- 100    # maximum number of days of delay between generations
w <- c(1, 1, .75, .5, .25, .25, 0)
w_mat <- makeWMat(max_kappa, max_delay, w)

sample_dates <- purrr::map_dbl(phylo_data$u_seq_all, ~ phylo_metadata$collection_date[phylo_metadata$seq_id == .])

seq_data <- tibble(
  date = sample_dates,
  seq = 0:(length(phylo_data$u_seq_all) - 1)    # remove 1 for zero-indexing in C
) %>% 
  mutate(date = as.Date(date, origin = "1970-01-01"),
         sugar_id = map_chr(phylo_data$u_seq_all, ~ phylo_metadata$sugar_id[phylo_metadata$seq_id == .]))

# Make vecotr with indext of sequence
seq_covar <- rep(-999, U)
time_seq <- rep(-999, U)

for (i in 1:nrow(seq_data)) {
  j <- which(sugar_ids_mapping$sugar_id == seq_data$sugar_id[i])
  seq_covar[j] <- seq_data$seq[i]
  time_seq[j] <- dateToTime(seq_data$date[i])
}

# Delay before which to consider sequence
seq_delay <- 10/365


gen <- append(phylo_data,
              list(
                w_mat = w_mat,
                mu = getMutationRate(),
                max_kappa = max_kappa,
                seq_covar = seq_covar,
                time_seq = time_seq,
                seq_delay = seq_delay
              ))

# F. Vaccination data --------------------------------------------------------

# Get vaccination data
vacc_data <- readRDS("generated_data/unified_vaccination_dates.rds") %>% 
  mutate(first_vacc_date = coalesce(d1_date, d2_date, d3_date)) %>% 
  filter(sugar_id %in% u_sugar_ids)

# Define vaccination times. Set to -1 when no vaccination date available
t_vacc <- rep(-1, U)

for (i in 1:U) {
  if (sugar_ids_mapping$sugar_id[i] %in% vacc_data$sugar_id) {
    t_vacc[i] <- dateToTime(vacc_data$first_vacc_date[vacc_data$sugar_id == sugar_ids_mapping$sugar_id[i]])
  }
}

# Fake date to define the susceptibility reduction equivalent to vaccination
t_vacc_fake <- dateToTime(as.Date("2021-01-30"))

vaccination <- list(
  t_vacc_fake = t_vacc_fake,
  t_vacc = t_vacc
)

# G. Test results data ---------------------------------------------------------

data <- test_data_proc %>% 
  full_join(
    sero_data_proc %>% 
      group_by(sugar_id, date) %>% 
      slice_max(sero_result) %>% 
      ungroup()) %>% 
  full_join(seq_data) %>% 
  distinct() %>% 
  mutate(participant = map_chr(sugar_id, ~ u_sugar_ids_num[which(u_sugar_ids == .)]),
         time = dateToTime(date)) %>% 
  filterRunDataForPomp(terms = filter_terms) %>%
  ungroup() %>% 
  select(participant, time, n_pcr = n_test, pcr_pos = n_pos, sero = sero_result, seq) %>% 
  distinct() %>% 
  arrange(participant, time) 

# Set pcr for fake serologies
data <- data %>% 
  mutate(n_pcr = case_when(
    time == min(data$time) & !is.na(sero) ~ 1,
    T ~ n_pcr),
    pcr_pos = case_when(
      time == min(data$time) & !is.na(sero) ~ 0,
      T ~ pcr_pos
    )
  )

# Add age cat as data
age_cat <- sugar_ids_mapping %>% 
  left_join(households %>% select(sugar_id, age_cat)) %>% 
  mutate(age_cat = case_when(is.na(age_cat) & str_detect(sugar_id, "-a-") ~ 1,
                             is.na(age_cat) & str_detect(sugar_id, "-c-") ~ 0,
                             age_cat == "adult" ~ 1,
                             TRUE ~ 0))

# Add sequence prevalence -------------------------------------------------

# Load predicted proportion of SARS-CoV-2 variants in Geneva from SP4 analysis
prob_variants <- readRDS("data/09_other_data/GE_variant_probs.rds") %>% 
  # Combine all Omicron variants
  mutate(variant = str_extract(variant, "Alpha|Delta|Omicron") %>% str_to_lower()) %>% 
  filter(!is.na(variant)) %>% 
  group_by(week, variant) %>% 
  summarise(prob = sum(mean)) %>% 
  ungroup() %>% 
  filter(week > "2021-01-01") %>% 
  mutate(variant = str_c("prob_", variant)) %>% 
  #!! shift forward by 10 days to capture change in transmission
  mutate(week = week + 10)

covar <- prob_variants %>% 
  pivot_wider(values_from = "prob",
              names_from = "variant") %>% 
  mutate(time = dateToTime(week)) %>% 
  select(-week) %>% 
  arrange(time) %>% 
  mutate(participant = "001") %>% 
  select(participant, time, unique(prob_variants$variant))

# Expand to cover all times 
max_time_data <- max(data$time)
max_time_covar <- max(covar$time)

covar <- bind_rows(
  covar,
  covar %>% 
    slice_max(time) %>% 
    mutate(time = max_time_data + getDeltaT())
)

testthat::expect_lt(max_time_data, max(covar$time))

# expand to all unit names
covar <- map_df(sugar_ids_mapping$id,
                function(x) {
                  covar %>%
                    mutate(participant = x)
                }) %>%
  select(participant, time , contains("prob")) %>%
  arrange(time, participant)



# G1. Add time of sequences ----------------------------------------------------
# -1 indicates that no sequences were taken for that participant
# gen$time_seq <- covar %>% 
#   group_by(participant) %>% 
#   arrange(seq_covar) %>% 
#   slice(1) %>% 
#   mutate(time = case_when(!is.na(seq_covar) ~ time,
#                           T ~ -1)) %>% 
#   pull(time)


# Time varying serology sensitivity ---------------------------------------
# Ensure consistent draws of sensitivity parameters
set.seed(1532497)

# Load fit results of statistical model [origin: XX_fit_stat_model.R]
beta_draws_data <- readRDS("generated_data/time_sens_beta_draws.rds")

# Extract the draws themselves
time_sens_coefs <- beta_draws_data$beta_draws %>% 
  select(contains("beta")) %>%
  # Take only 50 draws for computation time
  sample_n(50) %>% 
  as.matrix()

# Setup auxilliary data used to infer sensitivity in statistical model run

T_max <- beta_draws_data$T_max
t_log_raw <- log(1:T_max)
t_log_mean <- mean(t_log_raw)
t_log_sd <- sd(t_log_raw)


sero_sens <- list(
  t_log_mean = t_log_mean,
  t_log_sd = t_log_sd,
  time_sens_coefs = time_sens_coefs
)

# H. Globals -----------------------------------------------------------------

test_perf <- list(
  sens_pcr = .9,
  spec_pcr = .99,
  sens_sero = .85,
  spec_sero = .99
)

liks <- list(
  use_test_lik = 1,
  use_sero_lik = 1,
  use_gen_lik = as.numeric(opt$use_genlik)
)

# Get natural history parameters
nat_hist <- getNatHistParams()


# Seroprev at baseline and endline
seroprev <- getSeroprevDraws(n_draws = 100)

globals_parlist <- list(
  test_perf = test_perf,
  adj_epi = adj_epi,
  gen = gen,
  liks = liks,
  nat_hist = nat_hist,
  masks_info = masks_info,
  vaccination = vaccination,
  sero_sens = sero_sens,
  seroprev = seroprev
)

# I. Initial parameters ------------------------------------------------------

init_params = c(
  "lambda_wt" = .05,        # community FOI
  "lambda_a" = 1,       
  "lambda_d" = 1,       
  "lambda_o" = 1,        
  "lambdaA_wt" = 15,        # same group FOI
  "lambdaA_a" = 15,      
  "lambdaA_d" = 15,        
  "lambdaA_o" = 15,       
  "lambdaB_wt" = 3,        # same class FOI
  "lambdaB_a" = 0,       
  "lambdaB_d" = 0,        
  "lambdaB_o" = 0,       
  "S_0" = .9,        # probability of susceptible at baseline
  "C_0" = .1         # probability of susceptible at baseline
)  


base_pars <- names(init_params)
fixed_params <- getFixedParams(base_pars = base_pars,
                               terms = filter_terms) 

# J. Blocks for ibpf ---------------------------------------------------------
# Make blocks based on groups

block_list <- makeBlockList(full_sugar_ids = full_sugar_ids,
                            u_sugar_ids = u_sugar_ids,
                            u_sugar_ids_num = u_sugar_ids_num)

# Compute statistics of numbers of children and adults by group and 
# expand to all participants to give as data
block_ages <- makeBlockAges(block_list = block_list,
                            age_cat = age_cat)

# Add to data
data <- data %>%
  # Add fake data at first time point
  bind_rows(
    map_df(
      seq_along(getSeroprevalenceDates()), 
      function(x) {
        block_ages %>% 
          mutate(time = dateToTime(getSeroprevalenceFakeDates()[x]))
      }
    )
  ) %>%
  distinct() %>%
  arrange(participant, time) %>% 
  group_by(participant, time) %>% 
  summarise(n_pcr = getNonNA(n_pcr),
            pcr_pos = getNonNA(pcr_pos),
            seq = getNonNA(seq),
            sero = getNonNA(sero),
            age_cat = getNonNA(age_cat),
            n_block = getNonNA(n_block)
  ) %>% 
  ungroup() %>%
  arrange(time, participant)

# Map from units to blocks
map_participant_block <- makeMapParticipantBlock(block_list = block_list)

# Update
globals_parlist$seroprev$map_participant_block <- map_participant_block
globals_parlist$seroprev$map_participant_age_cat <- age_cat$age_cat

# K. Save data ---------------------------------------------------------------

# Save all model run data (to be committed to github for yggdrasil)
saveRDS(
  list(
    adj_list_df = adj_list_df,
    covar = covar,
    group_composistions = group_composistions,
    sugar_ids_mapping = sugar_ids_mapping,
    full_sugar_ids = full_sugar_ids,
    phylo_metadata = phylo_metadata,
    globals_parlist = globals_parlist,
    data = data,
    init_params = init_params,
    fixed_params = fixed_params,
    block_list = block_list,
    opt = opt,
    U = U,
    times = data$time
  ),
  file = run_file
)

cat("---- Saved data to", run_file, "\n")

# Save test/sero data to separate file for internal report (not committed)
saveRDS(
  list(
    test_data_proc = test_data_proc,
    sero_data_proc = sero_data_proc
  ),
  file = str_replace(run_file, "\\.rds", "_test_data.rds")
)

# Make data report
if(opt$make_data_report) {
  rmarkdown::render(
    'notebooks/model_data_report_template.Rmd', 
    output_file = markdown_file, 
    output_dir = "notebooks/data_reports/", 
    params = list(run_file = run_file) 
  )
}

