# This script makes the final figures for the paper


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(posterior)
library(cmdstanr)
library(here)
library(cowplot)
library(tidytree)
library(ggtree)
library(treeio)

dir("analysis/", pattern = "utils", full.names = T) %>% 
  walk(~source(.))

Sys.setlocale("LC_ALL","C")
Sys.setenv("REDO_DATA" = "FALSE")

# Load data ---------------------------------------------------------------

## Epi data ----
# All available data epi data
all_data <- getDataCombined() 

# Outbreak metadata
outbreak_metadata <- getOutbreakMetadata() %>% 
  group_by(school, outbreak_id) %>% 
  summarise(date = min(visite_1)) %>% 
  arrange(date) %>% 
  correctOutbreakID()

# Epi data by outbreak
outbreak_data <- getDataOutbreak() %>% 
  correctOutbreakID()


## Contextual data ----
# Predicted proportion of SARS-CoV-2 variants in Geneva from SP4 analysis
prob_variants <- readRDS("data/09_other_data/GE_variant_probs.rds") %>% 
  # Combine all Omicron variants
  mutate(variant = str_extract(variant, "Alpha|Delta|Omicron") %>% str_to_lower()) %>% 
  filter(!is.na(variant)) %>% 
  group_by(week, variant) %>% 
  summarise(prob = sum(mean),
            lo = sum(q5),
            hi = sum(q95)) %>% 
  ungroup() %>% 
  filter(week > "2021-01-01")

# Geneva cantonal epi data
gva <- read_csv("data/09_other_data/GVA_covid19_data.csv") %>% 
  janitor::clean_names() %>% 
  mutate(date = str_c(semaine, "-1") %>% as.Date("%y-%U-%u"))

# Vacation dates for plots
vacations <- getVacations()

## Genetic data ----
# Genetic metadata
phylo_metadata <- getPhyloMetadata() %>% 
  filter(coverage > .96) %>% 
  mutate(hh = !(sugar_id %in% school_groups$sugar_id)) %>% 
  addEpiWeek(date_col = "collection_date")


## Inference results ----
# Dynamic parameter estimates
param_estimates <- readRDS("generated_data/all_param_estimates.rds")

# State-space trajectories from smoothing distribution 
traj_stats <- readRDS("generated_data/pomp_traj_4e535e_sb_sa_sc_sd_20_21_21_22_gen_genlik_true_refined/pomp_traj_4e535e_sb_sa_sc_sd_20_21_21_22_gen_genlik_true_stats.rds") %>% 
  # Remove data for school SB which was not part of the study in the first year
  filter(!(school == "SB" & date < "2021-08-01")) 


## Scenario results ----
# Get the results directories of simulations ignoring the default scenario
scenario_dirs <- dir("generated_data/simulations/", full.names = T) %>% 
  str_subset("default", negate = TRUE)

# Statistics on cases averted
effect_stats <- computeScenarioStats(scenario_dirs = scenario_dirs, redo = FALSE) %>% 
  mutate(variant = str_to_title(variant) %>% str_replace("_", "-"),
         variant = factor(variant, levels = c("Omicron-like", "Delta-like", "Alpha-like")),
         setting = case_when(setting == "preschool" ~ "pre-school",
                             TRUE ~ "primary school"),
         scenario = case_when(str_detect(variant, "Omi") ~ "Very frequent",
                              str_detect(variant, "Del") ~ "Rare",
                              TRUE ~ "Very rare"),
         scenario =  factor(scenario, levels = c("Very rare", "Rare", "Very frequent")))



# Figure 1: Context -------------------------------------------------------
# Contextualize the study and provide data on outbreaks

## Data wrangling for figure ----
# Combine contextual data
context_data <-   bind_rows(
  prob_variants %>% 
    select(date = week, value = prob, lo, hi, variant) %>%
    mutate(what = "prob_variant"),
  gva %>% 
    select(date, value = nombre_cas_covid_19_ge) %>% 
    mutate(what = "weekly_cases")
) %>% 
  # Make labels for figure
  mutate(
    what_label = case_when(
      what == "prob_variant" ~ "Proportion\nSARS-CoV-2\nvariant",
      T ~ "Weekly\nconfirmed\ncases")
  ) 

# Use study period bounds to set limits on plots
date_limits <- getRunDateBounds(getAllTerms())

# Combine metadata
metadata <- bind_rows(
  outbreak_metadata,
  getBaselineDates() %>% 
    mutate(school = getSchoolUMapping()[school],
           date = date_min) %>% 
    select(school, date, term) %>% 
    filter(school != "SE")
) %>% 
  makeSetting()

# Compute statistics on outbreak epi data
outbreak_stats <- outbreak_data %>% 
  mutate(visit = case_when(what == "sero" & date < visite_2 ~ "v1",
                           what == "sero" & date > visite_2 ~ "v3",
                           TRUE  ~ "v1")) %>% 
  group_by(school, outbreak_id, visit, what)  %>% 
  computeTestStats() %>% 
  makeSetting() %>% 
  inner_join(outbreak_metadata) %>% 
  mutate(dat_date2 = case_when(what == "sero" & visit == "v3" ~ dat_date,
                               TRUE ~ date))

# Compute non-outbreak test statistics
test_data_non_outbreak <- getDataCombined() %>% 
  filter(what == "pcr/antigen") %>% 
  left_join(getOutbreakMetadata()) %>% 
  filter(is.na(visite_1) | date < visite_1 | date > visite_3) %>% 
  distinct(sugar_id, date, school, result) %>% 
  addEpiWeek() %>% 
  makeSetting() %>% 
  group_by(setting, epiweek_date) %>% 
  computeTestStats() %>% 
  mutate(what = "pcr/antigen")

# Baseline/endline stats
baseline_endline_data <- getDataCombined() %>% 
  filter(what == "sero") %>%
  inner_join(
    bind_rows(
      getBaselineDates() %>% 
        mutate(school = getSchoolUMapping()[school],
               when = "baseline"),
      getEndlineDates() %>% 
        mutate(school = getSchoolUMapping()[school],
               when = "endline")
    )
  ) %>% 
  # Keep only data within baseline dates
  filter(date >= date_min - 1, date <= date_max + 1)  %>% 
  makeSetting() %>% 
  group_by(setting, term, school, when, what) %>% 
  computeTestStats() %>% 
  mutate(what = "anti-S capillary\nserology") 

# Combine test data for plot
test_stats_combined <- bind_rows(
  outbreak_stats %>%
    addEpiWeek(date_col = "dat_date2"), 
  baseline_endline_data %>% 
    addEpiWeek(date_col = "dat_date"),
  test_data_non_outbreak
) %>% 
  mutate(
    what = case_when(
      what == "sero" ~ "anti-S capillary\nserology",
      what == "pcr/antigen" ~ "RT-PCR/antigenic",
      T ~ what)
  ) %>% 
  group_by(epiweek_date, what, setting) %>% 
  summarise(n_neg = sum(n) - sum(n_pos),
            n_pos = sum(n_pos)) %>% 
  pivot_longer(cols = c("n_neg", "n_pos"),
               values_to = "value",
               names_to = "result") %>% 
  mutate(variant = map_chr(epiweek_date, ~ getVariant(.))) %>% 
  mutate(result = factor(result, levels = c("n_neg", "n_pos"), labels = c("negative", "positive")),
         result2 = case_when(result == "positive" ~ variant,
                             TRUE ~ result),
         result2 = factor(result2, levels = c("negative", "alpha", "delta", "omicron"))) %>% 
  filter(result2 != "wildtype") 

ob_lwd <- .25
ob_lty <- c(2:5)

## Top panel: Context with SARS-CoV-2 variant prevalence and GVA epi curve ----
p_context <- context_data %>% 
  ggplot(aes(x = date, y = value)) +
  # Timing of outbreaks
  geom_vline(data = outbreak_metadata, 
             aes(xintercept = date, lty = school),
             lwd = ob_lwd) +
  # SARS-CoV-2 VOC prevalence
  geom_line(data = filter(context_data, 
                          what == "prob_variant",
                          !(what == "prob_variant" & value < 1e-2),
                          !(variant == "omicron" & date < "2021-08-01")),
            aes(color = variant), size = .8) +
  # GVA Epi data
  geom_bar(data = filter(context_data, what == "weekly_cases"),
           stat = "identity", fill = "darkgray") +
  facet_grid(what_label ~ ., scales = "free_y", switch = "y") +
  coord_cartesian(xlim = date_limits) +
  scale_fill_manual(values = getSchoolColors()) +
  scale_color_manual(values = getVariantColors(), getSchoolColors()) +
  scale_shape_manual(values = c(21, 24, 23)) +
  scale_linetype_manual(values = ob_lty) +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.switch.pad.grid = unit(1, units = "lines")) +
  labs(y = "", shape = "SARS-CoV-2\nvariant") +
  guides(color = "none", fill = "none", lty = "none") +
  scale_x_date(date_breaks = "2 months", date_labels = "%Y-%b")

## Middle panel: outbreak timings ----
p_outbreak <- metadata %>% 
  mutate(outbreak_id = factor(outbreak_id, levels = rev(levels(outbreak_metadata$outbreak_id)))) %>% 
  ggplot(aes(x = date, y = school, color = school))  +
  geom_line(lwd = .4) +
  # Outbreak timings
  geom_vline(data = outbreak_metadata, 
             aes(xintercept = date, color = school, lty = school),
             lwd = ob_lwd) +
  # Outbreak numbers
  geom_label(data = filter(metadata, !is.na(outbreak_id)),
             aes(label = outbreak_id), 
             size = 3) +
  # Filler
  geom_point(data = filter(metadata, is.na(outbreak_id)),
             aes(pch = term, fill = school),
             color = "white",
             size = 3) +
  facet_grid(setting ~ ., scales = "free", switch = "y") +
  coord_cartesian(xlim = date_limits) +
  scale_color_manual(values = getSchoolColors()) +
  scale_fill_manual(values = getSchoolColors()) +
  scale_linetype_manual(values = ob_lty) +
  scale_shape_manual(values = c(22, 25)) +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.switch.pad.grid = unit(1, units = "lines")) +
  labs(y = "Outbreak setting") +
  guides(color = "none", fill = "none", lty = "none", 
         shape = guide_legend("Baseline\nschool term", 
                              override.aes = list(color = "black")),
         label = guide_legend("Outbreak ID")) +
  scale_x_date(date_breaks = "2 months", date_labels = "%Y-%b")


## Bottom panel: SARS-CoV-2 test results ---- 
p_stats <- test_stats_combined %>% 
  ggplot(aes(x = epiweek_date)) +
  # Test results
  geom_bar(stat = "identity", aes(y = value, fill = result2)) +
  # Filler
  geom_point(data = tibble(y = -40, x = as.Date("2021-03-30"), what = "RT-PCR/antigenic"), alpha = 0,
             aes(x = x, y = y)) +
  # Outbreak timing
  geom_vline(data = outbreak_metadata, 
             aes(xintercept = date, color = school, lty = school),
             lwd = ob_lwd) +
  # Sequence sampling dates
  geom_dotplot(data = phylo_metadata %>% 
                 filter(!hh) %>% 
                 mutate(what = "RT-PCR/antigenic"),
               aes(x = collection_date, color = school),
               binwidth = 1,
               fill = "white",
               dotsize = 1.3, method = "histodot", stackdir = "down",
               stroke = .5,
               stackgroups = TRUE) +
  facet_grid(what ~ ., switch = "y") +
  scale_color_manual(values = getSchoolColors()) +
  scale_linetype_manual(values = ob_lty) +
  theme_bw() +
  theme(strip.placement = "outside",
        panel.grid.minor.y = element_blank(),
        strip.switch.pad.grid = unit(1, units = "lines"))  +
  coord_cartesian(xlim = date_limits) +
  labs(x = "Date", y = "Weekly number of tests", 
       size = "# samples",
       color = "School/preschool", lty = "School/preschool",
       fill = "Test result") +
  scale_x_date(date_breaks = "2 months", date_labels = "%Y-%b") +
  scale_fill_manual(values = c( "lightgray", getVariantColors()))#c("lightgray", "black"))



## Assemble Figure 1 ----
p_fig1 <- plot_grid(
  p_context +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()),
  p_outbreak +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()),
  p_stats,
  align = "v",
  axis = "lr",
  labels = "auto",
  ncol = 1,
  rel_heights = c(.5, .55, 1)
)

# Save
ggsave(
  p_fig1, 
  filename = "figures/figure_1_context_data.png",
  width = 10,
  height = 8,
  dpi = 300
)


# Figure 2: phylogenetic trees --------------------------------------------


# Get metadata
phylo_meta <- getPhyloMetadata() %>% 
  left_join(households) %>% 
  # Put seq_id first to use to annotate ggtree
  select(seq_id, all_of(colnames(.))) %>% 
  addOutbreakID(date_col = "collection_date") %>% 
  group_by(seq_id) %>% 
  slice(1) %>% 
  mutate(variant = map_chr(collection_date, ~ getVariant(.)),
         school = factor(school, levels = getSchoolUMapping())) %>% 
  ungroup()  %>% 
  correctOutbreakID()

# Make zoom tree by outbreak
u_outbreaks <- unique(phylo_meta$outbreak_id[!is.na(phylo_meta$outbreak_id)])

all_variants <- c("alpha", "delta", "omicron")

variant_plots <- map(all_variants, function(x) {
  p <- plot_variant_tree(x, 
                         phylo_meta = phylo_meta)
  p
})

variant_plot_dims <- map_df(1:length(variant_plots), function(x) {
  p <- variant_plots[[x]]
  if (!is.null(p)) {
    tibble(
      plot = x,
      variant = all_variants[x],
      x_min = p$scales$scales[[2]]$limits[1],
      x_max = p$scales$scales[[2]]$limits[2],
      y_min = 0,
      y_max = 1
    )
  }
}) %>% 
  mutate(x_size = as.numeric(x_max-x_min),
         y_size = y_max-y_min)


walk(1:length(variant_plots), function(x) {
  p <- variant_plots[[x]] +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  if (x != 1) {
    p <- p +
      guides(fill = "none", color = "none")
  }
  
  if (!is.null(p)) {
    w_scale <- variant_plot_dims$x_size[variant_plot_dims$plot == x]/max(variant_plot_dims$x_size)
    this_width  <- 7 * variant_plot_dims$y_size[variant_plot_dims$plot == x]/max(variant_plot_dims$y_size)
    h_scale <- 1#variant_plot_dims$y_size[variant_plot_dims$plot == x]/max(variant_plot_dims$y_size)
    this_height <- 9 * h_scale
    
    ggsave(p, filename = str_glue("figures/tree_variant_{variant_plot_dims$variant[variant_plot_dims$plot == x]}.pdf"), 
           width = this_width, 
           height = this_height#, 
           # scale = 1/min(w_scale, h_scale)#,
           # dpi = 500
    )
  }
})

variant_plots_for_fig <- map(1:length(variant_plots), function(x) {
  if (x == 3) {
    variant_plots[[x]] +
      scale_color_manual(values = getSchoolColors()) +
      scale_fill_manual(values = getSchoolColors()) +
      scale_x_date(date_labels = "%Y-%b") +
      theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines"),
            title = element_text(hjust = 1),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  } else {
    variant_plots[[x]] +
      guides(color = "none", fill = "none")  +
      scale_x_date(date_labels = "%Y-%b")+
      theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines"),
            title = element_text(hjust = 1),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  }  
})

# Combine plots of variants
p_by_variant <- plot_grid(plotlist = variant_plots_for_fig,
                          nrow = 1, 
                          rel_widths = c(1.2, 2, 3.5),
                          labels = "auto")


outbreak_plots <- map(sort(u_outbreaks), function(x) {
  p <- try(plot_outbreak_tree(x, 
                              phylo_meta = phylo_meta,
                              what = "phylo"))
  
  if (!inherits(p, "try-error")) {
    p
  } else {
    NULL
  }
})

outbreak_plot_dims <- map_df(1:length(outbreak_plots), function(x) {
  p <- outbreak_plots[[x]]
  if (!is.null(p)) {
    tibble(
      plot = x,
      outbreak = sort(u_outbreaks)[x],
      x_min = p$coordinates$limits$x[1],
      x_max = p$coordinates$limits$x[2]
    )
  }
}) %>% 
  mutate(x_size = sqrt(as.numeric(x_max-x_min))#,
         # y_size = sqrt(y_max-y_min)
  )


walk(1:length(outbreak_plots), function(x) {
  p <- outbreak_plots[[x]]
  if (!is.null(p) & x != 4) {
    w_scale <- 1#outbreak_plot_dims$x_size[outbreak_plot_dims$plot == x]/max(outbreak_plot_dims$x_size)
    this_width  <- 6 * w_scale
    h_scale <- 1#outbreak_plot_dims$y_size[outbreak_plot_dims$plot == x]/max(outbreak_plot_dims$y_size)
    this_height <- 8 * h_scale
    
    ggsave(p, filename = str_glue("figures/tree_outbreak_{outbreak_plot_dims$outbreak[outbreak_plot_dims$plot == x]}.pdf"), 
           width = this_width, 
           height = this_height#, 
           # scale = 1/min(w_scale, h_scale)#,
           # dpi = 500
    )
  }
})


outbreak_plots_for_fig <- outbreak_plots[map_lgl(outbreak_plots, ~ !is.null(.))]
outbreak_variants <- c("alpha", rep("delta", 2), rep("omicron", 3))

outbreak_plots_by_variant <-  map(
  unique(outbreak_variants), 
  function(v) {
    plot_grid(plotlist = outbreak_plots_for_fig[map_lgl(outbreak_variants, ~ . == v)],
              nrow = 1) +
      theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines")) +
      ggtitle(v)
  })

p_outbreak_trees <- plot_grid(
  plotlist = outbreak_plots_by_variant,
  nrow = 1,
  labels = c("d", "e", "f"),
  rel_widths = table(outbreak_variants)
)

p_fig2 <- plot_grid(
  p_by_variant,
  p_outbreak_trees,
  ncol = 1,
  rel_heights = c(1.4, 1)
) +
  theme(plot.background = element_rect(fill = "white", color = "white"))


ggsave(p_fig2, 
       filename = "figures/fiture_2_new.png",
       width = 12,
       height = 9,
       dpi = 400)

# Figure 3: Inference results ---------------------------------------------

## Data wrangling for figure ----
# Model parameter estimates
param_estimates_forplot <- param_estimates %>% 
  filter(gen_like == "with_gen") %>% 
  computeProbFromParams() %>% 
  mutate(where_label = case_when(
    where == "community" ~ "community\n[/week]",
    T ~ "school\n[/infector/infectious period]"),
    what = factor(what, levels = c("community", "within-group", "between-group")
    )
  ) %>% 
  mutate(valid = !is.infinite(lo) & !is.infinite(hi)) %>% 
  filter(valid) 


# Function to make labels for panel c
makeWhatLabel <- function(df) {
  df %>% 
    mutate(what_label = case_when(what == "C" ~ "Cumulative exposures\n[-]", 
                                  what == "hazard" ~ "Infection hazard\n[/day]", 
                                  T ~ "Ratio of within-schoolto\ncommunity hazard\n[-]"))
}

pd <- position_dodge(width = .3)


## Panel a: community infectious pressure parameters ----
p_fig3a <- param_estimates_forplot %>% 
  filter(where == "community") %>% 
  ggplot(aes(x = variant, y = mean_prob, 
             ymin = lo_prob, ymax = hi_prob, color = variant,
             alpha = setting)) +
  geom_point(aes(pch = setting), position = pd, size = 2) +
  geom_errorbar(aes(group = setting),  lwd = .3,
                width = 0, position = pd)  +
  scale_y_log10() +
  coord_cartesian(ylim = c(1e-3, .3)) +
  scale_color_manual(values = getVariantColors()) +
  scale_alpha_manual(values = c(1, .5, .5)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        plot.margin = unit(c(1, 1, 1, 1), units = "lines"),
        legend.position = c(.8, .3),
        legend.key.width = unit(1, units = "lines"),
        legend.key.height = unit(.5, units = "lines"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8)) +
  labs(y = "Community infection hazard\n[probability per week]") +
  ggtitle("Community transmission\n") +
  guides(color = "none", pch = "none", alpha = "none")

## Panel b: community infectious pressure parameters ----
p_fig3b <- param_estimates_forplot %>% 
  filter(where != "community") %>% 
  ggplot(aes(x = variant, y = mean_prob, 
             ymin = lo_prob, ymax = hi_prob, color = variant,
             alpha = setting)) +
  geom_point(aes(pch = setting), position = pd) +
  geom_errorbar(aes(group = setting), lwd = .3,
                width = 0, position = pd)  +
  scale_y_log10() +
  coord_cartesian(ylim = c(1e-3, .3)) +
  scale_color_manual(values = getVariantColors()) +
  scale_alpha_manual(values = c(1, .5, .5)) +
  facet_grid(. ~ what) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()) +
  labs(y = "Within-school infection hazard\n[probability per infectious period]") +
  ggtitle("Within-school transmission") +
  theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines"))



# Combine parameter estimates
p_estimates_main2 <- plot_grid(
  p_fig3a,
  p_fig3b,
  labels = "auto",
  align = "h",
  axis = "tb",
  rel_widths = c(1, 1.75)
)


## Panel c: cumulative infections and hazard

# Cumulative infections
p_C <- traj_stats %>% 
  filter(age_cat == "child", what == "C") %>% 
  ggplot(aes(x = date, y = mean)) +
  geom_rect(data = vacations %>% 
              mutate(TL = TL + 2,
                     TR = TR + 3), 
            aes(xmin = TL, xmax = TR, ymin = 0, ymax = Inf),
            inherit.aes = F, 
            alpha = .10) +
  geom_vline(data = outbreak_metadata, 
             aes(xintercept = date, color = school, lty = school),
             lwd = .3) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = school), alpha = .15) +
  geom_line(aes(color = school), lwd = .5) +
  geom_hline(data = tribble(
    ~what, ~y,
    "C", 1,
    "hazard", NA,
    "ratio", 1
  ) %>% 
    makeWhatLabel(), 
  aes(yintercept = y), col = "darkgray", lwd = .3, lty = 1) +
  scale_fill_manual(values = getSchoolColors(), drop = F) +
  scale_color_manual(values = getSchoolColors(), drop = F) +
  scale_linetype_manual(values = ob_lty) +
  labs(x = "date", y = "Cumulative attack rate [-]") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  guides(fill = "none", color = "none", lty = "none")


# Figure of hazard with shading
hazard_data <-  traj_stats %>% 
  filter(age_cat == "child", what %in% c("hazard", "ratio")) %>% 
  select(what, school, date, mean) %>% 
  group_by(what, school, date) %>% 
  slice_max(mean, with_ties = F) %>% 
  ungroup() %>% 
  pivot_wider(values_from = c("mean"),
              names_from = c("what")) %>%
  mutate(hazard = hazard/100) %>% 
  mutate(hz_community = hazard/(1+ratio),
         hz_school = hazard - hz_community) %>% 
  filter(!(school == "SB" & date < "2021-08-01"))

p_hazard <- hazard_data %>% 
  ggplot(aes(x = date, fill = school)) +
  geom_rect(data = vacations %>%
              #   inner_join(tribble(~institution, ~school,
              #                      "school", "SA",
              #                      "school", "SB",
              #                      "preschool", "SC",
              #                      "preschool", "SD")) %>% 
              mutate(TL = TL + 2,
                     TR = TR + 3),
            aes(xmin = TL, xmax = TR, ymin = 0, ymax = Inf),
            inherit.aes = F, 
            alpha = .1) +
  geom_ribbon(aes(ymin = 0, ymax = hz_community), alpha = 1) +
  geom_ribbon(aes(ymin = hz_community, ymax = hazard), alpha = 1)  +
  # geom_ribbon(
  #   data = hazard_data %>% 
  #     group_by(date) %>% 
  #     slice_max(hazard, with_ties = FALSE) %>% 
  #     ungroup() %>% 
  #     filter(!is.na(date)),
  #   inherit.aes = FALSE,
  #   aes(x = date, ymin = hz_community, ymax = hazard), 
  #   alpha = .35, fill = "white") +
  ggpattern::geom_ribbon_pattern(
    data = hazard_data %>%
      group_by(date) %>%
      slice_max(hazard, with_ties = FALSE) %>%
      ungroup() %>%
      filter(!is.na(date)),
    inherit.aes = FALSE,
    aes(x = date, ymin = hz_community, ymax = hazard),
    alpha = 0.35, 
    fill = "white",
    pattern_fill = "white", 
    pattern_color = "white",
    pattern = "stripe", 
    pattern_spacing = .025,
    pattern_density = .02,
    pattern_angle = 45,
    pattern_alpha = .5)  +
  geom_segment(data = outbreak_metadata, 
               inherit.aes = FALSE,
               aes(x = date, xend = date, y = Inf, yend = 0, color = school, lty = school),
               lwd = .3) +
  # facet_wrap(~school, ncol = 4, dir = "v") +
  scale_linetype_manual(values = ob_lty) +
  theme_bw() +
  scale_fill_manual(values = getSchoolColors()) +
  scale_color_manual(values = getSchoolColors()) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  guides(fill = "none", color = "none", lty = "none") +
  labs(y = "infection hazard\n[probability per day]")

# p_hazard

# Incidence
# Population in the canton of Geneva: https://statistique.ge.ch/domaines/apercu.asp?dom=01_01
pop_gva <- 517802

rate_comp <- context_data %>% 
  filter(what == "weekly_cases") %>% 
  select(date, gva_cases = value) %>% 
  arrange(date) %>% 
  mutate(
    # Adjust for under-reporting
    gva_cases = gva_cases * 5,
    # Compute rates
    gva_rates = gva_cases/pop_gva*1e5
  ) %>% 
  addEpiWeek() %>% 
  group_by(epiweek_date) %>% 
  summarise(gva_rates = mean(gva_rates)) %>% 
  inner_join(
    traj_stats %>% 
      filter(age_cat == "child", what == "hazard") %>% 
      select(what, school, date, mean) %>% 
      group_by(what, school, date) %>% 
      slice_max(mean, with_ties = F) %>% 
      ungroup() %>% 
      mutate(school_rates = (1-(1-mean/100)^7) * 1e5) %>% 
      select(school, date, school_rates) %>% 
      addEpiWeek() %>% 
      group_by(school, epiweek_date) %>% 
      summarise(school_rates = mean(school_rates))
  ) %>% 
  mutate(variant = map_chr(epiweek_date, ~ getVariant(.))) %>% 
  filter(!is.na(variant))

rate_comp_by_variant <- rate_comp %>% 
  group_by(variant, school) %>% 
  summarise(school_rates = mean(school_rates),
            gva_rates = mean(gva_rates))

multipliers <- tibble(intercept = log10(c(.1, .2, .5, 1, 2, 5, 10, 50))) %>% 
  mutate(label = str_c("x ", formatC(10^(intercept), format = "f", digits = 1)) %>% 
           str_remove("\\.0"),
         label = factor(label) %>% 
           forcats::fct_reorder(intercept))

# https://stackoverflow.com/questions/50413812/log-axis-labels-in-ggplot2-show-only-necessary-digits
plain <- function(x,...) {
  formatC(x, ..., big.mark = ",", format = "f", drop0trailing = TRUE)
}

p_rate_comp <- rate_comp %>% 
  ggplot(aes(x = gva_rates, y = school_rates, color = variant)) +
  geom_abline(data = multipliers,
              aes(intercept = intercept, slope = 1), lty = 2, lwd = .3,
              color = "gray") +
  geom_abline(color = "darkgray") +
  ggtext::geom_richtext(data = multipliers, 
                        aes(x = 2e3, y = 2e3*10^(intercept), label = label),
                        angle = 45, 
                        inherit.aes = FALSE,
                        color = "darkgray",
                        label.colour = "white",
                        size = 3) +
  geom_point(alpha = 1, aes(pch = school), size = 1.25) +
  geom_path(alpha = .35, lwd = .5, aes(group = school)) +
  geom_point(data = rate_comp_by_variant,
             size = 6, pch = 23, fill = "white", alpha = .9) +
  geom_point(data = rate_comp_by_variant,
             size = 3,
             aes(pch = school), alpha = .9) +
  theme_bw() +
  scale_x_log10(labels = plain) +
  scale_y_log10(labels = plain) +
  scale_color_manual(values = getVariantColors()) +
  # scale_linetype_manual(values = c("dotted", "dotdash", "dashed", "solid", "dashed", "dotdash", "dotted", "dotted")) +
  coord_equal() +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Incidence rates in general population\n[infections per 100,000 population per week]",
       y = "Incidence rates in study primary and pre-schools\n[infections per 100,000 population per week]")


p_fig3_2 <- plot_grid(
  p_estimates_main2,
  plot_grid(
    plot_grid(
      p_C +
        theme(plot.margin = unit(c(1, 0, 0, 1), units = "lines")) +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank()) +
        ggtitle("Attack rates and hazard partition"),
      p_hazard +
        theme(plot.margin = unit(c(0, 0, 1, 1), units = "lines")),
      nrow = 2,
      labels = c("c", "d"),
      align = "v",
      axis = "lr"
    ),
    p_rate_comp +
      theme(panel.border = element_blank(),
            axis.line = element_line())  +
      ggtitle("School vs. general population incidence") +
      theme(plot.margin = unit(c(1, 0, 0, 0), units = "lines")),
    nrow = 1, 
    # align = "v",
    # axis = "tb",
    rel_widths = c(1, 1.65),
    labels = c(NA_character_, "e")
  ) +
    theme(plot.background = element_rect(color = "white", fill = "white")),
  labels = NULL,
  ncol = 1,
  rel_heights = c(1, 1.5),
  align = "v",
  axis = "lr"
)

ggsave(p_fig3_2, 
       filename = "figures/figure_3_inference_v2.png",
       width = 12,
       height = 10)


# Figure 4: Scenario simulations ------------------------------------------

# New attempt

varplots <- map(rev(c("Omicron-like", "Delta-like", "Alpha-like")), function(v) {
  col <- case_when(str_detect(v, "Omi")  ~ getVariantColors()[3],
                   str_detect(v, "Del")  ~ getVariantColors()[2],
                   T ~  getVariantColors()[1])
  
  effect_stats %>%
    filter(setting == "primary school", variant == v) %>% 
    ggplot() +
    geom_tile(aes(x = 100*(1-w), y = 100*(1-b), fill = 100*(1-median))) +
    scale_fill_gradient(low = "white", high = col, limits = c(0, 100)) +
    theme_bw()
}) 


varplots2 <- map(seq_along(varplots), function(x){
  if (x == 1) {
    varplots[[x]] +
      geom_segment(aes(x = 0, y = 0, xend = 99, yend = 99),
                   color = "black", lwd = .4,
                   arrow = arrow(length = unit(0.2, "cm"))) +
      geom_segment(aes(x = 0, y = 0, xend = 99, yend = 99),
                   color = "white", lwd = .2,
                   arrow = arrow(length = unit(0.2, "cm"))) +
      geom_point(data = tibble(
        x = seq(0, 100, by = 10),
        y = seq(0, 100, by = 10)
      ), #%>% 
      # filter(!(x ==0 & y == 0)),
      inherit.aes = FALSE,
      aes(x = x, y = y), fill = "white", color = "black", pch = 21, stroke = .3) +
      # guides(fill = "none") +
      labs(x = NULL, 
           y = "Between-class transmission\nreduction [%]") +
      # coord_equal() +
      theme(legend.position = "top",
            legend.title.position = "left", 
            legend.title = element_text(size = 10, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      # scale_x_continuous(position = "top") +
      guides(fill = guide_colorbar(title = "Averted    \ninfections [%]    ")) +
      ggtitle("Very rare introductions")
  } else if (x == 3) {
    varplots[[x]] +
      # guides(fill = "none") +
      labs(x = NULL, 
           y = NULL) +
      # coord_equal() +
      theme(legend.position = "top",
            legend.title.position = "top",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.background = element_blank(), 
            plot.title = element_text(hjust = 0.5)) +
      # scale_x_continuous(position = "bottom") +
      guides(fill = guide_colorbar(title = NULL)) +
      ggtitle("Very frequent introductions")
  }  else {
    varplots[[x]] +
      # guides(fill = "none") +
      labs(x = "Within-class transmission reduction [%]", 
           y = NULL) +
      # coord_equal() +
      theme(legend.position = "top",
            legend.title.position = "top",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            plot.title = element_text(hjust = 0.5)) +
      # scale_x_continuous(position = "bottom") +
      guides(fill = guide_colorbar(title = NULL)) +
      ggtitle("Rare introductions")
  }
})

p_grids <- plot_grid(plotlist = varplots2, 
                     nrow = 1,
                     rel_widths = c(1.15, 1, 1),
                     align = "h",
                     axis = "tblr")

diag_dat <- effect_stats %>% 
  filter(b == w, date == max(date), setting == "primary school") %>%
  bind_rows(tibble(w = 1, median = 1, scenario = levels(effect_stats$scenario)))

p_diag2 <- effect_stats %>% 
  filter(setting == "primary school", date == max(date)) %>%
  # bind_rows(tibble(w = 1, median = 0)) %>% 
  mutate(g = str_c(variant, b)) %>% 
  # bind_rows(tibble(w = 1, median = 0)) %>% 
  ggplot(aes(x = 100*(1-w), y = 100*(1-median))) +
  geom_line(aes(color = scenario, group = g), alpha = .4, lty = 2) +
  geom_line(data = diag_dat, aes(color = scenario)) +
  geom_point(data = diag_dat, aes(color = scenario), pch = 21, fill = "white",
             size = 2.5) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(),
        legend.title = element_text(size = 10, hjust = 1)) +
  scale_color_manual(values = getVariantColors()) +
  scale_x_continuous(breaks = seq(0, 100, by = 25)) +
  scale_y_continuous(breaks = seq(0, 100, by = 25)) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100))  +
  labs(y = "Averted infections [%]",
       x = "Within-class transmission reduction [%]",
       color = "Pathogen introductions\nfrom the community")


p_hist2 <- effect_stats %>% 
  filter(date == max(date), setting == "primary school") %>% 
  ggplot(aes(x = 100*(1-median))) +
  geom_histogram(aes(fill = scenario)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line()) +
  scale_fill_manual(values = getVariantColors()) +
  coord_flip() +
  # scale_y_continuous(position = "right") +
  # scale_x_continuous(position = "top") +
  labs(x = "Averted infections [%]", y = "# of scenarios")


p_scenarios_v2 <- plot_grid(
  p_grids +
    theme( plot.margin = unit(c(0, .5, 0, .5), units = "lines")),
  plot_grid(
    p_diag2 +
      theme(legend.position = "bottom",
            plot.margin = unit(c(1, 0, 1, 1.4), units = "lines")),
    p_hist2 + 
      guides(fill = "none") +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = unit(c(1, 1, 1, 1), units = "lines"),
            plot.background = element_blank()) +
      scale_x_continuous(breaks = seq(0, 100, by = 25)),#+
    # theme(plot.margin = unit(c(0, -4, 0, 0), units = "lines")),
    nrow = 1,
    align = "h",
    axis = "tb",
    rel_widths = c(2, 1),
    labels = c("b", "c")
  ),
  labels = c("a", NULL),
  ncol = 1,
  rel_heights = c(1, 1.5)
) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave(p_scenarios_v2, 
       filename = "figures/figure_4_new_v4.png",
       width = 9, height = 9.5, dpi = 300)


# Supplementary figures ---------------------------------------------------

## Prior sensitivity --------
prop_false_negs <- readRDS(here("generated_data/prop_false_negs.rds"))
pos_dat_stats <- readRDS(here("generated_data/pos_dat_stats.rds"))
X_sens_comb_df <- readRDS(here("generated_data/X_sens_comb_df.rds"))
dt_sens_pred <- readRDS(here("generated_data/dt_sens_pred.rds"))
genquant <- readRDS("generated_data/stan_model_genquant_non_pooled.rds")


T_max <- 500
t_log_raw <- log(1:T_max)
t_log_mean <- mean(t_log_raw)
t_log_sd <- sd(t_log_raw)

t_log_mat <- matrix(NA, nrow = length(dt_sens_pred), ncol = 4)

for(i in 1:length(dt_sens_pred)) {
  t_log_mat[i, 1] = 1;
  t_log_mat[i, 2] = (log(dt_sens_pred[i]) - t_log_mean)/t_log_sd;
  t_log_mat[i, 3] = t_log_mat[i, 2]^2;
  t_log_mat[i, 4] = t_log_mat[i, 2]^3;
}


inv_mat <- solve(t_log_mat[c(1, nrow(t_log_mat)), c(2, 4)])
beta <- matrix(NA, nrow = 4, ncol = 1)

sens_sims <- map_df(1:1000, function(i) {
  beta[1] = rnorm(1, 2, .5);
  beta[3] = rnorm(1, -1, .5);
  y1 = rnorm(1, -3.5, .75);
  y500 = rnorm(1, 0, .5);
  res <- inv_mat %*% matrix(
    c(y1 - t_log_mat[1, c(1, 3)] %*% beta[c(1, 3)],
      y500 - t_log_mat[nrow(t_log_mat), c(1, 3)] %*% beta[c(1, 3)]),
    ncol = 1
  )
  
  beta[c(2,4)] <- res
  
  # beta[4] = (y - t_log_mat[1, 1:3] %*% beta[1:3])/t_log_mat[1, 4];
  tibble(
    time = dt_sens_pred,
    sim = i,
    sens = as.numeric(inv_logit(t_log_mat %*% beta))
  )  
})


combined_sens <- bind_rows(
  sens_sims %>% 
    mutate(set = "prior"),
  genquant$draws("sens_pred") %>%
    as_draws() %>% 
    as_draws_df() %>% 
    as_tibble() %>% 
    pivot_longer(cols = contains("sens"),
                 names_to = "variable",
                 values_to = "sens") %>% 
    rename(sim = .draw) %>% 
    mutate(time = dt_sens_pred[as.numeric(str_extract(variable, "(?<=\\[)[0-9]+(?=,)"))],
           comb_row = as.numeric(str_extract(variable, "(?<=,)[0-9]+(?=\\])"))) %>%
    inner_join(X_sens_comb_df) %>% 
    mutate(set = "posterior")
) %>% 
  mutate(set = factor(set, levels = c("prior", "posterior")))


p_sens_traj <- combined_sens %>% 
  ggplot(aes(x = time, y = sens)) +
  geom_hex(aes(fill = after_stat(log(count)))) +
  geom_line(data = combined_sens %>% 
              filter(sim %in% sample(1:250, 50)),
            aes(group = sim),
            alpha = .2, color = "white",
            lwd = .35)+
  scale_fill_viridis_c() +
  facet_wrap(~set, ncol = 1) +
  theme_bw() +
  labs(x = "Time post infection [days]", y = "Sensitivity") +
  guides(fill = "none")


p_sens_data <- genquant$summary("sens_pred") %>%
  mutate(dt = dt_sens_pred[as.numeric(str_extract(variable, "(?<=\\[)[0-9]+(?=,)"))],
         comb_row = as.numeric(str_extract(variable, "(?<=,)[0-9]+(?=\\])"))) %>%
  inner_join(X_sens_comb_df) %>% 
  mutate(age_cat = "combined",
         group = "estimated this study") %>% 
  rename(prop = mean, lo = q5, hi = q95) %>% 
  # Add the emperical estimates for participants with test date
  bind_rows(
    prop_false_negs %>% 
      filter(what == "earliest_pos_date") %>% 
      mutate(group = "empirical this study") %>% 
      select(age_cat, group, dt = mean_delay, prop, lo, hi)
  ) %>% 
  # Add the control data from Michielin et al.
  bind_rows(
    pos_dat_stats
  ) %>% 
  mutate(group = factor(group, levels = c("Michielin et al.",
                                          "empirical this study",
                                          "estimated this study"))) %>% 
  ggplot(aes(x = dt, color = group)) +
  geom_point(alpha = .9, aes(pch = age_cat, y = prop, size = group)) +
  geom_errorbar(aes(ymin = lo, ymax = hi, group = age_cat), width = 0, alpha = .3) +
  theme_bw() +
  scale_color_manual(values = c("orange", "blue", "#707070")) +
  scale_size_manual(values = c(2.2, 2.2, 1)) +
  labs(x = "Time post infection [days]", y = "Sensitivity",
       shape = "Age class",
       color = "",
       size = "")


p_sens <- plot_grid(
  p_sens_traj,
  p_sens_data,
  labels = "auto",
  rel_widths = c(1, 2.1)
)

ggsave(p_sens, filename = "figures/supfig_sero_sens.png", 
       width = 11, height = 5.5, dpi = 300)



## Param estimates ----------------------

# Comparison by location
pd <- position_dodge(width = .2)
p_estimates <- param_estimates %>% 
  computeProbFromParams() %>% 
  mutate(valid = !is.infinite(lo) & !is.infinite(hi)) %>% 
  filter(valid) %>% 
  mutate(param = case_when(param == "lambda" ~ "community",
                           param == "lambdaA" ~ "school\nwithin-group",
                           T ~ "school\nbetween-group"),
         gen_like = case_when(gen_like == "with_gen" ~ "with genetic\nlikelihood",
                              T ~ "w/o genetic\nlikelihood")) %>% 
  ggplot(aes(x = variant, y = mean_prob, ymin = lo_prob, ymax = hi_prob, color = setting)) +
  geom_point(aes(pch = setting), position = pd) +
  geom_errorbar(width = 0, position = pd, alpha = .6) +
  # scale_color_manual(values = getVariantColors()) +
  facet_grid(param ~ gen_like) +
  # scale_y_log10() +
  theme_bw() +
  labs(x = "Period of\nSARS-CoV-2 variant", 
       y = "Probability of infection\n[community:weekly/school:per infectious period]")


ggsave(p_estimates, filename = "figures/all_param_estimates_supp.png", 
       width = 9, height = 5.5, dpi = 300)

p_estimates_v2 <- param_estimates %>%
  mutate(valid = !is.infinite(lo) & !is.infinite(hi)) %>%
  filter(valid) %>%
  ggplot(aes(x = variant, y = mle, ymin = lo, ymax = hi, color = gen_like)) +
  geom_point(aes(pch = gen_like), position = pd) +
  geom_errorbar(width = 0, position = pd, alpha = .6) +
  facet_grid(param ~ setting, scales = "free") +
  # scale_y_log10() +
  theme_bw()


ggsave(p_estimates_v2, filename = "figures/param_estimates_combined_by_setting.png",
       width = 10, height = 6, dpi = 300)


## Life-history parameters ---- 
nathist_estimates <- readRDS("generated_data/nathist_review_param.rds") %>% 
  filter(parameter != "serial interval",
         variant %in% c("alpha", "delta", "omicron"),
         str_detect(reference, "Manica", negate = T)) %>% 
  {
    x <- .
    bind_rows(x, 
              x %>%
                filter(str_detect(reference, "Har"), 
                       parameter == "generation time",
                       variant == "delta") %>% 
                mutate(estimate_mean = 0.9 * estimate_mean,
                       variant = "omicron")
    )
  } %>% 
  group_by(reference, variant, parameter) %>% 
  group_modify(function(x, y) {
    vals <- seq(0.01, 20, by = .1)
    gamma_dens <- dgamma(vals, shape = x$estimate_shape, scale = x$estimate_scale)
    
    tibble(vals = vals,
           dens = gamma_dens,
           what = "estimated") 
  }) %>% 
  ungroup()

nathist_erlang <- readRDS("generated_data/review_nathist_erlang.rds") %>% 
  filter(parameter == "incubation period",
         variant %in% c("alpha", "delta", "omicron"),
         str_detect(reference, "Man", negate = T))  %>% 
  group_by(reference, variant, parameter) %>% 
  group_modify(function(x, y) {
    vals <- seq(0.01, 20, by = .1)
    gamma_dens <- dgamma(vals, shape = x$erlang_shape, scale = x$erlang_scale)
    
    tibble(vals = vals,
           dens = gamma_dens,
           what = "erlang") 
  }) %>% 
  ungroup()

readRDS("generated_data/variant_nathis_erlang_param.rds") %>% 
  group_by(variant) %>% 
  group_modify(function(x, y) {
    vals <- seq(0.01, 20, by = .1)
    gamma_dens <- dgamma(vals, shape = x$k_I, scale = x$scale_I)
    tibble(vals = vals,
           dens = gamma_dens,
           what = "erlang",
           parameter = "symptomatic") 
  })

gentime_fit <- readRDS("generated_data/variant_nathis_erlang_fits.rds") %>% 
  group_by(variant) %>% 
  group_modify(function(x, y) {
    vals <- seq(0.01, 20, by = .1)
    scale <- x$v_g/x$e_g
    shape <- x$e_g/scale
    gamma_dens <- dgamma(vals, shape = shape, scale = scale)
    tibble(vals = vals,
           dens = gamma_dens,
           what = "fit",
           parameter = "generation time") 
  })
  

p_nathist <- bind_rows(nathist_erlang,
          nathist_estimates,
          gentime_fit) %>%  
  mutate(what = ifelse(what == "estimated", "Erlang approximation", "Estimate from literature")) %>% 
  ggplot(aes(x = vals, y = dens)) +
  geom_line(aes(color = variant, lty = what)) +
  facet_grid(parameter ~ .) +
  theme_bw() +
  scale_color_manual(values = getVariantColors()) +
  scale_linetype_manual(values = c(2, 1)) +
  labs(x = "days", y = "probability density",
       linetype = "Estimate")


ggsave(p_nathist, filename = "figures/supfig_natural_history_parameters.png",
       width = 7, height = 4.5, dpi = 300)


## Serologies ----

sero_long <- readRDS(here("generated_data/sero_long.rds")) %>% 
  mutate(event_type = map_chr(outbreak_id, ~ makeEventType(., simplify = F, add_terms = TRUE)),
         event_type = factor(event_type, levels = names(getEventTypeColors(with_terms = T))))   %>% 
  mutate(school = factor(school, levels = getAllSchools()),
         sugar_id = factor(sugar_id),
         sugar_id = forcats::fct_reorder2(sugar_id, school, group),
         school = factor(school, levels = getAllSchools() %>% sort())) %>% 
  filter(str_detect(event_type, "outbreak|baseline|start")) %>% 
  mutate(result_txt = ifelse(result == 0,  "negative", "positive"))

p_sero <- sero_long %>% 
  ggplot(aes(x = date, y = sugar_id)) +
  geom_line(aes(group = sugar_id), lwd = .2, alpha = .2) +
  geom_point(data = sero_long %>% filter(is.na(result)),
             aes(color = result_txt, alpha = !is.na(result))) +
  geom_point(data = sero_long %>% filter(result == 0),
             aes(color = result_txt, alpha = !is.na(result))) +
  geom_point(data = sero_long %>% filter(result == 1),
             aes(color = result_txt, alpha = !is.na(result))) +
  theme_bw() +
  # scale_color_manual(values = getEventTypeColors(with_terms = T)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_alpha_manual(values = c(.2, 1)) +
  facet_grid(school ~ ., scales = "free", space = "free_y") +
  labs(alpha = "Serology available", y = "participant",
       color = "Serology result") +
  scale_color_manual(values = c("blue", "red"))


ggsave(p_sero, filename = here("figures/supfig_all_serologies.png"), 
       width = 10, height = 10, dpi = 300)

## Tests ----

proc_test_data <- readRDS(here("generated_data/test_data_long.rds")) %>% 
  mutate(school = factor(school, levels = getAllSchools()),
         sugar_id = factor(sugar_id),
         sugar_id = forcats::fct_reorder2(sugar_id, school, group),
         outbreak_id = str_c("outbreak-", outbreak_id) %>% 
           factor(levels = str_c("outbreak-", 1:11))) %>% 
  arrange(term, group) %>% 
  mutate(source = case_when(source == "argos" ~ "ARGOS",
                            source == "visites" ~ "SEROCoV-Schools\n(visites)",
                            T ~ "SEROCoV-Schools\n(questionnaires)")) %>%
  filter(test_date < "2022-02-28") 

p_tests <- proc_test_data %>% 
  ggplot(aes(x = test_date, y = sugar_id)) +
  geom_line(aes(group = sugar_id), lwd = .2, alpha = .2) +
  geom_point(data = proc_test_data %>% filter(test_result == "negative"),
             aes(color = test_result)) +
  geom_point(data = proc_test_data %>% filter(test_result == "positive"),
             aes(color = test_result)) +
  facet_grid(school ~ ., scale = "free_y", space = "free_y") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(y = "participant", color = "RT-PCR/RADT result", x = "date") +
  scale_color_manual(values = c("blue", "red"))
  

ggsave(p_tests, filename = here("figures/supfig_all_tests.png"), 
       width = 10, height = 10, dpi = 300)

p_sero_tests <- plot_grid(
  p_sero +
    theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines"),
          legend.position = "bottom"), 
  p_tests +
    theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines"),
          legend.position = "bottom"),
  nrow = 1, 
  labels = "auto"
)


ggsave(p_sero_tests, filename = here("figures/supfig_all_tests_sero.png"), 
       width = 13, height = 9, dpi = 300)

## Statistical model fit ---
genquant <- readRDS(here("generated_data/stan_model_genquant_non_pooled.rds"))
event_comb <- readRDS(here("generated_data/event_comb.rds"))

p_stat_fit <- genquant$summary("gen_sero_event_obs") %>%
  bind_cols(event_comb, .) %>%
  mutate(event_type = str_replace(event_type, "debut", "start"),
         what = str_replace(what, "debut", "start"),
         what = case_when(str_detect(what, "base|start") ~ "baseline/start_school_year",
                          T ~ what),
         what = factor(what, levels = c("baseline/start_school_year", "outbreak", "last_visit")),
         event_type = factor(event_type, levels = names(getEventTypeColors(with_terms = TRUE))),
         group2 = str_c(outbreak_id, group)) %>% 
  filter(n_obs > 0) %>% 
  ggplot(aes(x = n_pos, group = event_id, color = event_type)) +
  geom_abline(lwd = .5, lty = 2) +
  geom_errorbar(aes(ymin = q5, ymax = q95), width = 0, alpha = .5, lwd = .5) +
  geom_line(aes(y = mean, group = group2), alpha = .2) +
  geom_point(aes(y = mean, group = group2, size = n_obs), alpha = .5) +
  scale_color_manual(values = getEventTypeColors(with_terms = TRUE)) +
  facet_wrap(what ~ ., ncol = 1) +
  labs(x = "# observed positive serologies",
       y = "# simulated positive serologies",
       color = "Event", 
       size = "# of serologies") +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_y_continuous(breaks = seq(0, 12, by = 2)) +
  theme_bw()

ggsave(p_stat_fit, filename = "figures/supfig_stat_fit.png",
       width = 7, height = 8, dpi = 300)

## Fidelity dynamic model ----

model_fidelity <- readRDS(file = "generated_data/model_fidelity.rds")

pd <- position_dodge(width = .05)

what_dict <- c("sero" = "serology",
               "pcr/antigen" = "RT-PCR/RADT")

p_fidelity <- model_fidelity %>% 
  mutate(outbreak_id = str_c("outbreak-", outbreak_id) %>% 
           factor(levels = str_c("outbreak-", 1:11)),
         what = what_dict[what]) %>% 
  ggplot(aes(x = n_pos, y = mean, group = date)) +
  geom_abline(lty = 2, lwd = .5) +
  geom_point(aes(size = n, color = outbreak_id), position = pd, alpha = .6) +
  geom_errorbar(aes(ymin = q025, ymax = q975, color = outbreak_id), width = 0, 
                position = pd, alpha = .4) +
  facet_grid(outbreak_id ~ what) +
  theme_bw() +
  scale_color_manual(values = getEventTypeColors()) +
  labs(x = "# observed positive samples", 
       y = "# simulated positive samples",
       color = "Outbreak",
       size = "# of samples")

ggsave(p_fidelity, filename = "figures/supfig_dyn_fit.png",
       width = 6, height = 12, dpi = 300)

## Stringency index ----
oxdata <- read_csv("data/03_study_info/OxCGRT_compact_national_v1.csv") %>% 
  janitor::clean_names()  %>% 
  filter(country_code == "CHE") %>% 
  mutate(date = as.Date(as.character(date), format = "%Y%m%d")) 

p_tile <- oxdata %>% 
  select(-contains("flag"), -country_name, -country_code, -region_name, -region_code, -jurisdiction,
         -contains("summary"), -contains("majority"), -contains("fiscal"), -contains("investment"), -contains("confirmed"), -contains("support"),
         -contains("index"), -contains("v2"), -population_vaccinated, -contains("debt")) %>% 
  pivot_longer(contains("_")) %>% 
  filter(date > "2021-01-01", date < "2022-02-01") %>% 
  ggplot(aes(x = date, y = name)) +
  geom_tile(aes(fill = as.factor(value))) +
  scale_fill_viridis_d() +
  theme_bw() +
  labs(y = "NPI", fill = "stringency score")

p_stringency <- oxdata %>% 
  select(date, contains("index"), -contains("economic"))  %>% 
  pivot_longer(contains("index")) %>% 
  filter(date > "2021-01-01", date < "2022-02-01") %>% 
  ggplot(aes(x = date, y = value)) +
  geom_step(aes(lty = name)) +
  theme_bw() +
  scale_linetype_manual(values = c(3, 2, 1)) +
  labs(y = "stringency index", lty = "category")

p_npis <- cowplot::plot_grid(
  p_stringency +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()),
  p_tile,
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = "auto",
  rel_heights = c(.4, 1)
)

ggsave(p_npis, filename = "figures/supfig_npi_stringency.png",
       width = 11, height = 8, dpi = 300)
