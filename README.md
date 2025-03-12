# Evolving infectious disease dynamics shape school-based intervention effectiveness

This repository contains code for the analysis of the SeroCoV-Schools study presented in the manuscript:

Perez-Saez J.\*, Bellon M.\*, Lessler J., Berthelot J., Hodcroft E., Michielin G., , Pennacchio F, Lamour J., Laubscher F.,L’Huillier A. G., Posfay-Barbe K. M., Maerkl S. J., Guessous I., Azman A. S., Eckerle I.\#, Stringhini S.\#, Lorthe E.\#; for the SEROCoV-Schools study group, *Evolving infectious disease dynamics shape school-based intervention effectiveness*.

## Data protection

Data in this study is protected by the ethics committee authorization of the research project. 
When necessary code used in data processing for the production of the analysis data have been masked.
Synthetic data to run the analysis pipeline will be made available.

## Repository structure
Code to run the analysis steps are located in folder `analysis` and numbered in the order they are intended to be run.

- Step 1: Inference with the statistical modeling framework (`01_run_statistical_model.R`)
- Setp 2: Inference wth the dynamical modeling framework
  - 2.1: Computation of natural history parameters by VOC using published estimates (`02_compute_variant_nathist_params.R`)
  - 2.2: Prepration of data for particle filtering using the `spatPomp` package (`03_prepare_data_for_spatpomp.R`)
  - 2.3: Run particle filter to get MLE (`04_run_param_search.R`) and confidence intervals with profiling (`05_run_param_profile.R`)
  - 2.4: Get epidemic trajectories from smoothing distribution (`06_draw_from_smoothing_distribution.R`) and compute summary statistics (`07_compute_smoothing_stats.R`)
- Step 3: Scenario simulations
  - 3.1: Make configs to run scenario simulations (`08_make_scenario_simulations.R`)
  - 3.2: Run simulations (`09_run_scenario_simulations.R`)
- Step 4: Making figures (`10_make_final_figures.R`) and stats (`11_extract_param_intervals.R`) for manuscript

