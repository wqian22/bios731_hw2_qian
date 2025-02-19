# bios731_hw2_qian

## Project Directory

This project runs simulations for BIOS 731 Homework 2.


## Workflow

1. Run the file `simulations/run_simulations.R`. This file runs simulations and sources the following files:
  * `source/01_simulate_data.R`: contains function used to simulate data
  * `source/02_apply_methods.R`: contains function for calculating negative log-likelihood, applying Newton method and MM-algorithm
Running this file outputs the simulation results to the `results` folder.

2. Run the Rmarkdown file `analysis/HW2_report.Rmd`. This file will pull together the simulation results in `results` folder and generated some tables and plots.
