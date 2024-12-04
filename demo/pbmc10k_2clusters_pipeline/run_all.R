#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
while (sink.number() > 0) {
  sink()
  sink(file = NULL)
}

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
current_path <- file.path(demo_path, "pbmc10k_2clusters_pipeline")
output_path <- file.path(current_path, "output")
local_data <- file.path(current_path, "data")

# Check if the folder exists, if not, create it
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
if (!dir.exists(local_data)) {
  dir.create(local_data, recursive = TRUE)
}

set.seed(123)
redo_flag <- F

source(file.path(current_path, "step1_clean_data.R"))
source(file.path(current_path, "step2_run_scregclust.R"))
source(file.path(current_path, "step3_analysis_of_scregclust_runs.R"))
source(file.path(current_path, "step4_run_bisc.R"))
source(file.path(current_path, "step5_run_bcplaid.R"))
source(file.path(current_path, "step6_make_tables.R"))
