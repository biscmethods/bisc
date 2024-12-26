# script that runs all scripts in this folder

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
current_path <- file.path(demo_path, "simulated_data_to_stats_runseeds_pipeline")
output_path <- file.path(current_path, "output")
output_path_rds <- file.path(current_path, "rds")

# Check if the folder exists, if not, create it
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
if (!dir.exists(output_path_rds)) {
  dir.create(output_path_rds, recursive = TRUE)
}

set.seed(123)
redo_flag <- F

source(file.path(current_path, "step1_generate_data.R"))
source(file.path(current_path, "step2_run_biclustbiclust.R"))
source(file.path(current_path, "step3_run_bisc.R"))
source(file.path(current_path, "step4_run_biclustbiclust_many_times.R"))
source(file.path(current_path, "step5_run_bisc_many_times.R"))
source(file.path(current_path, "step6_boxplot.R"))
source(file.path(current_path, "step7_make_tables.R"))
source(file.path(current_path, "step8_plot_network_temp.R"))
