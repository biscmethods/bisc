# script that runs all scripts in this folder

#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
sink()

# options(warn=2)  # To convert warning messages into error messages which display row of error. For debugging.

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- here::here("demo/simulated_data_to_stats_pipeline")
path_data <- here::here('data')


redo_flag <- FALSE

source(file.path(output_path, "step1_generate_data.R"))
source(file.path(output_path, "step2_run_biclustbiclust.R"))
source(file.path(output_path, "step3_run_biclust_screg.R"))
source(file.path(output_path, "step4_run_biclustbiclust_many_times.R"))
source(file.path(output_path, "step5_run_biclust_screg_many_times.R"))
