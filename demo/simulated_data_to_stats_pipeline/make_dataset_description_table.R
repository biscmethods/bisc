format_dataset_description <- function(scenario_list, description, n_cell_clusters, n_target_gene_clusters, n_target_genes, n_regulator_genes, n_cells, regulator_means, coefficient_means, coefficient_sds, disturbed_fraction, testing_penalization_data_gen, trivial_regulator_networks, seed, plot_heatmap = TRUE, plot_suffix = "simple-scenario-1",num = 1) {
  output <- sprintf("Scenario %i & %s & %d & %s & %d & %d & %s & %s & %s \\\\",
                    num,
                    description,
                    n_cell_clusters,
                    paste(as.character(n_target_gene_clusters), collapse = ", "),
                    n_target_genes,
                    n_regulator_genes,
                    # sum(n_cells),
                    paste0("(",paste0(n_cells, collapse = ", "), ")"),
                    paste(
                      lapply(coefficient_means, function(x) {
                        sprintf("(%s)", paste(as.character(x), collapse = ", "))
                      }),
                      collapse = ", "
                    ),
                    paste(
                      lapply(coefficient_sds, function(x) {
                        sprintf("(%s)", paste(as.character(x), collapse = ", "))
                      }),
                      collapse = ", "
                    )
                    )

  cat(output, "\n")
}

sink(paste0(output_path, "/", "scenario_table_info"))

format_dataset_description(scenario_list = scenarios,
                description = "Simple",
                n_cell_clusters = 2,
                n_target_gene_clusters = c(4, 2),  # Number of target gene clusters in each cell cluster
                n_target_genes = 100,
                n_regulator_genes = 6,
                n_cells = c(200, 200)*2,
                regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                coefficient_means = list(c(1, 3, 5, 7), c(10, 20)),  # For generating dummy data, coefficient means in each cell cluster
                coefficient_sds = list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01)),
                disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                testing_penalization_data_gen = c(0.1, 0.5),
                trivial_regulator_networks = TRUE,
                seed = 1234,
                plot_heatmap=TRUE,
                plot_suffix="simple-scenario-1",
                num = 1)


format_dataset_description(scenario_list = scenarios,
                             description = "Simple",
                             n_cell_clusters = 3,
                             n_target_gene_clusters = c(4, 2, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 120,
                             n_regulator_genes = 8,
                             n_cells = c(220, 220, 220)*2,
                             regulator_means = c(0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3, 5, 7), c(10, 12), c(22, 24)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01)),
                             disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5, 0.5),
                             trivial_regulator_networks = TRUE,
                             seed = 1234,
                             num = 2)

format_dataset_description(scenario_list = scenarios,
                             description = "Simple",
                             n_cell_clusters = 4,
                             n_target_gene_clusters = c(4, 2, 2, 3),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 120,
                             n_regulator_genes = 11,
                             n_cells = c(220, 220, 220, 220)*2,
                             regulator_means = c(0, 0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3, 5, 7), c(10, 12), c(22, 24), c(33, 35, 37)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01, 0.01)),
                             disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5, 0.5, 0.5),
                             trivial_regulator_networks = TRUE,
                             seed = 1234,
                             num = 3)

format_dataset_description(scenario_list = scenarios,
                             description = "Simple",
                             n_cell_clusters = 4,
                             n_target_gene_clusters = c(4, 3, 2, 3),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 160,
                             n_regulator_genes = 12,
                             n_cells = c(220, 220, 220, 220)*2,
                             regulator_means = c(0, 0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3, 5, 7), c(15, 17, 13), c(22, 24), c(33, 35, 37)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01, 0.01)),
                             disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5, 0.5, 0.5),
                             trivial_regulator_networks = TRUE,
                             seed = 1234,
                             num = 4)

format_dataset_description(scenario_list = scenarios,
                             description = "Simple",
                             n_cell_clusters = 4,
                             n_target_gene_clusters = c(2, 2, 2, 5),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 240,
                             n_regulator_genes = 11,
                             n_cells = c(250, 250, 250, 250)*2,
                             regulator_means = c(0, 0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3), c(15, 17), c(22, 24), c(33, 35, 37, 39, 41)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01), c(0.01, 0.01, 0.01, 0.01, 0.01)),
                             disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5, 0.5, 0.5),
                             trivial_regulator_networks = TRUE,
                             seed = 1234,
                             num = 5)

format_dataset_description(scenario_list = scenarios,
                             description = "Simple",
                             n_cell_clusters = 2,
                             n_target_gene_clusters = c(4, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 100,
                             n_regulator_genes = 6,
                             n_cells = c(200, 200)*2,
                             regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3, 5, 7), c(10, 20)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.06, 0.03, 0.04, 0.02), c(0.05, 0.02)),
                             disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.1, 0.5),
                             trivial_regulator_networks = TRUE,
                             seed = 1234,
                             num = 6)

format_dataset_description(scenario_list = scenarios,
                             description = "Simple",
                             n_cell_clusters = 2,
                             n_target_gene_clusters = c(4, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 100,
                             n_regulator_genes = 6,
                             n_cells = c(400, 200)*2,
                             regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3, 5, 7), c(10, 20)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.1, 0.03, 0.04, 0.02), c(0.05, 0.02)),
                             disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.1, 0.5),
                             trivial_regulator_networks = TRUE,
                             seed = 1234,
                             num = 7)

format_dataset_description(scenario_list = scenarios,
                             description = "Simple",
                             n_cell_clusters = 2,
                             n_target_gene_clusters = c(7, 3),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 100,
                             n_regulator_genes = 10,
                             n_cells = c(600, 300)*2,
                             regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3, 5, 7, 9, 11, 13), c(22, 20, 24)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.1, 0.03, 0.04, 0.02, 0.01, 0.01, 0.01), c(0.05, 0.02, 0.01)),
                             disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.1, 0.5),
                             trivial_regulator_networks = TRUE,
                             seed = 1234,
                             num = 8)

format_dataset_description(scenario_list = scenarios,
                             description = "Simple",
                             n_cell_clusters = 3,
                             n_target_gene_clusters = c(2, 2, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 120,
                             n_regulator_genes = 6,
                             n_cells = c(300, 220, 260)*2,
                             regulator_means = c(0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3), c(10, 12), c(22, 24)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.02, 0.06), c(0.04, 0.03), c(0.02, 0.04)),
                             disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5, 0.5),
                             trivial_regulator_networks = TRUE,
                             seed = 1234,
                             num = 9)

format_dataset_description(scenario_list = scenarios,
                             description = "Simple",
                             n_cell_clusters = 3,
                             n_target_gene_clusters = c(2, 10, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 220,
                             n_regulator_genes = 14,
                             n_cells = c(300, 2200, 260)*2,
                             regulator_means = c(0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3), c(10, 12, 14, 16, 18, 20, 22, 24, 26, 28), c(51, 53)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.02, 0.06), c(0.04, 0.03, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01), c(0.02, 0.04)),
                             disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5, 0.5),
                             trivial_regulator_networks = TRUE,
                             seed = 1234,
                             num = 10)



# Generate multiple complicated examples --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

format_dataset_description(scenario_list = scenarios,
                             description = "Complex",
                             n_cell_clusters = 2,
                             n_target_gene_clusters = c(3, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 300,
                             n_regulator_genes = 20,
                             n_cells = c(600, 400)*2,
                             regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 2, 3), c(18, 20)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.01, 0.01, 0.01), c(0.03, 0.01)),
                             disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5),
                             trivial_regulator_networks = FALSE,
                             seed = 1234,
                             num = 11)

format_dataset_description(scenario_list = scenarios,
                             description = "Complex",
                             n_cell_clusters = 2,
                             n_target_gene_clusters = c(4, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 300,
                             n_regulator_genes = 20,
                             n_cells = c(500, 300)*2,
                             regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3, 5, 7), c(18, 20)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.01, 0.01, 0.01, 0.01), c(0.01, 0.01)),
                             disturbed_fraction = 0.1,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5),
                             trivial_regulator_networks = FALSE,
                             seed = 123,
                             num = 12)

format_dataset_description(scenario_list = scenarios,
                             description = "Complex",
                             n_cell_clusters = 2,
                             n_target_gene_clusters = c(3, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 300,
                             n_regulator_genes = 20,
                             n_cells = c(600, 400)*2,
                             regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 2, 3), c(18, 20)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.01, 0.01, 0.01), c(0.03, 0.01)),
                             disturbed_fraction = 0.3,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5),
                             trivial_regulator_networks = FALSE,
                             seed = 124,
                             num = 13)

format_dataset_description(scenario_list = scenarios,
                             description = "Complex",
                             n_cell_clusters = 3,
                             n_target_gene_clusters = c(3, 3, 3),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 300,
                             n_regulator_genes = 20,
                             n_cells = c(500, 300, 500)*2,
                             regulator_means = c(0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 2, 3), c(22, 20, 24), c(10,12,14)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.03, 0.01, 0.03), c(0.02, 0.04, 0.01), c(0.02, 0.04, 0.01)),
                             disturbed_fraction = 0.5,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5, 0.5),
                             trivial_regulator_networks = FALSE,
                             seed = 1236,
                             num = 14)

format_dataset_description(scenario_list = scenarios,
                             description = "Complex",
                             n_cell_clusters = 4,
                             n_target_gene_clusters = c(3, 3, 3, 5),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 500,
                             n_regulator_genes = 30,
                             n_cells = c(500, 300, 500, 300)*2,
                             regulator_means = c(0, 0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 2, 3), c(22, 20, 24), c(10,12,14), c(40,42,44,46,48)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.03, 0.01, 0.03), c(0.02, 0.04, 0.01), c(0.02, 0.04, 0.01), c(0.01, 0.02, 0.012, 0.011, 0.005)),
                             disturbed_fraction = 0.101,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5, 0.5, 0.5),
                             trivial_regulator_networks = FALSE,
                             seed = 126,
                             num = 15)

format_dataset_description(scenario_list = scenarios,
                             description = "Complex",
                             n_cell_clusters = 2,
                             n_target_gene_clusters = c(2, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 100,
                             n_regulator_genes = 8,
                             n_cells = c(340, 340)*2,
                             regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(5, 3), c(18, 20)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.01, 0.01), c(0.03, 0.05)),
                             disturbed_fraction = 0.23,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5),
                             trivial_regulator_networks = FALSE,
                             seed = 1224,
                             num = 16)

format_dataset_description(scenario_list = scenarios,
                             description = "Complex",
                             n_cell_clusters = 2,
                             n_target_gene_clusters = c(4, 4),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 120,
                             n_regulator_genes = 12,
                             n_cells = c(440, 440)*2,
                             regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(15, 13, 14, 10), c(31, 30, 28, 29)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.03, 0.012, 0.011, 0.015), c(0.03, 0.05, 0.01, 0.01)),
                             disturbed_fraction = 0.33,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5),
                             trivial_regulator_networks = FALSE,
                             seed = 124,
                             num = 17)

format_dataset_description(scenario_list = scenarios,
                             description = "Complex",
                             n_cell_clusters = 2,
                             n_target_gene_clusters = c(2, 3),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 130,
                             n_regulator_genes = 15,
                             n_cells = c(500, 500)*2,
                             regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(15, 13), c(30, 28, 29)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.03, 0.012), c(0.03, 0.05, 0.01)),
                             disturbed_fraction = 0.23,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5),
                             trivial_regulator_networks = FALSE,
                             seed = 1214,
                             num = 18)

format_dataset_description(scenario_list = scenarios,
                             description = "Complex",
                             n_cell_clusters = 4,
                             n_target_gene_clusters = c(2, 3, 2, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 130,
                             n_regulator_genes = 15,
                             n_cells = c(500, 500, 400, 400)*2,
                             regulator_means = c(0, 0, 0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(15, 13), c(30, 28, 29), c(43, 45), c(20,22)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.03, 0.012), c(0.03, 0.05, 0.01), c(0.01, 0.01), c(0.01, 0.01)),
                             disturbed_fraction = 0.13,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5, 0.5, 0.5),
                             trivial_regulator_networks = FALSE,
                             seed = 1241,
                             num = 19)

format_dataset_description(scenario_list = scenarios,
                             description = "Complex",
                             n_cell_clusters = 2,
                             n_target_gene_clusters = c(4, 2),  # Number of target gene clusters in each cell cluster
                             n_target_genes = 100,
                             n_regulator_genes = 8,
                             n_cells = c(220, 220)*2,
                             regulator_means = c(0, 0), # For generating dummy data, regulator mean in each cell cluster
                             coefficient_means = list(c(1, 3, 5, 7), c(25, 27)),  # For generating dummy data, coefficient means in each cell cluster
                             coefficient_sds = list(c(0.03, 0.02, 0.01, 0.01), c(0.01, 0.01)),
                             disturbed_fraction = 0.56,  # Value between 0 and 1. How large portion of cells should move to other cell clusters.
                             testing_penalization_data_gen = c(0.5, 0.5),
                             trivial_regulator_networks = FALSE,
                             seed = 631,
                             num = 20)
while (sink.number() > 0) {
  sink()
  sink(file = NULL)
}

