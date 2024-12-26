# !/usr/bin/Rscript

library(ggplot2) # To access a color palette
library(scales)  # For additional palettes
library(dplyr)
# df_si <- data.frame(
#   simple = logical(),
#   n_cell_clusters = integer(),
#   n_target_genes = integer(),
#   n_regulator_genes = integer(),
#   n_total_cells = integer(),
#   stringsAsFactors = FALSE
# )
#
#
# for(i_scenario in seq_along(scenarios)){
#   cs <- scenarios[[i_scenario]]
#   new_row <- data.frame(simple = cs$description == "Simple",  # chr "Simple"
#                         n_cell_clusters = cs$n_cell_clusters,  # 2
#                         n_target_gene_clusters = cs$n_target_gene_clusters,  # num [1:2] 4 2
#                         n_target_genes = cs$n_target_genes,  # 100
#                         n_regulator_genes = cs$n_regulator_genes,  # 6
#                         n_cells = cs$n_cells,  # int [1:2] 200 202
#                         n_total_cells = cs$n_total_cells  # 398
#   )
#
#   df_si <- rbind(df_si, new_row)
# }

# plotta antalet celler i ett genkluster mot converged?

# !/usr/bin/Rscript




# Plot Rand index vs BIC ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_mp <- data.frame(
  i_scenario = integer(),
  RI_cellclusters = numeric(),
  RI_genemodules = numeric(),
  RI_biclusters = numeric(),
  n_iterations = integer(),
  converged = logical(),
  silhouette_measure = numeric(),
  davies_bouldin_index = numeric(),
  BIC = numeric(),
  simple = logical(),
  n_cell_clusters = integer(),
  n_target_genes = integer(),
  n_regulator_genes = integer(),
  n_total_cells = integer(),
  cell_targetgene_ratio = numeric(),
  n_modules_cells_ratio = numeric(),
  stringsAsFactors = FALSE
)
for(i_scenario_results in seq_along(bisc_results_list)){
  cs <- scenarios[[i_scenario_results]]
  new_row <- data.frame(simple = cs$description == "Simple",  # chr "Simple"
                        n_cell_clusters = cs$n_cell_clusters,  # 2
                        n_target_gene_clusters = cs$n_target_gene_clusters,  # num [1:2] 4 2
                        n_target_genes = cs$n_target_genes,  # 100
                        n_regulator_genes = cs$n_regulator_genes,  # 6
                        # n_cells = cs$n_cells,  # int [1:2] 200 202
                        n_total_cells = cs$n_total_cells,  # 398
                        cell_targetgene_ratio = cs$n_total_cells/cs$n_target_genes,
                        n_modules_cells_ratio = cs$n_cells/cs$n_target_gene_clusters
  )

  for(i_current_results in seq_along(bisc_results_list[[i_scenario_results]]$bisc_results)){
    print(i_scenario_results)
    current_results <- bisc_results_list[[i_scenario_results]]$bisc_results[[i_current_results]]
    current_RI <- bisc_results_list[[i_scenario_results]]$RI[[i_current_results]]

    if(length(current_results) > 1 || !is.na(current_results)){
      new_df <- data.frame(i_scenario = i_scenario_results,
                           RI_cellclusters = current_RI$RI_cell_clustering_bisc,
                           RI_genemodules = mean(current_RI$RI_gene_clustering_bisc),
                           RI_biclusters = current_RI$RI_biclust_bisc,
                           n_iterations = current_results$n_iterations,
                           converged = current_results$converged,
                           silhouette_measure = current_results$silhouette_measure,
                           davies_bouldin_index =  current_results$davies_bouldin_index,
                           BIC = Rmpfr::asNumeric(current_results$BIC[[current_results$n_iterations]])
      )
      df_mp <- rbind(df_mp, cbind(new_df,new_row))
    }
  }
}
# df_mp <- df_mp[df_mp$converged,]
# df_mp <- df_mp %>%
#   group_by(i_scenario) %>%
#   mutate(BIC = (BIC - min(BIC)) / (max(BIC) - min(BIC))) %>%
#   ungroup() %>%
#   # Convert i_scenario to factor with correct ordering
#   mutate(i_scenario = factor(i_scenario, levels = sort(unique(i_scenario))))

df_means <- df_mp %>%
  group_by(i_scenario) %>%
  summarise(across(everything(), mean, na.rm = TRUE))


panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

upper.panel<-function(x, y){
    points(x,y)
    abline(lm(y~x), col='red')
}

pairs(df_means[,c(2:6,9,11:17)],
      lower.panel = panel.cor,
      upper.panel = upper.panel)



# df_mp <- df_mp[df_mp$converged,]
# df_mp <- df_mp[,c(1,2,4,7)]
# df_mp <- df_mp %>%
#   group_by(i_scenario) %>%
#   mutate(BIC = (BIC - min(BIC)) / (max(BIC) - min(BIC))) %>%
#   ungroup() %>%
#   # Convert i_scenario to factor with correct ordering
#   mutate(i_scenario = factor(i_scenario, levels = sort(unique(i_scenario))))

p <- ggplot(df_mp, aes(x = n_regulator_genes, y = rand_index, color = converged)) +
  geom_point(size = 4, alpha = 0.5) +
  facet_wrap(~ i_scenario, scales = "free", labeller = as_labeller(function(x) {
    ifelse(as.numeric(x) <= 10,
           paste("Simple scenario", x),
           paste("Complex scenario", as.numeric(x) - 10))
  })) +
  labs(
    title = "Rand Index vs Normalised BIC by Scenario",
    x = "Normalised BIC (within each scenario)",
    y = "Rand Index"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"),
                     name = "bisc converge")
print(p)
ggsave(file.path(output_path, paste0("RI_vs_BIC_per_scenario.pdf")), p, width = 12, height = 8)

