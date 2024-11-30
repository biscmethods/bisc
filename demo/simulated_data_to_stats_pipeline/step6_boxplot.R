# !/usr/bin/Rscript
library(tidyverse)
library(hrbrthemes)
library(viridis)

string_counts <- table(unlist(lapply(scenarios, function(x) x$description)))
types <-  rownames(string_counts)
n_types <- length(types)

cells <- list()
cells$Simple <- c(length=unname(string_counts["Simple"]))
cells$Complex <- c(length=unname(string_counts["Complex"]))

biclust <- list()
biclust$Simple <- c(length=unname(string_counts["Simple"]))
biclust$Complex <- c(length=unname(string_counts["Complex"]))

genes <- list()
genes$Simple <- c()
genes$Complex <- c()


cells_biclust <- list()
cells_biclust$Simple <- c(length=unname(string_counts["Simple"]))
cells_biclust$Complex <- c(length=unname(string_counts["Complex"]))

biclust_biclust <- list()
biclust_biclust$Simple <- c(length=unname(string_counts["Simple"]))
biclust_biclust$Complex <- c(length=unname(string_counts["Complex"]))

genes_biclust <- list()
genes_biclust$Simple <- c()
genes_biclust$Complex <- c()

i_simple <- 0
i_complex <- 0

null2NA <- function(x){
  ifelse(is.null(x), NA, x)
}

for(i in seq_along(bisc_results_list)){
  print(paste0('Now running outer iteration ', i))
  print(bisc_results_list[[i]]$RIs[[1]]$RI_cell_clustering_bisc)

  for(i_seed in seq_along(bisc_results_list[[i]]$RIs)){
    if(length(bisc_results_list[[i]]$bisc_results[[i_seed]])>1 || !is.na(bisc_results_list[[i]]$bisc_results[[i_seed]])){
      if(bisc_results_list[[i]]$bisc_results[[i_seed]]$converged){
        new_cell_val <- null2NA(bisc_results_list[[i]]$RIs[[i_seed]]$RI_cell_clustering_bisc)
        new_biclust_val <- null2NA(bisc_results_list[[i]]$RIs[[i_seed]]$RI_biclust_bisc)
        new_genes_val <- bisc_results_list[[i]]$RIs[[i_seed]]$RI_gene_clustering_bisc
        cells[[scenarios[[i]]$description]] <- c(cells[[scenarios[[i]]$description]], new_cell_val)
        biclust[[scenarios[[i]]$description]] <- c(biclust[[scenarios[[i]]$description]], new_biclust_val)
        genes[[scenarios[[i]]$description]] <- c(genes[[scenarios[[i]]$description]], new_genes_val)
      }
      new_cell_biclust_val <- null2NA(biclustbiclust_results_list[[i]][[i_seed]]$RI_cell_clustering_biclustbiclust)
      new_biclust_biclust_val <- null2NA(biclustbiclust_results_list[[i]][[i_seed]]$RI_biclust_biclustbiclust)
      new_genes_biclust_val <- as.numeric(strsplit(trimws( null2NA(biclustbiclust_results_list[[i]][[i_seed]]$RI_gene_clustering_biclustbiclust_all)), " ")[[1]])
      cells_biclust[[scenarios[[i]]$description]] <- c(cells_biclust[[scenarios[[i]]$description]], new_cell_biclust_val)
      biclust_biclust[[scenarios[[i]]$description]] <- c(biclust_biclust[[scenarios[[i]]$description]], new_biclust_biclust_val)
      genes_biclust[[scenarios[[i]]$description]] <- c(genes_biclust[[scenarios[[i]]$description]], new_genes_biclust_val)
    }
  }
}

# Create the data frame
cell_data <- bind_rows(
  tibble(method="bisc", type = "Simple", value = unlist(cells$Simple)),
  tibble(method="bisc", type = "Complex", value = unlist(cells$Complex)),
  tibble(method="BCPlaid", type = "Simple", value = unlist(cells_biclust$Simple)),
  tibble(method="BCPlaid", type = "Complex", value = unlist(cells_biclust$Complex))
)

biclust_data <- bind_rows(
  tibble(method="bisc", type = "Simple", value = unlist(biclust$Simple)),
  tibble(method="bisc", type = "Complex", value = unlist(biclust$Complex)),
  tibble(method="BCPlaid", type = "Simple", value = unlist(biclust_biclust$Simple)),
  tibble(method="BCPlaid", type = "Complex", value = unlist(biclust_biclust$Complex))
)

gene_data <- bind_rows(
  tibble(method="bisc", type = "Simple", value = unlist(genes$Simple)),
  tibble(method="bisc", type = "Complex", value = unlist(genes$Complex)),
  tibble(method="BCPlaid", type = "Simple", value = unlist(genes_biclust$Simple)),
  tibble(method="BCPlaid", type = "Complex", value = unlist(genes_biclust$Complex))
)

constructed_plot <- ggplot(cell_data, aes(x = interaction(method, type), y = value, fill = type)) +
  geom_violin(position = "dodge") +
  stat_summary(fun = median, geom = "crossbar", width = 0.6, color = "cyan") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black",  size=0.6, width = 0.3, alpha = 0.5, height = 0) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  ggtitle("Cell clusters rand index comparison") +
  xlab("Method.ScenarioType") +
  ylab("RI") +
  ylim(0, 1)
print(constructed_plot)
ggplot2::ggsave(file.path(output_path,  paste0("boxplot_cell_cluster_RI_comparison.pdf")), constructed_plot, width = 8, height = 4)


constructed_plot <- ggplot(biclust_data, aes(x = interaction(method, type), y = value, fill = type)) +
  geom_violin(position = "dodge") +
  stat_summary(fun = median, geom = "crossbar", width = 0.6, color = "cyan") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black",  size=0.6, width = 0.3, alpha = 0.5, height = 0) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  ggtitle("Biclust rand index comparison") +
  xlab("Method.ScenarioType") +
  ylab("RI") +
  ylim(0, 1)
print(constructed_plot)
ggplot2::ggsave(file.path(output_path,  paste0("boxplot_biclust_RI_comparison.pdf")), constructed_plot, width = 8, height = 4)

constructed_plot <- ggplot(gene_data, aes(x = interaction(method, type), y = value, fill = type)) +
  geom_violin(position = "dodge") +
  stat_summary(fun = median, geom = "crossbar", width = 0.6, color = "cyan") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black", size=0.6, width = 0.3, alpha = 0.5, height = 0) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  ggtitle("Gene modules rand index comparison") +
  xlab("Method.ScenarioType") +
  ylab("RI") +
  ylim(0, 1)
print(constructed_plot)
ggplot2::ggsave(file.path(output_path,  paste0("boxplot_gene_modules_RI_comparison.pdf")), constructed_plot, width = 8, height = 4)
