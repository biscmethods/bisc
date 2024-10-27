# !/usr/bin/Rscript

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

for(i in 1:length(bisc_results_list)){
  print(paste0('Now running outer iteration ', i))
  print(bisc_results_list[[i]]$RIs[[1]]$RI_cell_clustering_bisc)

  # cells[i] <- bisc_results_list[[i]]$RIs$RI_cell_clustering_bisc

  if(scenarios[[i]]$description=="Simple"){
    i_simple <- i_simple +1
    cells[[scenarios[[i]]$description]][i_simple] <- bisc_results_list[[i]]$RIs[[1]]$RI_cell_clustering_bisc
    biclust[[scenarios[[i]]$description]][i_simple] <- bisc_results_list[[i]]$RIs[[1]]$RI_biclust_bisc
    genes[[scenarios[[i]]$description]] <- c(genes[[scenarios[[i]]$description]], bisc_results_list[[i]]$RIs[[1]]$RI_gene_clustering_bisc)

    cells_biclust[[scenarios[[i]]$description]][i_simple] <- biclustbiclust_results_list[[i]]$RI_cell_clustering_biclustbiclust
    biclust_biclust[[scenarios[[i]]$description]][i_simple] <- biclustbiclust_results_list[[i]]$RI_biclust_biclustbiclust
    genes_biclust[[scenarios[[i]]$description]] <- c(genes_biclust[[scenarios[[i]]$description]],
                                                     as.numeric(strsplit(trimws(null2NA(biclustbiclust_results_list[[i]]$RI_gene_clustering_biclustbiclust_all)), " ")[[1]]))
  }else
  {
    i_complex <- i_complex +1
    cells[[scenarios[[i]]$description]][i_complex]   <- null2NA(bisc_results_list[[i]]$RIs[[1]]$RI_cell_clustering_bisc)
    biclust[[scenarios[[i]]$description]][i_complex] <- null2NA(bisc_results_list[[i]]$RIs[[1]]$RI_biclust_bisc)
    genes[[scenarios[[i]]$description]]              <- c(genes[[scenarios[[i]]$description]], bisc_results_list[[i]]$RIs[[1]]$RI_gene_clustering_bisc)

    cells_biclust[[scenarios[[i]]$description]][i_complex]   <- null2NA(biclustbiclust_results_list[[i]]$RI_cell_clustering_biclustbiclust)
    biclust_biclust[[scenarios[[i]]$description]][i_complex] <- null2NA(biclustbiclust_results_list[[i]]$RI_biclust_biclustbiclust)
    genes_biclust[[scenarios[[i]]$description]]              <- c(genes_biclust[[scenarios[[i]]$description]],
                                                                  as.numeric(strsplit(trimws( null2NA(biclustbiclust_results_list[[i]]$RI_gene_clustering_biclustbiclust_all)), " ")[[1]]))
  }

}



library(tidyverse)
library(hrbrthemes)
library(viridis)

# Create the data frame
cell_data <- bind_rows(
  tibble(method="BS", type = "Simple", value = unlist(cells$Simple)),
  tibble(method="BS", type = "Complex", value = unlist(cells$Complex)),
  tibble(method="BB", type = "Simple", value = unlist(cells_biclust$Simple)),
  tibble(method="BB", type = "Complex", value = unlist(cells_biclust$Complex))
)

biclust_data <- bind_rows(
  tibble(method="BS", type = "Simple", value = unlist(biclust$Simple)),
  tibble(method="BS", type = "Complex", value = unlist(biclust$Complex)),
  tibble(method="BB", type = "Simple", value = unlist(biclust_biclust$Simple)),
  tibble(method="BB", type = "Complex", value = unlist(biclust_biclust$Complex))
)

gene_data <- bind_rows(
  tibble(method="BS", type = "Simple", value = unlist(genes$Simple)),
  tibble(method="BS", type = "Complex", value = unlist(genes$Complex)),
  tibble(method="BB", type = "Simple", value = unlist(genes_biclust$Simple)),
  tibble(method="BB", type = "Complex", value = unlist(genes_biclust$Complex))
)


plot_height <- 400
plot_width <- 750

constructed_plot <- ggplot(cell_data, aes(x = interaction(method, type), y = value, fill = type)) +
  geom_violin(position = "dodge") +
  stat_summary(fun = median, geom = "crossbar", width = 0.6, color = "cyan") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black",  size=0.6, width = 0.3, alpha = 0.5, height = 0) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)
  ) +
  ggtitle("Cell clusters rand index comparison") +
  xlab("BB = biclustbiclust, BS = bisc") +
  ylab("RI") +
  ylim(0, 1)
print(constructed_plot)
png(file.path(output_path, paste0("boxplot_cell_cluster_RI_comparison.png")), width = plot_width, height = plot_height, units = "px")
print(constructed_plot)
dev.off()

constructed_plot <- ggplot(biclust_data, aes(x = interaction(method, type), y = value, fill = type)) +
  geom_violin(position = "dodge") +
  stat_summary(fun = median, geom = "crossbar", width = 0.6, color = "cyan") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black",  size=0.6, width = 0.3, alpha = 0.5, height = 0) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)
  ) +
  ggtitle("Biclust rand index comparison") +
  xlab("BB = biclustbiclust, BS = bisc") +
  ylab("RI") +
  ylim(0, 1)
print(constructed_plot)
png(file.path(output_path, paste0("boxplot_biclust_RI_comparison.png")), width = plot_width, height = plot_height, units = "px")
print(constructed_plot)
dev.off()

constructed_plot <- ggplot(gene_data, aes(x = interaction(method, type), y = value, fill = type)) +
  geom_violin(position = "dodge") +
  stat_summary(fun = median, geom = "crossbar", width = 0.6, color = "cyan") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black", size=0.6, width = 0.3, alpha = 0.5, height = 0) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)
  ) +
  ggtitle("Gene modules rand index comparison") +
  xlab("BB = biclustbiclust, BS = bisc") +
  ylab("RI") +
  ylim(0, 1)
print(constructed_plot)
png(file.path(output_path, paste0("boxplot_gene_modules_RI_comparison.png")), width = plot_width, height = plot_height, units = "px")
print(constructed_plot)
dev.off()
