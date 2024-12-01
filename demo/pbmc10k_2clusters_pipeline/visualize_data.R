#!/usr/bin/Rscript
rm(list = ls())

library(here)  # To work with paths
library(Seurat)
library(scregclust)
library(data.table)
# BiocManager::install('GenomicRanges', force=TRUE)
library("glmGamPoi")

# Get absolute path where script is located, by using relative paths.
demo_path <- here::here("demo")
R_path <- here::here("R")
output_path <- file.path(demo_path,'pbmc10k_2clusters_pipeline', 'output')
path_data <- here::here('data')

path_pbmc10k_sctransform <- file.path(path_data, "env_data_pbmc10k_sctransform.RData")
path_labels <- file.path(path_data, "GEX Graph-Based10k.csv")

# create output folder
if (!dir.exists(output_path)) {
  # If the folder doesn't exist, create it
  dir.create(output_path, recursive = TRUE)
  cat("Folder created:", output_path, "\n")
} else {
  cat("Folder already exists:", output_path, "\n")
}


plot_height <- 400
plot_width <- 750

load(path_pbmc10k_sctransform)
str(d)
d <- t(d)

labels <- read.csv(path_labels)


cells <- cbind(Barcode  = rownames(d), cellnum = 1:nrow(d))


key <- merge(labels, cells, by = 'Barcode')
key <- key[order(as.numeric(key$cellnum)),]
head(key)

library(Rtsne)

# Assume you have a data matrix 'data'
# Perform t-SNE
tsne_result <- Rtsne(d, verbose = T, partial_pca = F)

# Extract the t-SNE coordinates
tsne_x <- tsne_result$Y[, 1]
tsne_y <- tsne_result$Y[, 2]

tsne <- cbind(Barcode = key$Barcode, cluster =key$GEX.Graph.based, tsne_x, tsne_y)

# Create a scatterplot of the t-SNE embedding
library(ggplot2)
ggplot(tsne, aes(x = tsne_x, y = tsne_y, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_discrete(name = "Cluster") +
  labs(x = "tSNE Dimension 1", y = "tSNE Dimension 2") +
  theme_minimal() +
  theme(
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) -> p
p

png(file.path(output_path, paste0("tsne_pbmc_2cluster.png")), width = plot_width, height = plot_height, units = "px")
print(p)
dev.off()

library(umap)

# Perform UMAP
umap_result <- umap(d)

# Extract the UMAP coordinates
umap_x <- umap_result$layout[, 1]
umap_y <- umap_result$layout[, 2]

# Combine the UMAP coordinates with additional information
umap_data <- cbind(Barcode = key$Barcode, cluster = key$GEX.Graph.based, umap_x, umap_y)

# Create the scatterplot using ggplot2 without axis labels
ggplot(umap_data, aes(x = umap_x, y = umap_y, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_discrete(name = "Cluster") +
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) -> p
p

png(file.path(output_path, paste0("umap_pbmc_2cluster.png")), width = plot_width, height = plot_height, units = "px")
print(p)
dev.off()



gex_umap_path <- file.path(path_data, "GEX UMAP-Projection.csv")
gex_labels_path <- file.path(path_data, "GEX Graph-Based.csv")


gex_umap <- read.csv(gex_umap_path)
gex_labels <- read.csv(gex_labels_path)
gex_data <- merge(gex_umap, gex_labels[, c("Barcode", "GEX.Graph.based")], by = "Barcode", all.x = TRUE)

# Create the scatterplot using ggplot2 without axis labels
ggplot(gex_data, aes(x = X.Coordinate, y = Y.Coordinate, color = GEX.Graph.based)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_discrete(name = "Cluster") +
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) -> p
p

png(file.path(output_path, paste0("umap_pbmc_complete.png")), width = plot_width, height = plot_height, units = "px")
print(p)
dev.off()


