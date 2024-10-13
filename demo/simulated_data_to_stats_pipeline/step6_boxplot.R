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
genes$Simple <- c(length=unname(string_counts["Simple"]))
genes$Complex <- c(length=unname(string_counts["Complex"]))



i_simple <- 0
i_complex <- 0
for(i in 1:length(biclust_screg_results_list)){
  print(paste0('Now running outer iteration ', i))
  print(biclust_screg_results_list[[i]]$RIs[[1]]$RI_cell_clustering_biclustscreg)

  # cells[i] <- biclust_screg_results_list[[i]]$RIs$RI_cell_clustering_biclustscreg

  if(scenarios[[i]]$description=="Simple"){
    i_simple <- i_simple +1
    cells[[scenarios[[i]]$description]][i_simple] <- biclust_screg_results_list[[i]]$RIs[[1]]$RI_cell_clustering_biclustscreg
    biclust[[scenarios[[i]]$description]][i_simple] <- biclust_screg_results_list[[i]]$RIs[[1]]$RI_biclust_biclustscreg
    genes[[scenarios[[i]]$description]][i_simple] <- median(biclust_screg_results_list[[i]]$RIs[[1]]$RI_gene_clustering_biclustscreg)
  }else
  {
    i_complex <- i_complex +1
    cells[[scenarios[[i]]$description]][i_complex] <- biclust_screg_results_list[[i]]$RIs[[1]]$RI_cell_clustering_biclustscreg
    biclust[[scenarios[[i]]$description]][i_complex] <- biclust_screg_results_list[[i]]$RIs[[1]]$RI_biclust_biclustscreg
    genes[[scenarios[[i]]$description]][i_complex] <- median(biclust_screg_results_list[[i]]$RIs[[1]]$RI_gene_clustering_biclustscreg)
  }

}



library(tidyverse)
library(hrbrthemes)
library(viridis)

# Create the data frame
cell_data <- bind_rows(
  tibble(type = "Simple", value = unlist(cells$Simple)),
  tibble(type = "Complex", value = unlist(cells$Complex))
)

biclust_data <- bind_rows(
  tibble(type = "Simple", value = unlist(biclust$Simple)),
  tibble(type = "Complex", value = unlist(biclust$Complex))
)

gene_data <- bind_rows(
  tibble(type = "Simple", value = unlist(genes$Simple)),
  tibble(type = "Complex", value = unlist(genes$Complex))
)

# Plot
cell_data %>%
  ggplot( aes(x=type, y=value, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Cell clusters RI result for biclust_screg") +
  xlab("")

biclust_data %>%
  ggplot( aes(x=type, y=value, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Biclust RI result for biclust_screg") +
  xlab("")


gene_data %>%
  ggplot( aes(x=type, y=value, fill=type)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Gene modules RI result for biclust_screg (median)") +
  xlab("")


