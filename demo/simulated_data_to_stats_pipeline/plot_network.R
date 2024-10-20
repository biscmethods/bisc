# This plots which regulators are affecting which target gene modules in the generated data
# Basically we want to translate the guessed labels to the true labels in our guesses to start with
# (in a cell cluster we guess which target genes belong to which gene module)
# After that we move the rows in e.g. true_S (or the corresponding matrix in our guess) so that
# the rows are in the same order in the guess and the genereated data.
# Then we can plot comapritive networks graphs?
# OR we can calculate the accuracy of our guesses by comparing the true labels to the guessed labels
# E.g. by ROC.





# Network plot of generated data, which regulator genes affect each target gene module ----------------------------------------------------------------------------------------------------------------------------------------------------


# Load required libraries
library(igraph)
library(ggplot2)
library(ggraph)

# Create a function to process the data and create an edge list
create_edge_list <- function(data) {
  edge_list <- data.frame()
  for (cluster in 1:length(data)) {
    for (module in 1:nrow(data[[cluster]])) {
      for (gene in 1:ncol(data[[cluster]])) {
        if (data[[cluster]][module, gene] != 0) {
          edge_list <- rbind(edge_list, data.frame(
            from = paste0("Gene_", gene),
            to = paste0("Module_", module, "_Cluster_", cluster)
          ))
        }
      }
    }
  }
  return(edge_list)
}


edge_list <- create_edge_list(scenarios[[11]]$generated_data$true_S)

# Create the graph
g <- graph_from_data_frame(edge_list, directed = FALSE)

# Set node types (gene or module)
V(g)$type <- ifelse(grepl("Gene", V(g)$name), "Gene", "Module")

# Create the plot
p <- ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(aes(color = type, size = type)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("Gene" = "lightblue", "Module" = "lightgreen")) +
  scale_size_manual(values = c("Gene" = 5, "Module" = 7)) +
  theme_void() +
  theme(legend.position = "bottom")

# Display the plot
print(p)

# Save the plot as a PNG file
ggsave("network_graph.png", p, width = 12, height = 10, dpi = 300)
