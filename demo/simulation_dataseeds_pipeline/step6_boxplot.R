# !/usr/bin/Rscript

library(RColorBrewer)
library(tidyverse)
library(hrbrthemes)
library(viridis)

# Plot Rand index vs max_loglikelihood ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df_mp <- data.frame(
  i_scenario = integer(),
  rand_index = numeric(),
  n_iterations = integer(),
  converged = logical(),
  silhouette_measure = numeric(),
  davies_bouldin_index = numeric(),
  BIC = numeric(),
  stringsAsFactors = FALSE
)
bics <- c()
i_counter <- 0
for(i_current_scenario_type in seq(20)){
  for(i in seq(10)){

    i_counter <- i_counter+1
    for(i_seed in seq_along(bisc_results_list[[i_counter]]$RIs)){
      current_results <- bisc_results_list[[i_counter]]$bisc_results[[i_seed]]
      bics <- c(bics, NA)
      if(length(current_results) > 1 || !is.na(current_results)){
        new_df <- data.frame(i_scenario = i_current_scenario_type,
                             rand_index = current_results$rand_index,
                             n_iterations = current_results$n_iterations,
                             converged = current_results$converged,
                             silhouette_measure = current_results$silhouette_measure,
                             davies_bouldin_index =  current_results$davies_bouldin_index,
                             BIC = Rmpfr::asNumeric(current_results$max_loglikelihood[[current_results$n_iterations]])
        )
        df_mp <- rbind(df_mp, new_df)
        bics[length(bics)] <- Rmpfr::asNumeric(current_results$max_loglikelihood[[current_results$n_iterations]])
      }
    }
  }
}
# df_mp <- df_mp[df_mp$converged,]
df_mp <- df_mp[,c(1,2,4,7)]
df_mp$i_scenario <- as.numeric(df_mp$i_scenario)
df_mp$i_scenario <- as.factor(df_mp$i_scenario)




# Initialize an empty vector to store the normalized BIC values
df_mp$LL <- numeric(nrow(df_mp))

# Loop through each unique scenario and normalize the BIC values
unique_scenarios <- unique(df_mp$i_scenario)
for (scenario in unique_scenarios) {
  scenario_data <- df_mp[df_mp$i_scenario == scenario, ]
  min_bic <- min(scenario_data$BIC)
  max_bic <- max(scenario_data$BIC)
  range_bic <- max_bic - min_bic

  df_mp$LL[df_mp$i_scenario == scenario] <- (scenario_data$BIC - min_bic) / range_bic
}

df_mp <- tibble::as_tibble(df_mp)

df_mp$BIC <- df_mp$LL
bics[!is.na(bics)] <- df_mp$BIC

p <- ggplot(df_mp, aes(x = BIC, y = rand_index, color = converged)) +
  geom_point(size = 4, alpha = 0.5) +
  facet_wrap(~ i_scenario, scales = "free", labeller = as_labeller(function(x) {
    ifelse(as.numeric(x) <= 10,
           paste("Simple scenario", x),
           paste("Complex scenario", as.numeric(x) - 10))
  })) +
  labs(
    # title = "", # "Rand Index vs Normalised BIC by Scenario",
    x = "Normalised BIC (within each scenario)",
    y = "Rand Index"
  ) +
  guides(color = guide_legend(title = "Converged")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),  # X-axis label size
        axis.title.y = element_text(size = 16),  # Y-axis label size
        axis.text.x = element_text(size = 12),   # X-axis tick labels
        axis.text.y = element_text(size = 12),   # Y-axis tick labels
        legend.text = element_text(size = 12),   # Legend text size
        legend.title = element_text(size = 16)   # Legend title size
  )

# scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"),
#                    name = "Bisc converge")
# print(p)
ggsave(file.path(output_path, paste0("RI_vs_BIC_per_scenario.pdf")), p, width = 12, height = 8)




# Other plots -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


string_counts <- table(unlist(lapply(scenarios, function(x) x$description)))
types <-  rownames(string_counts)
n_types <- length(types)

cells <- list()
cells$Simple <- c()
cells$Complex <- c()

cells_genes <- list()
cells_genes$Simple <- c()
cells_genes$Complex <- c()

biclust <- list()
biclust$Simple <- c()
biclust$Complex <- c()

biclust_genes <- list()
biclust_genes$Simple <- c()
biclust_genes$Complex <- c()

genes <- list()
genes$Simple <- c()
genes$Complex <- c()

cells_var <- list()
cells_var$Simple <- c()
cells_var$Complex <- c()

biclust_var <- list()
biclust_var$Simple <- c()
biclust_var$Complex <- c()

genes_var <- list()
genes_var$Simple <- c()
genes_var$Complex <- c()

cells_biclust <- list()
cells_biclust$Simple <- c()
cells_biclust$Complex <- c()

biclust_biclust <- list()
biclust_biclust$Simple <- c()
biclust_biclust$Complex <- c()

genes_biclust <- list()
genes_biclust$Simple <- c()
genes_biclust$Complex <- c()

cells_biclust_var <- list()
cells_biclust_var$Simple <- c()
cells_biclust_var$Complex <- c()

biclust_biclust_var <- list()
biclust_biclust_var$Simple <- c()
biclust_biclust_var$Complex <- c()

genes_biclust_var <- list()
genes_biclust_var$Simple <- c()
genes_biclust_var$Complex <- c()

scenario_id_bisc_converged <- list()
scenario_id_bisc_converged$Simple <- c()
scenario_id_bisc_converged$Complex <- c()

scenario_id_all <- list()
scenario_id_all$Simple <- c()
scenario_id_all$Complex <- c()

scenario_id_bcplaid_genes <- list()
scenario_id_bcplaid_genes$Simple <- c()
scenario_id_bcplaid_genes$Complex <- c()

scenario_id_bisc_converged_genes <- list()
scenario_id_bisc_converged_genes$Simple <- c()
scenario_id_bisc_converged_genes$Complex <- c()

bics_cells <- c()
bics_cells$Simple <- c()
bics_cells$Complex <- c()

bics_genes <- c()
bics_genes$Simple <- c()
bics_genes$Complex <- c()

genes_module_id <- c()
genes_module_id$Simple <- c()
genes_module_id$Complex <- c()

biclust_genes_module_id <- c()
biclust_genes_module_id$Simple <- c()
biclust_genes_module_id$Complex <- c()

converged <- c()
converged$Simple <- c()
converged$Complex <- c()

converged_genes <- c()
converged_genes$Simple <- c()
converged_genes$Complex <- c()

null2NA <- function(x){
  ifelse(is.null(x), NA, x)
}

i_counter <- 0
for(i_current_scenario_type in seq(20)){
  temp_genes_var <- data.frame()
  temp_cells_var <- c()
  temp_biclust_var <- c()
  temp_biclust_genes_var <- data.frame()
  temp_biclust_cells_var <- c()
  temp_biclust_biclust_var <- c()
  for(i in seq(10)){

    print(paste0('Now running outer iteration ', i_current_scenario_type, i ))
    # print(bisc_results_list[[i]]$RIs[[1]]$RI_cell_clustering_bisc)


    i_counter <- i_counter+1
    for(i_seed in seq_along(bisc_results_list[[i_counter]]$RIs)){

      current_scenario_type <- scenarios[[i_counter]]$description
      if(length(bisc_results_list[[i_counter]]$bisc_results[[i_seed]])>1 || !is.na(bisc_results_list[[i_counter]]$bisc_results[[i_seed]])){

        new_cell_val <- null2NA(bisc_results_list[[i_counter]]$RIs[[i_seed]]$RI_cell_clustering_bisc)
        new_biclust_val <- null2NA(bisc_results_list[[i_counter]]$RIs[[i_seed]]$RI_biclust_bisc)
        new_genes_val <- bisc_results_list[[i_counter]]$RIs[[i_seed]]$RI_gene_clustering_bisc
        cells[[current_scenario_type]] <- c(cells[[current_scenario_type]], new_cell_val)
        biclust[[current_scenario_type]] <- c(biclust[[current_scenario_type]], new_biclust_val)
        genes[[current_scenario_type]] <- c(genes[[current_scenario_type]], new_genes_val)

        cells_genes[[current_scenario_type]] <- c(cells_genes[[current_scenario_type]], rep(new_cell_val, length(new_genes_val)))
        biclust_genes[[current_scenario_type]] <- c(biclust_genes[[current_scenario_type]], rep(new_biclust_val, length(new_genes_val)))


        scenario_id_bisc_converged[[current_scenario_type]] <- c(scenario_id_bisc_converged[[current_scenario_type]], i_current_scenario_type)
        scenario_id_bisc_converged_genes[[current_scenario_type]] <- c(scenario_id_bisc_converged_genes[[current_scenario_type]], rep(i_current_scenario_type, length(new_genes_val)))
        bics_cells[[current_scenario_type]] <- c(bics_cells[[current_scenario_type]], bics[i_counter])
        bics_genes[[current_scenario_type]] <- c(bics_genes[[current_scenario_type]], rep(bics[i_counter], length(new_genes_val)))
        genes_module_id[[current_scenario_type]] <- c(genes_module_id[[current_scenario_type]], seq(length(new_genes_val)))
        converged[[current_scenario_type]] <- c(converged[[current_scenario_type]], bisc_results_list[[i_counter]]$bisc_results[[i_seed]]$converged)
        converged_genes[[current_scenario_type]] <- c(converged_genes[[current_scenario_type]], rep(bisc_results_list[[i_counter]]$bisc_results[[i_seed]]$converged, length(new_genes_val)))

        # Only count variance if converged
        if(bisc_results_list[[i_counter]]$bisc_results[[i_seed]]$converged){
          temp_genes_var <- rbind(temp_genes_var,  as.data.frame(t(new_genes_val)))
          temp_cells_var <- c(temp_cells_var, new_cell_val)
          temp_biclust_var <- c(temp_biclust_var, new_biclust_val)
        }

        new_cell_biclust_val <- null2NA(biclustbiclust_results_list[[i_counter]][[i_seed]]$RI_cell_clustering_biclustbiclust)
        new_biclust_biclust_val <- null2NA(biclustbiclust_results_list[[i_counter]][[i_seed]]$RI_biclust_biclustbiclust)
        new_genes_biclust_val <- as.numeric(strsplit(trimws( null2NA(biclustbiclust_results_list[[i_counter]][[i_seed]]$RI_gene_clustering_biclustbiclust_all)), " ")[[1]])
        cells_biclust[[current_scenario_type]] <- c(cells_biclust[[current_scenario_type]], new_cell_biclust_val)
        biclust_biclust[[current_scenario_type]] <- c(biclust_biclust[[current_scenario_type]], new_biclust_biclust_val)
        genes_biclust[[current_scenario_type]] <- c(genes_biclust[[current_scenario_type]], new_genes_biclust_val)

        temp_biclust_genes_var <- rbind(temp_biclust_genes_var,  as.data.frame(t(new_genes_biclust_val)))
        temp_biclust_cells_var <- c(temp_biclust_cells_var, new_cell_biclust_val)
        temp_biclust_biclust_var <- c(temp_biclust_biclust_var, new_biclust_biclust_val)
        scenario_id_all[[current_scenario_type]] <- c(scenario_id_all[[current_scenario_type]], i_current_scenario_type)
        scenario_id_bcplaid_genes[[current_scenario_type]] <- c(scenario_id_bcplaid_genes[[current_scenario_type]], rep(i_current_scenario_type, length(new_genes_biclust_val)))
        biclust_genes_module_id[[current_scenario_type]] <- c(biclust_genes_module_id[[current_scenario_type]], seq(length(new_genes_biclust_val)))
      }
    }
  }

  if(nrow(temp_genes_var)!=0){
    genes_var[[current_scenario_type]] <- c(genes_var[[current_scenario_type]], sapply(temp_genes_var, var))
  }
  if(length(temp_cells_var)!=0){
    cells_var[[current_scenario_type]] <- c(cells_var[[current_scenario_type]], var(temp_cells_var))
  }
  if(length(temp_biclust_var)!=0){
    biclust_var[[current_scenario_type]] <- c(biclust_var[[current_scenario_type]], var(temp_biclust_var))
  }

  if(nrow(temp_biclust_genes_var)!=0){
    genes_biclust_var[[current_scenario_type]] <- c(genes_biclust_var[[current_scenario_type]], sapply(temp_biclust_genes_var, var))
  }
  if(length(temp_biclust_cells_var)!=0){
    cells_biclust_var[[current_scenario_type]] <- c(cells_biclust_var[[current_scenario_type]], var(temp_biclust_cells_var))
  }
  if(length(temp_biclust_biclust_var)!=0){
    biclust_biclust_var[[current_scenario_type]] <- c(biclust_biclust_var[[current_scenario_type]], var(temp_biclust_biclust_var))
  }
}

# Create the data frame

var_data <- bind_rows(
  tibble(method="bisc", scenario = "Simple", type="Cells", value = unlist(cells_var$Simple)),
  tibble(method="bisc", scenario = "Complex", type="Cells", value = unlist(cells_var$Complex)),
  tibble(method="BCPlaid", scenario = "Simple", type="Cells", value = unlist(cells_biclust_var$Simple)),
  tibble(method="BCPlaid", scenario = "Complex", type="Cells", value = unlist(cells_biclust_var$Complex)),
  tibble(method="bisc", scenario = "Simple", type="Biclust", value = unlist(biclust_var$Simple)),
  tibble(method="bisc", scenario = "Complex", type="Biclust", value = unlist(biclust_var$Complex)),
  tibble(method="BCPlaid", scenario = "Simple", type="Biclust", value = unlist(biclust_biclust_var$Simple)),
  tibble(method="BCPlaid", scenario = "Complex", type="Biclust", value = unlist(biclust_biclust_var$Complex)),
  tibble(method="bisc", scenario = "Simple", type="Genes", value = unlist(genes_var$Simple)),
  tibble(method="bisc", scenario = "Complex", type="Genes", value = unlist(genes_var$Complex)),
  tibble(method="BCPlaid", scenario = "Simple", type="Genes", value = unlist(genes_biclust_var$Simple)),
  tibble(method="BCPlaid", scenario = "Complex", type="Genes", value = unlist(genes_biclust_var$Complex))
)



cell_data <- bind_rows(
  tibble(method="bisc",    type = "Simple",  scenario = unlist(scenario_id_bisc_converged$Simple),  converged=converged$Simple, bic = unlist(bics_cells$Simple), value = unlist(cells$Simple)),
  tibble(method="bisc",    type = "Complex", scenario = unlist(scenario_id_bisc_converged$Complex), converged=converged$Complex, bic = unlist(bics_cells$Complex), value = unlist(cells$Complex)),
  tibble(method="BCPlaid", type = "Simple",  scenario = unlist(scenario_id_all$Simple),             value = unlist(cells_biclust$Simple)),
  tibble(method="BCPlaid", type = "Complex", scenario = unlist(scenario_id_all$Complex),            value = unlist(cells_biclust$Complex))
)

cell_data <- cell_data %>%
  group_by(scenario) %>%
  mutate(lowest_bic = ifelse(method == "bisc" & bic == min(bic[method == "bisc"]), TRUE, FALSE)) %>%
  ungroup()

biclust_data <- bind_rows(
  tibble(method="bisc",    type = "Simple",  scenario = unlist(scenario_id_bisc_converged$Simple), converged=converged$Simple,  bic = unlist(bics_cells$Simple), value = unlist(biclust$Simple)),
  tibble(method="bisc",    type = "Complex", scenario = unlist(scenario_id_bisc_converged$Complex), converged=converged$Complex,  bic = unlist(bics_cells$Complex), value = unlist(biclust$Complex)),
  tibble(method="BCPlaid", type = "Simple",  scenario = unlist(scenario_id_all$Simple),   value = unlist(biclust_biclust$Simple)),
  tibble(method="BCPlaid", type = "Complex", scenario = unlist(scenario_id_all$Complex),  value = unlist(biclust_biclust$Complex))
)

biclust_data <- biclust_data %>%
  group_by(scenario) %>%
  mutate(lowest_bic = ifelse(method == "bisc" & bic == min(bic[method == "bisc"]), TRUE, FALSE)) %>%
  ungroup()

gene_data <- bind_rows(
  tibble(method="bisc", type = "Simple", cells_ri=cells_genes$Simple, biclust_ri =biclust_genes$Simple, scenario = unlist(scenario_id_bisc_converged_genes$Simple), converged=converged_genes$Simple, module_id = genes_module_id$Simple, bic = unlist(bics_genes$Simple), value = unlist(genes$Simple)),
  tibble(method="bisc", type = "Complex", cells_ri=cells_genes$Complex, biclust_ri =biclust_genes$Complex, scenario = unlist(scenario_id_bisc_converged_genes$Complex), converged=converged_genes$Complex, module_id = genes_module_id$Complex, bic = unlist(bics_genes$Complex), value = unlist(genes$Complex)),
  tibble(method="BCPlaid", type = "Simple", scenario = unlist(scenario_id_bcplaid_genes$Simple), module_id = biclust_genes_module_id$Simple, value = unlist(genes_biclust$Simple)),
  tibble(method="BCPlaid", type = "Complex", scenario = unlist(scenario_id_bcplaid_genes$Complex), module_id = biclust_genes_module_id$Complex, value = unlist(genes_biclust$Complex))
)

gene_data <- gene_data %>%
  group_by(scenario) %>%
  mutate(lowest_bic = ifelse(method == "bisc" & bic == min(bic[method == "bisc"]), TRUE, FALSE)) %>%
  ungroup()


plot_height <- 400
plot_width <- 750


#---------------------------------

constructed_plot <- ggplot(var_data, aes(x = interaction(scenario, type), y = value, fill = method)) +
  geom_boxplot() +
  geom_point(aes(color = method), position=position_jitterdodge(), alpha=0.3, show.legend = FALSE) +
  theme_minimal(base_size = 16) +
  labs(fill = "Method") +
  # ggtitle("Rand index variance over different run seeds") +
  xlab("ScenarioType.Type") +
  ylab("Variance") +
  ylim(0, 0.064) +
  theme(axis.title.x = element_text(size = 16),  # X-axis label size
        axis.title.y = element_text(size = 16),  # Y-axis label size
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),   # X-axis tick labels
        axis.text.y = element_text(size = 12),   # Y-axis tick labels
        legend.text = element_text(size = 12),   # Legend text size
        legend.title = element_text(size = 16)   # Legend title size
  )
# print(constructed_plot)
ggsave(file.path(output_path, paste0("variance_dataseeds.pdf")), constructed_plot, width = 8, height = 8)


# Plot RI_vs_BIC_per_scenario_for_bisc_biclusters.pdf
p <- ggplot(biclust_data[biclust_data$method=="bisc",], aes(x = bic, y = value, color = converged)) +
  geom_point(size = 4, alpha = 0.5) +
  facet_wrap(~ scenario, scales = "free", labeller = as_labeller(function(x) {
    ifelse(as.numeric(x) <= 10,
           paste("Simple scenario", x),
           paste("Complex scenario", as.numeric(x) - 10))
  })) +
  labs(
    # title = "", # "Rand Index vs Normalised BIC by Scenario",
    x = "Normalised BIC (within each scenario)",
    y = "Rand Index"
  ) +
  guides(color = guide_legend(title = "Converged")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),  # X-axis label size
        axis.title.y = element_text(size = 16),  # Y-axis label size
        axis.text.x = element_text(size = 12),   # X-axis tick labels
        axis.text.y = element_text(size = 12),   # Y-axis tick labels
        legend.text = element_text(size = 12),   # Legend text size
        legend.title = element_text(size = 16)   # Legend title size
  )


# print(p)
ggsave(file.path(output_path, paste0("RI_vs_BIC_per_scenario_for_bisc_biclusters.pdf")), p, width = 12, height = 8)

gene_data$module_id <- as.factor(gene_data$module_id)

# Plot RI_vs_BIC_per_scenario_for_bisc_gene_modules.pdf
p <- ggplot(gene_data[gene_data$method=="bisc",], aes(x = bic, y = value, color=converged, shape  = module_id)) +
  geom_point(size = 4, alpha = 0.5) +
  facet_wrap(~ scenario, scales = "free", labeller = as_labeller(function(x) {
    ifelse(as.numeric(x) <= 10,
           paste("Simple scenario", x),
           paste("Complex scenario", as.numeric(x) - 10))
  })) +
  labs(
    x = "Normalised BIC (within each scenario)",
    y = "Rand Index"
  ) +
  guides(shape = guide_legend(title = "Module ID"),
         color = guide_legend(title = "Converged")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)
  )

# print(p)
ggsave(file.path(output_path, paste0("RI_vs_BIC_per_scenario_for_bisc_gene_modules.pdf")), p, width = 12, height = 8)

# -------------------------------

gene_data$gene_module_ri <- gene_data$value
bisc_gene_data <- gene_data[gene_data$method=="bisc",]
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

pairs(bisc_gene_data[,c(3,4,11)],
      lower.panel = panel.cor,
      upper.panel = upper.panel)

#--------------
p <- ggplot(gene_data[gene_data$method=="bisc",], aes(x = cells_ri, y = gene_module_ri, color=converged, shape  = module_id)) +
  geom_point(size = 4, alpha = 0.3) +
  facet_wrap(~ scenario, scales = "free", labeller = as_labeller(function(x) {
    ifelse(as.numeric(x) <= 10,
           paste("Simple scenario", x),
           paste("Complex scenario", as.numeric(x) - 10))
  })) +
  labs(
    x = "Cell clustering rand index",
    y = "Gene module rand index"
  ) +
  guides(shape = guide_legend(title = "Module ID"),
         color = guide_legend(title = "Converged")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)
  )

# print(p)
ggsave(file.path(output_path, paste0("CellRI_vs_GeneRI_per_scenario_for_bisc.pdf")), p, width = 12, height = 8)

#--------------
p <- ggplot(gene_data[gene_data$method=="bisc",], aes(x = biclust_ri, y = gene_module_ri, color=converged, shape  = module_id)) +
  geom_point(size = 4, alpha = 0.5) +
  facet_wrap(~ scenario, scales = "free", labeller = as_labeller(function(x) {
    ifelse(as.numeric(x) <= 10,
           paste("Simple scenario", x),
           paste("Complex scenario", as.numeric(x) - 10))
  })) +
  labs(
    x = "Biclust clustering rand index",
    y = "Gene module rand index"
  ) +
  guides(shape = guide_legend(title = "Module ID"),
         color = guide_legend(title = "Converged")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)
  )

# print(p)
ggsave(file.path(output_path, paste0("BiclustRI_vs_GeneRI_per_scenario_for_bisc.pdf")), p, width = 12, height = 8)


#--------------
# Plot RI_vs_BIC_per_scenario_for_bisc_gene_modules.pdf
p <- ggplot(gene_data[gene_data$method=="bisc",], aes(x = cells_ri, y = biclust_ri, color=converged)) +
  geom_point(size = 4, alpha = 0.5) +
  facet_wrap(~ scenario, scales = "free", labeller = as_labeller(function(x) {
    ifelse(as.numeric(x) <= 10,
           paste("Simple scenario", x),
           paste("Complex scenario", as.numeric(x) - 10))
  })) +
  labs(
    x = "Cell clustering rand index",
    y = "Biclustering rand index"
  ) +
  guides(color = guide_legend(title = "Converged")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)
  )

# print(p)
ggsave(file.path(output_path, paste0("CellRI_vs_BiclustRI_per_scenario_for_bisc.pdf")), p, width = 12, height = 8)



# Cell clustering RI results per scenario
# cell_data[is.na(cell_data$converged),'converged'] <- TRUE
# cell_data <- cell_data[cell_data$converged,]
p <- ggplot(cell_data, aes(x = method, y = value, color = method)) +
  geom_boxplot() +
  facet_wrap(~ scenario, scales = "free", labeller = as_labeller(function(x) {
    ifelse(as.numeric(x) <= 10,
           paste("Simple scenario", x),
           paste("Complex scenario", as.numeric(x) - 10))
  })) +
  geom_point(
    data = subset(cell_data, method == "bisc" & lowest_bic),
    aes(shape = "Best LL"),
    color = "red",
    size = 3,
    alpha = 0.4
  ) +
  labs(
    x = "Method",
    y = "Rand Index"
  ) +
  guides(color = guide_legend(title = "Method"),
         shape = guide_legend(title = "")) +
  scale_shape_manual(values = c("Best LL" = 17)) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),  # X-axis label size
        axis.title.y = element_text(size = 16),  # Y-axis label size
        axis.text.x = element_blank(),           # Remove X-axis tick labels
        axis.ticks.x = element_blank(),          # Remove X-axis ticks
        axis.text.y = element_text(size = 12),   # Y-axis tick labels
        legend.text = element_text(size = 12),   # Legend text size
        legend.title = element_text(size = 16)   # Legend title size
  )

print(p)
ggsave(file.path(output_path, paste0("boxplot_cell_cluster_RI_comparison_dataseeds.pdf")), p, width = 10, height = 7)


# Biclustering RI results per scenario
# biclust_data[is.na(biclust_data$converged),'converged'] <- TRUE
# biclust_data <- biclust_data[biclust_data$converged,]
p <- ggplot(biclust_data, aes(x = method, y = value, color = method)) +
  geom_boxplot() +
  facet_wrap(~ scenario, scales = "free", labeller = as_labeller(function(x) {
    ifelse(as.numeric(x) <= 10,
           paste("Simple scenario", x),
           paste("Complex scenario", as.numeric(x) - 10))
  })) +
  geom_point(
    data = subset(biclust_data, method == "bisc" & lowest_bic == TRUE),
    aes(shape = "Best LL"),
    color = "red",
    size = 3,
    alpha = 0.4
  ) +
  labs(
    # title = "", # "Rand Index vs Normalised BIC by Scenario",
    x = "Method",
    y = "Rand Index"
  ) +
  guides(color = guide_legend(title = "Method"),
         shape = guide_legend(title = "")) +
  scale_shape_manual(values = c("Best LL" = 17)) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),  # X-axis label size
        axis.title.y = element_text(size = 16),  # Y-axis label size
        axis.text.x = element_blank(),           # Remove X-axis tick labels
        axis.ticks.x = element_blank(),          # Remove X-axis ticks
        axis.text.y = element_text(size = 12),   # Y-axis tick labels
        legend.text = element_text(size = 12),   # Legend text size
        legend.title = element_text(size = 16)   # Legend title size
  )

# print(p)
ggsave(file.path(output_path, paste0("boxplot_biclust_RI_comparison_dataseeds.pdf")), p, width = 10, height = 7)

# Gene module clustering RI results per scenario
# gene_data[is.na(gene_data$converged),'converged'] <- TRUE
# gene_data <- gene_data[gene_data$converged,]
p <- ggplot(gene_data, aes(x = method, y = value, color = method)) +
  geom_boxplot() +
  geom_point(
    data = subset(gene_data, method == "bisc" & lowest_bic == TRUE),
    aes(shape = "Best LL"),
    color = "red",
    size = 3,
    alpha = 0.4
  ) +
  facet_wrap(~ scenario, scales = "free", labeller = as_labeller(function(x) {
    ifelse(as.numeric(x) <= 10,
           paste("Simple scenario", x),
           paste("Complex scenario", as.numeric(x) - 10))
  })) +
  labs(
    x = "Method",
    y = "Rand Index"
  ) +
  guides(
    color = guide_legend(title = "Method"),
    shape = guide_legend(title = "")
  ) +
  scale_shape_manual(values = c("Best LL" = 17)) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16)
  )
# print(p)
ggsave(file.path(output_path, paste0("boxplot_gene_modules_RI_comparison_dataseeds.pdf")), p, width = 10, height = 7)


#
# constructed_plot <- ggplot(cell_data, aes(x = interaction(method, type), y = value, fill = type)) +
#   geom_violin(position = "dodge") +
#   stat_summary(fun = median, geom = "crossbar", width = 0.6, color = "cyan") +
#   scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
#   geom_jitter(color = "black",  size=0.8, width = 0.3, alpha = 0.5, height = 0) +
#   theme_ipsum() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 11)
#   ) +
#   ggtitle("Cell clusters rand index comparison") +
#   xlab("Method.Scenario") +
#   ylab("RI") +
#   ylim(0, 1)
# print(constructed_plot)
# png(file.path(output_path, paste0("boxplot_cell_cluster_RI_comparison_dataseeds.png")), width = plot_width, height = plot_height, units = "px")
# print(constructed_plot)
# dev.off()
#
# constructed_plot <- ggplot(biclust_data, aes(x = interaction(method, type), y = value, fill = type)) +
#   geom_violin(position = "dodge") +
#   stat_summary(fun = median, geom = "crossbar", width = 0.6, color = "cyan") +
#   scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
#   geom_jitter(color = "black",  size=0.8, width = 0.3, alpha = 0.5, height = 0) +
#   theme_ipsum() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 11)
#   ) +
#   ggtitle("Biclust rand index comparison") +
#   xlab("Method.Scenario") +
#   ylab("RI") +
#   ylim(0, 1)
# print(constructed_plot)
# png(file.path(output_path, paste0("boxplot_biclust_RI_comparison_dataseeds.png")), width = plot_width, height = plot_height, units = "px")
# print(constructed_plot)
# dev.off()
#
# constructed_plot <- ggplot(gene_data, aes(x = interaction(method, type), y = value, fill = type)) +
#   geom_violin(position = "dodge") +
#   stat_summary(fun = median, geom = "crossbar", width = 0.6, color = "cyan") +
#   scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
#   geom_jitter(color = "black", size=0.8, width = 0.3, alpha = 0.5, height = 0) +
#   theme_ipsum() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 11)
#   ) +
#   ggtitle("Gene modules rand index comparison") +
#   xlab("Method.Scenario") +
#   ylab("RI") +
#   ylim(0, 1)
# print(constructed_plot)
# png(file.path(output_path, paste0("boxplot_gene_modules_RI_comparison_dataseeds.png")), width = plot_width, height = plot_height, units = "px")
# print(constructed_plot)
# dev.off()
