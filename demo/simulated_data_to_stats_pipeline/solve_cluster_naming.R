# Install and load required packages
if (!require(clue)) install.packages("clue")
library(clue)


# Find the correct mapping of what guesses from biclustbiclust and biclustscreg should match true generated data. For gene modules --------------------------------------------------------------------------------------------------------
# Correct = "The most advantageous labeling that produces the least amount of error"
get_mapping <- function(true_labels,
                        guessed_labels,
                        verbose=FALSE){
  # Function to find optimal cluster label mapping
  find_optimal_mapping <- function(true_labels, guessed_labels) {
    # Remove entries where guessed_labels is -1
    valid_indices <- which(guessed_labels != -1)
    true_labels_valid <- true_labels[valid_indices]
    guessed_labels_valid <- guessed_labels[valid_indices]

    # Get unique cluster labels
    true_clusters <- sort(unique(true_labels))
    guessed_clusters <- sort(unique(guessed_labels_valid))

    # Create confusion matrix
    conf_matrix <- table(true_labels_valid, guessed_labels_valid)

    # Ensure the confusion matrix is square by adding rows/columns of zeros if necessary
    max_clusters <- max(length(true_clusters), length(guessed_clusters))
    if (nrow(conf_matrix) < max_clusters) {
      conf_matrix <- rbind(conf_matrix, matrix(0, max_clusters - nrow(conf_matrix), ncol(conf_matrix)))
    }
    if (ncol(conf_matrix) < max_clusters) {
      conf_matrix <- cbind(conf_matrix, matrix(0, nrow(conf_matrix), max_clusters - ncol(conf_matrix)))
    }

    # Solve the assignment problem
    assignment <- solve_LSAP(conf_matrix, maximum = TRUE)

    # Create the mapping
    mapping <- data.frame(
      true_cluster = true_clusters,
      guessed_cluster = guessed_clusters[assignment],
      stringsAsFactors = FALSE
    )

    # Add -1 to the mapping
    mapping <- rbind(mapping, data.frame(true_cluster = NA, guessed_cluster = -1))

    return(mapping)
  }

  # Function to apply the mapping and calculate accuracy
  apply_mapping_and_calculate_accuracy <- function(true_labels, guessed_labels, mapping) {
    # Create a named vector for easy mapping
    translation <- setNames(mapping$true_cluster, mapping$guessed_cluster)

    # Translate guessed labels
    translated_guessed <- translation[as.character(guessed_labels)]

    # Calculate accuracy, ignoring -1 in guessed_labels
    valid_indices <- which(guessed_labels != -1)
    accuracy <- mean(true_labels[valid_indices] == translated_guessed[valid_indices], na.rm = TRUE)

    # Calculate the proportion of unassigned labels
    unassigned_prop <- mean(guessed_labels == -1)

    return(list(translated_guessed = translated_guessed,
                accuracy = accuracy,
                unassigned_proportion = unassigned_prop))
  }

  # Example usage
  # set.seed(123)  # for reproducibility

  optimal_mapping <- find_optimal_mapping(true_labels, guessed_labels)

  result <- apply_mapping_and_calculate_accuracy(true_labels, guessed_labels, optimal_mapping)

  # Replace NA in guess_cluster column with value in true cluster then remove rows with NA
  optimal_mapping[is.na(optimal_mapping$guessed_cluster),]$guessed_cluster <- optimal_mapping[is.na(optimal_mapping$guessed_cluster),]$true_cluster
  optimal_mapping_no_na <- na.omit(optimal_mapping)

  if(verbose){
    print(guessed_labels)
    print(optimal_mapping)
    print(optimal_mapping_no_na)
    print(unname(result$translated_guessed))
    print(true_labels)
    print(paste("Accuracy:", result$accuracy))
    print(paste("Proportion of unassigned labels:", result$unassigned_proportion))
  }
  return(optimal_mapping_no_na)
}


# Returns a df with ratio of correct guesses of regulator genes affecting a gene module ---------------------------------------------------------------------------------------------------------------------------------------------------
for(scenario_number in 1:length(scenarios)){
  guessed_S <- lapply(biclust_screg_results_list[[scenario_number]]$biclustscreg_results[[1]]$scregclust_final_result_models, t)
  guessed_S <- lapply(guessed_S, function(x) x * 1)
  true_S <- lapply(scenarios[[scenario_number]]$generated_data$true_S, abs)
  for(i_cluster in 1:scenarios[[scenario_number]]$n_cell_clusters){
    mapping <- get_mapping(true_labels = scenarios[[scenario_number]]$generated_data$true_target_gene_allocation[[i_cluster]],
                           guessed_labels = biclust_screg_results_list[[scenario_number]]$biclustscreg_results[[1]]$scregclust_final_result_module[[i_cluster]][scenarios[[scenario_number]]$ind_targetgenes]
    )
    guessed_S[[i_cluster]] <- guessed_S[[i_cluster]][mapping$guessed_cluster,]

  }
  # print(guessed_S)
  # print(true_S)


  # Function to count correct and incorrect guesses for each pair of matrices
  compare_matrices <- function(true_matrix, guess_matrix) {
    total_elements <- length(true_matrix)  # Total number of elements

    TP_ratio <- sum(true_matrix == 1 & guess_matrix == 1) / total_elements # True Positives
    FP_ratio <- sum(true_matrix == 0 & guess_matrix == 1) / total_elements # False Positives
    TN_ratio <- sum(true_matrix == 0 & guess_matrix == 0) / total_elements # True Negatives
    FN_ratio <- sum(true_matrix == 1 & guess_matrix == 0) / total_elements # False Negatives
    correct_ratio <- sum(true_matrix == guess_matrix)/total_elements

    return(list(correct_ratio=correct_ratio, TP_ratio = TP_ratio, FP_ratio = FP_ratio, TN_ratio = TN_ratio, FN_ratio = FN_ratio))
  }

  # Apply the function to each pair of matrices using mapply
  ratios_list <- mapply(compare_matrices, true_S, guessed_S, SIMPLIFY = FALSE)
  # Convert the list of ratios into a data frame
  if(scenario_number==1){
    ratios_df <- do.call(rbind, lapply(ratios_list, as.data.frame))
    ratios_df <- cbind(gene_module_name = 1:scenarios[[scenario_number]]$n_cell_clusters, ratios_df)
    ratios_df <- cbind(scenario = scenario_number, ratios_df)
  } else{
    ratios_df_tmp <- do.call(rbind, lapply(ratios_list, as.data.frame))
    ratios_df_tmp <- cbind(gene_module_name = 1:scenarios[[scenario_number]]$n_cell_clusters, ratios_df_tmp)
    ratios_df_tmp <- cbind(scenario = scenario_number, ratios_df_tmp)

    ratios_df <- rbind(ratios_df, ratios_df_tmp)
  }
}
# Print the results
print(ratios_df)


