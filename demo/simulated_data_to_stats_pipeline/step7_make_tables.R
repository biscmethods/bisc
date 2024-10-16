stats2table <- function(tibble, outfile = '/out_temp.tex', caption = 'Your caption here'){
  summary_stats <- summary(tibble$value)
  summary_df <- data.frame(Statistic = names(summary_stats), Value = as.numeric(summary_stats))

  library(xtable)
  latex_table <- xtable(summary_df)
  latex_code <- print(latex_table,
                      type = "latex",
                      include.rownames = FALSE,
                      booktabs = TRUE)

  # Insert the caption after \begin{table}[ht]
  latex_code_ <- gsub("\\\\begin\\{table\\}\\[ht\\]", paste0("\\\\begin{table}[ht]\n\\\\caption{", caption, "}"), latex_code)

  sink(paste0(output_path, "/", outfile))
  cat(latex_code_)
  sink()
}

# Example usage
stats2table(tibble = cell_data,    outfile = 'cell_data.tex',    caption = 'Rand Indexes on cell cluster allocation')
stats2table(tibble = gene_data,    outfile = 'gene_data.tex',    caption = 'Rand Indexes on gene module allocation')
stats2table(tibble = biclust_data, outfile = 'biclust_data.tex', caption = 'Rand Indexes on bicluster allocation')
