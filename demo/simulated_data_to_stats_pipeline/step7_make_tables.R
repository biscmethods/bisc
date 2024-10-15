stats2table <- function(tibble, outfile = '/out_temp.tex'){
  summary_stats <- summary(gene_data$value)
  summary_df <- data.frame(Statistic = names(summary_stats), Value = as.numeric(summary_stats))


  library(xtable)
  latex_table <- xtable(summary_df)
  latex_code <- print(latex_table, type = "latex", include.rownames = FALSE)


  sink(paste0(output_path,'/',outfile))
  cat(latex_code)
  sink()
}

stats2table(cell_data, 'cell_data.tex')
stats2table(gene_data, 'gene_data.tex')
stats2table(biclust_data, 'biclust_data.tex')
