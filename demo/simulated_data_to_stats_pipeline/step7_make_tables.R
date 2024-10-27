  library(dplyr)
  library(xtable)

  stats2table <- function(tibble, outfile = '/out_temp.tex', caption = 'Your caption here'){
    # Group by method and calculate summary statistics
    summary_stats <- tibble %>%
      group_by(method, type) %>%
      summarise(
        Min = min(value, na.rm = TRUE),
        Q1 = quantile(value, 0.25, na.rm = TRUE),
        Median = median(value, na.rm = TRUE),
        Mean = mean(value, na.rm = TRUE),
        Q3 = quantile(value, 0.75, na.rm = TRUE),
        Max = max(value, na.rm = TRUE)
      )

    # Convert to data frame for xtable
    summary_df <- as.data.frame(summary_stats)

    # Replace "BB" with "biclust::biclust()" and "BS" with "bisc"
    summary_df$method <- gsub("BB", "biclust::biclust()", summary_df$method)
    summary_df$method <- gsub("BS", "bisc", summary_df$method)

    # Create LaTeX table
    latex_table <- xtable(summary_df)
    latex_code <- print(latex_table,
                        type = "latex",
                        include.rownames = FALSE,
                        booktabs = TRUE)

    # Insert the caption after \begin{table}[ht]
    latex_code_ <- gsub("\\\\begin\\{table\\}\\[ht\\]",
                        paste0("\\\\begin{table}[ht]\n\\\\caption{", caption, "}\n\\\\vspace{1.5cm}"),
                        latex_code)

    # Fix headers
    latex_code_ <- gsub("method", "Method", latex_code_)
    latex_code_ <- gsub("type", '\\\\makecell{Regulator \\\\\\\\ Structure}', latex_code_)


    # Write to filey
    sink(paste0(output_path, "/", outfile))
    cat(latex_code_)
    sink()
  }

  # Example usage
  stats2table(tibble = cell_data,    outfile = 'cell_data.tex',    caption = 'Rand Indexes on cell cluster allocation')
  stats2table(tibble = gene_data,    outfile = 'gene_data.tex',    caption = 'Rand Indexes on gene module allocation')
  stats2table(tibble = biclust_data, outfile = 'biclust_data.tex', caption = 'Rand Indexes on bicluster allocation')
