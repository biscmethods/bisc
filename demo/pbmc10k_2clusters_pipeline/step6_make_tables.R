library(dplyr)
library(xtable)

outfile <- 'bisc_pbmc_cellclust_RI.tex'

caption = 'Rand Indexes on cell cluster allocation by convergence status'

# Assuming you have the 'converged_RIs' and 'nonconverged_RIs' vectors
RI_data <- data.frame(
  RI = c(converged_RIs, nonconverged_RIs),
  Converged = c(rep("Converged", length(converged_RIs)),
                rep("Nonconverged", length(nonconverged_RIs)))
)


# Calculate the summary statistics for each group
RI_summary <- RI_data %>%
  group_by(Converged) %>%
  summarize(
    Min = min(RI, na.rm = TRUE),
    Q1 = quantile(RI, 0.25, na.rm = TRUE),
    Median = median(RI, na.rm = TRUE),
    Mean = mean(RI, na.rm = TRUE),
    Q3 = quantile(RI, 0.75, na.rm = TRUE),
    Max = max(RI, na.rm = TRUE)
  )

# Create the LaTeX table
# latex_code <- kable(RI_summary, format = "latex", digits = 5, booktabs = TRUE) %>%
#   kable_styling(latex_options = c("striped", "hold_position")) %>%
#   column_spec(1, bold = TRUE)


latex_table <- xtable(RI_summary,  digits = 4)
latex_code <- print(latex_table,
                    type = "latex",
                    include.rownames = FALSE,
                    booktabs = TRUE)

latex_code_ <- gsub("\\\\begin\\{table\\}\\[ht\\]",
                    paste0("\\\\begin{table}[ht]\n\\\\caption{", # wtf
                           caption,
                           "}"#, "\n\\\\vspace{1.5cm}"
                           ),
                            latex_code
                    )

sink(paste0(output_path, "/", outfile))
cat(latex_code_)
sink()
