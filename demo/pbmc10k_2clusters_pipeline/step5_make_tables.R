library(dplyr)
library(xtable)

outfile <- 'bisc_pbmc_cellclust_RI.tex'

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
    Min = min(RI),
    `1st Qu.` = quantile(RI, 0.25),
    Median = median(RI),
    Mean = mean(RI),
    `3rd Qu.` = quantile(RI, 0.75),
    Max = max(RI)
  )

# Create the LaTeX table
latex_code_ <- kable(RI_summary, format = "latex", digits = 2, booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped", "hold_position")) %>%
  column_spec(1, bold = TRUE)


latex_table <- xtable(RI_summary)
latex_code <- print(latex_table,
                    type = "latex",
                    include.rownames = FALSE,
                    booktabs = TRUE)

sink(paste0(output_path, "/", outfile))
cat(latex_code)
sink()
