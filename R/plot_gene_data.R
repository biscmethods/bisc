#!/usr/bin/Rscript

library(ggplot2)  # To plot things #TODO: What things?

# Plot original data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_onex_to_oney_original_data <- function(dat=dat,
                                            plot_path=demo_path,
                                            filename="scatterplot_data.png") {


  if (length(ind_targetgenes) == 1 && length(ind_reggenes) == 1) {
    p <- ggplot(dat,
                aes(x = r1,
                    y = t1,
                    color = true_cell_cluster_allocation)
    )
    png(file.path(plot_path, filename))
    p + geom_point()
    dev.off()
  }
}