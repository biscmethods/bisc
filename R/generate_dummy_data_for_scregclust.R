#!/usr/bin/Rscript

#' Dummy data generation for scRegClust
#'
#' Generates dummy data that work with scRegClust.
#' Rows are cluster index. Cols are target gene index.
#' @param n_target_genes The number of target genes
#' @param n_regulator_genes The number of regulator genes
#' @param n_cells The number of cells
#' @param n_target_gene_clusters The number of target gene clusters
#' @param regulator_mean Mean expression of regulator genes
#' @param coefficient_sd Mean coefficients in true model, length n_target_gene_clusters
#' @return A list with
#'   \item{Z_t}{a}
#'   \item{Z_r}{a}
#'   \item{Pi}{n_target_gene_clusters x n_target_genes. True cell cluster allocation.}
#'   \item{R}{a}
#'   \item{S}{a}
#'   \item{B}{a}
#' @examples
#' res <- generate_dummy_data_for_scregclust(50, 30, 1000, 3, 1, c(1,20,100));
#' @export
generate_dummy_data_for_scregclust <- function(
  n_target_genes = 2193,  # Number of target genes               2193
  n_regulator_genes = 493,  # Number of regulator genes          493
  n_cells = 6000,  # Number of cells                             6000
  n_target_gene_clusters = 10,  # Number of target gene clusters 10
  regulator_mean = 0,  # Mean expression of regulator genes      0
  coefficient_mean = c(0.0417,
                       0.0343,
                       0.0576,
                       0.043 ,
                       0.0576,
                       0.0413,
                       0.0473,
                       0.0444,
                       0.0481,
                       -0.0139),
  coefficient_sd = c(0.0556,
                     0.037,
                     0.0638,
                     0.0466,
                     0.0761,
                     0.0471,
                     0.0468,
                     0.0611,
                     0.0623,
                     0.0394
  ),  # Mean coefficients/betas in true model, length n_target_gene_clusters
  make_regulator_network_plot = TRUE,
  plot_suffix = "vignette",
  testing_penalization = 0.1 # optional, for the screg run in the end
)
{
  # Check arguments
  checks <- function(x, one_element = TRUE, atomic = TRUE, numeric = TRUE, positive = TRUE, int = TRUE) {
    error_message <- ""
    if (one_element) {
      if (!(length(x) == 1L)) {
        add_message <- paste0(deparse(substitute(x)), " must be one number.")
        error_message <- paste(error_message, add_message, sep = "\n")
      }
    }

    if (atomic) {
      if (!is.atomic(x)) {
        add_message <- (paste0(deparse(substitute(x)), " must be atomic."))
        error_message <- paste(error_message, add_message, sep = "\n")
      }
    }

    if (numeric) {
      if (!is.numeric(x)) {
        add_message <- (paste0(deparse(substitute(x)), " must be numeric."))
        error_message <- paste(error_message, add_message, sep = "\n")
      }
    }

    if (positive) {
      if (!all(x > 0)) {
        add_message <- (paste0(deparse(substitute(x)), " must be >=0."))
        error_message <- paste(error_message, add_message, sep = "\n")
      }
    }

    if (int) {
      if (!(round(x) == x)) {
        add_message <- (paste0(deparse(substitute(x)), " must be an integer."))
        error_message <- paste(error_message, add_message, sep = "\n")
      }
    }

    if (length(error_message) > 1) {
      stop(error_message)
    }
  }


  checks(n_target_genes)
  checks(n_regulator_genes)
  checks(n_cells)
  checks(n_target_gene_clusters)
  checks(regulator_mean, one_element = TRUE, atomic = TRUE, numeric = TRUE, positive = TRUE, int = FALSE)
  checks(coefficient_sd, one_element = FALSE, atomic = TRUE, numeric = TRUE, positive = TRUE, int = FALSE)

  if (length(coefficient_sd) != n_target_gene_clusters) {
    stop(paste0("coefficient_sd has length ", length(coefficient_sd), ", but must have length n_target_gene_clusters: ", n_target_gene_clusters, "."))
  }

  # It's +1 because kmeans algo doesn't work otherwise
  if (!(n_target_genes >= (n_target_gene_clusters + 1))) {
    stop(paste0("It must be that: n_target_genes>=n_target_gene_clusters+1, but right now: ", n_target_genes, "<", n_target_gene_clusters))
  }



  # help functions
  diag_ <- function(vec)diag(x = vec, nrow = length(vec))
  # set.seed(3333) #To get a nice matrix
  Pi <- matrix(0, n_target_gene_clusters, n_target_genes)

  # Minst en etta per rad. Exakt en etta per kolumn.
  # To guarantee at least one target gene per cluster,
  # assign the first n_target_gene_clusters target
  # gene to different clusters
  for (i_target_gene_cluster in 1:n_target_gene_clusters) {
    for (i_target_gene in i_target_gene_cluster:min(i_target_gene_cluster,
                                                    ncol(Pi))) {
      # print(paste("Target gene cluster", i_target_gene_cluster, "Target gene", i_target_gene))
      Pi[i_target_gene_cluster, i_target_gene] <- 1
    }
  }

  # Set rest of matrix randomly
  if (ncol(Pi) >= (n_target_gene_clusters + 1)) {
    for (index in (n_target_gene_clusters + 1):n_target_genes) {
      random_cluster_index <- sample.int(n = n_target_gene_clusters, size = 1)
      Pi[random_cluster_index, index] <- 1
    }
  }

  if (!all(colSums(Pi) == 1)) {
    stop("True cluster allocation matrix is wrong. At least one target gene belongs to several clusters.")
  }

  if (!all(rowSums(Pi) >= 1)) {
    stop("True cluster allocation matrix is wrong. At least one target gene cluster doesn't have a target gene.")
  }

  # Binary matrix R ---------------------------------------------------------
  # The regulators (columns) that affect each cluster (rows, i in the manuscript).
  # Note that in the paper this is a vector with set of indexes,
  # not a binary matrix, see code for an example.
  # Note that regulators can affect any number of clusters.

  cat("Generating R\n")

  r_data <- rbinom(n_target_gene_clusters * n_regulator_genes, 1, 1 / n_target_gene_clusters)
  R <- matrix(data = r_data,
              nrow = n_target_gene_clusters,
              ncol = n_regulator_genes)  # Row i is an indicator version of R_i in manuscript

  # Check if any cluster (row) has no regulator, and if so assign one
  for (row in 1:n_target_gene_clusters) {
    if (sum(R[row,]) == 0) {
      R[row, sample.int(n_regulator_genes, 1)] <- 1
    }
  }

  # R

  # R[1,]  # Cluster 1 is affected by these regulators
  # sum(R[1,])  # R_1 in the manuscript is which regulators affect cluster 1

  R2R_i <- function(i) {
    which(R[i,] != 0)
  }

  # R2R_i(1) # R_1 in the manuscript is which regulators affect cluster 1

  # Matrix S ----------------------------------------------------------------
  # A n_target_gene_clusters x n_regulator_genes matrix with 1 or -1 if the regulator (columns)
  # of how (1 stimulating, -1 repressing) a regulator affects a cluster (rows),
  # 0 if it doesn't affect it.
  # This has the same information as the manuscript's s_i
  # set.seed(10) # to get a nice matrix

  cat("Generating S\n")

  s_data <- rbinom(n = n_target_gene_clusters * n_regulator_genes, 1, 0.8) * 2 - 1
  S <- R * matrix(data = s_data,
                  nrow = n_target_gene_clusters,
                  ncol = n_regulator_genes)  # Just randomize signs

  S2S_i <- function(i, mat = S) {
    mat[i, which(mat[i,] != 0)]
  }

  # S2S_i(2)  # Non-zero entries of this is s_i in manuscript

  # For each regulator j in module i, a mean value was
  # chosen uniformly at random between 0.01 and 0.1
  # Regulator_means<- runif(n_regulator_genes, 0.01, 0.1)
  # Regulator_means

  # Matrix Zr ---------------------------------------------------------------
  # n_cells x n_regulator_genes, cells are rows, regulator genes are columns
  # Just get some random expression for regulator genes for now
  #Z_r <- matrix(data = rnorm(n_cells * n_regulator_genes, mean = regulator_mean, sd = 0.1),
  #              nrow = n_cells,
  #              ncol = n_regulator_genes)

  cat("Generating Z_r\n")

  Z_r <- matrix(data = rnorm(n_cells * n_regulator_genes,
                             mean = regulator_mean,
                             sd   = 1),
                nrow = n_cells,
                ncol = n_regulator_genes)

  # Z_r <- sapply(1:n_regulator_genes, function(i)  rnorm(n_cells, mean = runif(1, 0.1,1), sd = 0.1) )

  # It should be OK that individual genes are below zero,
  # compare to data in vignette which is centered
  # if (!all(Z_r > 0)) {
  #   stop("Not all regulator genes generated where positive, try with a larger regulator gene mean.")
  # }

  # Z_r
  # dim(Z_r)

  # Array 𝚩 ---------------------------------------------------------------
  # Now we want to build 𝚩 and use 𝚩 to build Z_t.
  # For that we need some coefficients from our regression models.
  # These coefficients are stored in the array Beta, which has
  # dimension  n_regulator_genes x n_target_genes x n_target_gene_clusters with non-negative entries
  # in this we store the (in the manuscript only the non-zero) coefficients
  # describing how the regulator genes affect the target genes.

  # Beta <- array(data = abs(rnorm(n_regulator_genes * n_target_genes * n_target_gene_clusters,
  #                                 mean = coefficient_sd, sd = 0.1)),
  #               c(n_regulator_genes,n_target_genes,n_target_gene_clusters))
  # Coefficients should have mean zero

  cat("Generating Coefficients\n")

  Beta <- array(
    data = sapply(1:n_target_gene_clusters,
                  function(i) abs(rnorm(n_regulator_genes * n_target_genes,
                                        mean = coefficient_mean[i],
                                        sd = coefficient_sd[i]))),
    dim = c(n_regulator_genes, n_target_genes, n_target_gene_clusters)
  )


  if (!all(Beta > 0)) {
    stop("Not all betas generated where positive, try with a larger coefficient mean.")
  }

  # Beta <- array(data = unlist(
  #   sapply(1:n_target_gene_clusters, function(i) rnorm(n_regulator_genes*n_target_genes, mean = runif(1,min = 1, max = 1), sd = 0.1), simplify = F)
  # ), dim = c(n_regulator_genes,n_target_genes,n_target_gene_clusters) )

  # Make 𝚩 zero in appropriate spots
  for (clust in 1:n_target_gene_clusters) {
    Beta[, , clust] <- diag_(R[clust,]) %*% Beta[, , clust]
  }
  # Beta
  # dim(Beta)

  # In the manuscript the zero rows are just dropped

  Beta2Beta_i <- function(i) {
    matrix(data = Beta[R2R_i(i), , i], nrow = length(R2R_i(i)), ncol = n_target_genes)
  }

  # Beta2Beta_i(1) #Beta_i as in the manuscript, has dimension |R_i| x n_target_genes
  # Beta[,,1]

  # Matrix Z_t --------------------------------------------------------------
  # If j is one target cell, that cells expression should then be,
  # according to (1) in the manuscript.

  # Z_t could below be initialized to something nonzero, and in the next step the
  # right hand side could be added instead of just inserted, this would make the
  # initialisation similar to some baseline exposure, or intercept in the model.

  cat("Generating Target gene expressions\n")

  Z_t <- matrix(data = 0, nrow = n_cells, ncol = n_target_genes)

  # Produce a Beta-vector with signs and true cluster allocation (which means we zero out columns with Pi)
  # This is just to output the correct Betas for debugging
  Beta_with_signs <- vector("list", length = n_target_gene_clusters)
  for (i_target_gene_cluster in 1:n_target_gene_clusters) {
    Beta_with_signs[[i_target_gene_cluster]] <- (diag_(S[i_target_gene_cluster,]) %*% Beta[, , i_target_gene_cluster]) %*% diag_(Pi[i_target_gene_cluster,])
  }

  # i_target_gene_cluster = 1
  # i_target_gen = 1

  # Create Z_t
  for (i_target_gene_cluster in 1:n_target_gene_clusters) {
    for (i_target_gen in 1:n_target_genes) {

      diagonal_matrix <- diag_(S2S_i(i_target_gene_cluster))

        target_gene <-
          (
            # Gene expression of regulators of cluster i
            Z_r[, R2R_i(i_target_gene_cluster)]
            %*%
              # Signs for whether regulators are stimulating or repressing
              diagonal_matrix
            %*%
              # How much reg of cluster i affects target j
              Beta2Beta_i(i_target_gene_cluster)[, i_target_gen]
          )

      error_mean <- 0  # mean(target_gene) / 10

      # logic check

      tryCatch(
        {
          error_term <- rnorm(n_cells,
                              mean = error_mean,
                              sd = abs(mean(target_gene)) / 100)
        },
        error = function(cond) {
          NULL
        },
        warning = function(cond) {
          print("Warning generated building error term")
        },
        finally = NULL
      )

      Z_t[, i_target_gen] <-
        Z_t[, i_target_gen] +
          Pi[i_target_gene_cluster, i_target_gen] *                #  True cluster allocation, zero if Z_t[,j] is not in cluster i
            (
              target_gene + error_term
            )


    }

    # cat(paste0("building cluster ", i,"\n_cells"))
  }

  # Check if generated data gives rand index 1. If not stop execution

  cat("Checking if scRegClust can re-create objective clusters\n")


  # set.seed(8374)
  # fit <- scregclust(
  #   z,
  #   genesymbols,
  #   is_regulator,
  #   penalization = seq(0.1, 0.5, 0.1),
  #   n_cl = 10L,
  #   n_cycles = 50L, noise_threshold = 0.05, center = FALSE,
  #   sample_assignment = sample_assignment
  # )

  scregclust::scregclust(
    expression = rbind(t(Z_t), t(Z_r)),  # scRegClust wants this form: row x col: genes x cells.
    genesymbols = 1:(n_target_genes + n_regulator_genes),  # gene row numbers
    is_regulator = (1:(n_target_genes + n_regulator_genes) > n_target_genes) + 0,  # vector indicating which genes are regulators
    n_cl = n_target_gene_clusters,
    penalization = testing_penalization, #generated data is supposed to resemble results from this
    n_cycles     = 10,
    verbose = TRUE,
    min_cluster_size = 1,
    noise_threshold = 0,
    center = T
  ) -> scRegOut

  if(make_regulator_network_plot){

    cat("Plotting regulator network\n")

    demo_path <- here::here("demo")
    R_path <- here::here("R")
    output_path <- demo_path

    tryCatch(
      {
        p <- scregclust::plot_regulator_network(scRegOut$results[[1]]$output[[1]])
      },
      error = function(cond) {
        print("error generated building network plot")
      },
      warning = function(cond) {
        print("Warning generated building network plot")
      },
      finally = function(cond){
        png(file.path(output_path, paste0("screg_network_",plot_suffix ,".png")),
            width = 1024, height = 480, units = "px")
        print(p)
      dev.off()
      }
    )
  }

  # str(scRegOut)
  # str(scRegOut$results)
  # str(scRegOut$results[[1]])
  # str(scRegOut$results[[1]]$output)
  # scRegOut$results[[1]]$output[[1]]$coeffs
  # scRegOut$results[[1]]$output[[1]]$models
  # scRegOut$results[[1]]$output[[1]]$cluster[!is_regulator]
  #
  # beta2scregcoeffmat <- function(mat){
  #    mat[, colSums(mat != 0) > 0]
  # }
  #
  # ncol(beta2scregcoeffmat(Beta_with_signs[[1]]))
  # beta2scregcoeffmat(Beta_with_signs[[2]])
  # ncol(scRegOut$results[[1]]$output[[1]]$coeffs[[1]])
  # scRegOut$results[[1]]$output[[1]]$coeffs[[1]] - beta2scregcoeffmat(Beta_with_signs[[2]])


  #TODO: somehow compare Beta_with_signs and scRegOut$results[[1]]$output[[1]]$coeffs

  true_clust_allocation <- apply(X = Pi, MARGIN = 2, FUN = function(x) which(x == 1))
  # cat(true_clust_allocation)
  predicted_cluster_allocation <- scRegOut$results[[1]]$output[[1]]$cluster_all[1:n_target_genes]
  # cat(predicted_cluster_allocation)
  rand_index <- aricode::RI(true_clust_allocation, predicted_cluster_allocation)
  rand_index <- rand_index - sum(predicted_cluster_allocation == -1) / length(predicted_cluster_allocation)

  cat(paste0("Rand index: ", rand_index))

  if (rand_index < 0.85) {
    cat("scregclust couldn't find correct clusters in generated data. Rand index:", rand_index)
  }

  # This can probably be vectorized
  # For this we are omitting the variance terms.
  # Z_t, Z_r, S_i, and B_i as here will minimize (1) in the manuscript
  list(
    Z_t = Z_t,  # cells x genes
    Z_r = Z_r,  # cells x genes
    Pi = Pi,
    R = R,
    S = S,
    B = Beta_with_signs
  )
}

# runs only when script is run by itself

if (sys.nframe() == 0) {
  # ... do main stuff

  # n_target_genes = 50  # Number of target genes
  # n_regulator_genes = 30  # Number of regulator genes
  # n_cells = 10000 # Number of cells
  # n_target_gene_clusters = 3  # Number of target gene clusters
  # regulator_mean = 1  # Mean expression of regulator genes
  # coefficient_mean = c(5, 20, 100)  # Mean coefficients/betas in true model, length n_target_gene_clusters)
  # coefficient_sd = c(0.1, 0.1, 0.1)
  #
  set.seed(1234)
  res <- generate_dummy_data_for_scregclust()

  print(str(res))
}
