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
#' @param coefficient_mean Mean coefficients in true model, length n_target_gene_clusters
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
    n_target_genes = 50,  # Number of target genes
    n_regulator_genes = 30,  # Number of regulator genes
    n_cells  = 1000,  # Number of cells
    n_target_gene_clusters  = 3,  # Number of target gene clusters
    regulator_mean   = 1,  # Mean expression of regulator genes
    coefficient_mean = c(1,20,100)  # Mean coefficients in true model, length n_target_gene_clusters
)
{
  # Check arguments
  checks <- function(x, one_element=TRUE, atomic=TRUE, numeric=TRUE, positive=TRUE, int=TRUE){
    error_message <- ""
    if(one_element){
      if (!(length(x) == 1L)) {
        add_message <- paste0(deparse(substitute(x)), " must be one number.")
        error_message <- paste(error_message, add_message, sep="\n")
      }
    }

    if(atomic){
      if(!is.atomic(x)) {
        add_message <- (paste0(deparse(substitute(x)), " must be atomic."))
        error_message <- paste(error_message, add_message, sep="\n")
      }
    }

    if(numeric){
      if (!is.numeric(x)) {
        add_message <- (paste0(deparse(substitute(x)), " must be numeric."))
        error_message <- paste(error_message, add_message, sep="\n")
      }
    }

    if(positive){
      if (x < 0) {
        add_message <- (paste0(deparse(substitute(x)), " must be >=0."))
        error_message <- paste(error_message, add_message, sep="\n")
      }
    }

    if(int){
      if (!(round(x)==x)) {
        add_message <- (paste0(deparse(substitute(x)), " must be an integer."))
        error_message <- paste(error_message, add_message, sep="\n")
      }
    }

    if(length(error_message)>1){
      stop(error_message)
    }
  }


  checks(n_target_genes)
  checks(n_regulator_genes)
  checks(n_cells)
  checks(n_target_gene_clusters)
  checks(regulator_mean, one_element=TRUE, atomic=TRUE, numeric=TRUE, positive=FALSE, int=FALSE)
  checks(coefficient_mean, one_element=FALSE, atomic=TRUE, numeric=TRUE, positive=FALSE, int=FALSE)

  if (length(coefficient_mean) != n_target_gene_clusters) {
    stop(paste0("coefficient_mean has length ", length(coefficient_mean), ", but must have length n_target_gene_clusters: ",n_target_gene_clusters,"."))
  }

  # It's +1 because kmeans algo doesn't work otherwise
  if (!(n_target_genes >= (n_target_gene_clusters+1))) {
    stop(paste0("It must be that: n_target_genes>=n_target_gene_clusters+1, but right now: ", n_target_genes, "<", n_target_gene_clusters))
  }


  diag_ <- function(vec)diag(x = vec, nrow = length(vec))
  # set.seed(3333) #To get a nice matrix
  Pi <- matrix(0, n_target_gene_clusters, n_target_genes)

  # Minst en etta per rad. Exakt en etta per kolumn.
  # To guarantee at least one target gene per cluster,
  # assign the first n_target_gene_clusters target
  # gene to different clusters
  for(i_target_gene_cluster in 1:n_target_gene_clusters){
    for(i_target_gene in i_target_gene_cluster:min(i_target_gene_cluster, ncol(Pi))){
      # print(paste("Target gene cluster", i_target_gene_cluster, "Target gene", i_target_gene))
      Pi[i_target_gene_cluster, i_target_gene] <- 1
    }
  }

  # Set rest of matrix randomly
  if(ncol(Pi)>=(n_target_gene_clusters + 1)){
    for(index in (n_target_gene_clusters + 1):n_target_genes){
      random_cluster_index <- sample.int(n = n_target_gene_clusters, size = 1)
      Pi[random_cluster_index, index] <- 1
    }
  }

  if(!all(colSums(Pi)==1)){
    stop("True cluster allocation matrix is wrong. At least one target gene belongs to several clusters.")
  }

  if(!all(rowSums(Pi)>=1)){
    stop("True cluster allocation matrix is wrong. At least one target gene cluster doesn't have a target gene.")
  }

  # Binary matrix R ---------------------------------------------------------
  # The regulators (columns) that affect each cluster (rows, i in the manuscript).
  # Note that in the paper this is a vector with set of indexes,
  # not a binary matrix, see code for an example.
  # Note that regulators can affect any number of clusters.

  # set.seed(12345)  # To get a nice matrix
  r_data <- rbinom(n_target_gene_clusters * n_regulator_genes, 1, 1/n_target_gene_clusters)
  R <- matrix(data=r_data,
              nrow=n_target_gene_clusters,
              ncol=n_regulator_genes)  # Row i is an indicator version of R_i in manuscript

  # Check if any cluster (row) has no regulator, and if so assign one
  for(row in 1:n_target_gene_clusters){
    if(sum(R[row,]) == 0){
      R[row, sample.int(n_regulator_genes, 1)] <- 1
    }
  }

  # R

  # R[1,]  # Cluster 1 is affected by these regulators
  # sum(R[1,])  # R_1 in the manuscript is which regulators affect cluster 1

  R2R_i <- function(i) {
    which(R[i,]!= 0)
  }
  # R2R_i(1) # R_1 in the manuscript is which regulators affect cluster 1

  # Matrix S ----------------------------------------------------------------
  # A n_target_gene_clusters x n_regulator_genes matrix with 1 or -1 if the regulator (columns)
  # of how (1 stimulating, -1 repressing) a regulator affects a cluster (rows),
  # 0 if it doesn't affect it.
  # This has the same information as the manuscript's s_i
  # set.seed(10) # to get a nice matrix
  s_data <- rbinom(n = n_target_gene_clusters * n_regulator_genes, 1, 0.8)*2-1
  S <- R * matrix(data=s_data,
                  nrow=n_target_gene_clusters,
                  ncol=n_regulator_genes)  # Just randomize signs

  S2S_i <- function(i, mat=S) {
    mat[i, which(mat[i,]!= 0) ]
  }
  # S2S_i(2)  # Non-zero entries of this is s_i in manuscript

  # For each regulator j in module i, a mean value was
  # chosen uniformly at random between 0.01 and 0.1
  # Regulator_means<- runif(n_regulator_genes, 0.01, 0.1)
  # Regulator_means

  # Matrix Zr ---------------------------------------------------------------
  # n_cells x n_regulator_genes, cells are rows, regulator genes are columns
  # Just get some random expression for regulator genes for now
  Z_r <- matrix(data = rnorm(n_cells * n_regulator_genes, mean = regulator_mean, sd = 0.1),
                nrow = n_cells,
                ncol = n_regulator_genes)
  # Z_r <- sapply(1:n_regulator_genes, function(i)  rnorm(n_cells, mean = runif(1, 0.1,1), sd = 0.1) )

  # Z_r
  # dim(Z_r)

  # Array 𝚩 ---------------------------------------------------------------
  # Now we want to build 𝚩 and use 𝚩 to build Z_t.
  # For that we need some coefficients from our regression models.
  # These coefficients are stored in the array Beta, which has
  # dimension  n_regulator_genes x n_target_genes x n_target_gene_clusters with nonegative entries
  # in this we store the (in the manuscript only the non-zero) coefficients
  # describing how the regulator genes affect the target genes.

  # Beta <- array(data = abs(rnorm(n_regulator_genes * n_target_genes * n_target_gene_clusters,
  #                                 mean = coefficient_mean, sd = 0.1)),
  #               c(n_regulator_genes,n_target_genes,n_target_gene_clusters))

  Beta <- array(
    data = sapply(1:n_target_gene_clusters, function(i) rnorm(n_regulator_genes*n_target_genes, mean = coefficient_mean[i], sd = 0.1)),
    dim = c(n_regulator_genes,n_target_genes,n_target_gene_clusters)
  )

  # Beta <- array(data = unlist(
  #   sapply(1:n_target_gene_clusters, function(i) rnorm(n_regulator_genes*n_target_genes, mean = runif(1,min = 1, max = 1), sd = 0.1), simplify = F)
  # ), dim = c(n_regulator_genes,n_target_genes,n_target_gene_clusters) )

  # Make 𝚩 zero in appropriate spots
  for (clust in 1:n_target_gene_clusters){
    Beta[,,clust] <-  diag_(R[clust,]) %*% Beta[,,clust]
  }
  # Beta
  # dim(Beta)

  # In the manuscript the zero rows are just dropped

  Beta2Beta_i <- function(i){
    matrix(data = Beta[R2R_i(i),,i], nrow = length(R2R_i(i)), ncol = n_target_genes)
  }

  # Beta2Beta_i(1) #Beta_i as in the manuscript, has dimension |R_i| x n_target_genes
  # Beta[,,1]

  # Matrix Z_t --------------------------------------------------------------
  # If j is one target cell, that cells expression should then be,
  # according to (1) in the manuscript.

  # Z_t could below be initialized to something nonzero, and in the next step the
  # right hand side could be added instead of just inserted, this would make the
  # initialisation similar to some baseline exposure, or intercept in the model.

  Z_t <- matrix(data = 0, nrow = n_cells, ncol=n_target_genes)

  # Produce a Beta-vector with signs and true cluster allocation (which means we zero out columns with Pi)
  # This is just to output the correct Betas for debugging
  Beta_with_signs <- vector("list", length = n_target_gene_clusters)
  for(i_target_gene_cluster in 1:n_target_gene_clusters){
    Beta_with_signs[[i_target_gene_cluster]] <-  (diag_(S[i_target_gene_cluster,]) %*% Beta[,,i_target_gene_cluster])  %*% diag_(Pi[i_target_gene_cluster,])
  }

  # Todo: Vectorize this
  # Create Z_t
  for(i in 1:n_target_gene_clusters){
    for(j in 1:n_target_genes){
      Z_t[,j] <-
        Z_t[,j] +
        Pi[i,j] *                #  True cluster allocation, zero if Z_t[,j] is not in cluster i
        (
          Z_r[,R2R_i(i)] %*%     # Gene expression of regulators of cluster i
            diag_(S2S_i(i)) %*%  # Signs for whether regulators are stimulating or repressing
            Beta2Beta_i(i)[,j]   # How much reg of cluster i affects target j
        )


    }

    # cat(paste0("building cluster ", i,"\n_cells"))
  }

  # Check if generated data gives rand index 1. If not stop execution
  scregclust::scregclust(
    expression = rbind(t(Z_t), t(Z_r)),    #scRegClust wants this form
    genesymbols = 1:(n_target_genes+n_regulator_genes),               #gene row numbers
    is_regulator = (1:(n_target_genes+n_regulator_genes) > n_target_genes) + 0, #vector indicating which genes are regulators
    n_cl        = n_target_gene_clusters,
    penalization = 0.001,
    verbose = FALSE
  )-> scRegOut

  true_clust_allocation <- apply(X=Pi, MARGIN=2, FUN=function(x) which(x==1))
  predicted_cluster_allocation <- scRegOut$results[[1]]$output[[1]]$cluster[1:n_target_genes]
  rand_index <- aricode::RI(true_clust_allocation, predicted_cluster_allocation)

  if(rand_index!=1){
    stop("scregclust couldn't find correct clusters in generated data. Rand index:", rand_index)
  }

  # This can probably be vectorized
  # For this we are omitting the variance terms.
  # Z_t, Z_r, S_i, and B_i as here will minimize (1) in the manuscript
  list(Z_t = Z_t,
       Z_r = Z_r,
       Pi  = Pi,
       R   = R,
       S   = S,
       B = Beta_with_signs
  )
}

# runs only when script is run by itself
# || interactive()
if (sys.nframe() == 0){
  # ... do main stuff
  res <- generate_dummy_data_for_scregclust(50, 30, 1000, 3, 1, c(1,20,100))
  print(str(res))
}
