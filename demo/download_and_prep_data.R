if (!require(Matrix)) install.packages('Matrix')
if (!require(Seurat)) install.packages('Seurat')
if (!require(here)) install.packages('here')
if (!require(utils)) install.packages('utils')
if (!require(httr)) install.packages('httr')
if (!require(jsonlite)) install.packages('jsonlite')
if (!require(BiocManager)) install.packages('BiocManager')  # For glmGamPoi that is used inside Seurat to speed things up
BiocManager::install('glmGamPoi')

library(Matrix)
library(Seurat)  # To work with Neftel data

# This script will automatically download the Neftel2019 dataset and
# add it to biclust/data/Neftel2019.
# You can download Group1 and/or Group2 from https://uppsala.box.com/s/aqf3fyeb43i5jy8v4s527u0l262virue if it fails.
# This then process the datasets into rds and/or loads the datasets.

# Which dataset to load
group_1_flag <- TRUE
group_2_flag <- TRUE

# Folders and filenames
path_temp <- here::here('temp')
path_data <- here::here('data')
path_save_zip <- file.path(path_temp, "neftel2019.zip")
path_Neftel2019 <- file.path(path_temp, "Neftel2019")

# Create paths that don't exist
if (!file.exists(path_data)) {
  dir.create(path_data)
}
if (!file.exists(path_temp)) {
  dir.create(path_temp)
}

# Download data
if (!file.exists(path_save_zip)) {
  get_download_link_url <- 'https://uppsala.app.box.com/index.php?folder_id=174080811166&q[shared_item][shared_name]=aqf3fyeb43i5jy8v4s527u0l262virue&rm=box_v2_zip_shared_folder'
  response <- httr::GET(url = get_download_link_url)
  download_url <- jsonlite::fromJSON(rawToChar(response$content))$download_url
  # File is 453.1 MB so it takes a while
  options(timeout = max(3000, getOption("timeout")))
  utils::download.file(url = download_url, destfile = path_save_zip, quiet = TRUE, cacheOK = FALSE, mode = 'wb')
  if (!file.exists(path_Neftel2019)) {
    utils::unzip(zipfile = path_save_zip, exdir = path_temp)
  }
}

# Group 1 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if (group_1_flag) {
  path_Group1 <- file.path(path_Neftel2019, "Group1")
  path_mtx_file <- file.path(path_Group1, 'Exp_data_UMIcounts_10X.mtx')
  path_cell_file <- file.path(path_Group1, 'Cells_10X.txt')
  path_gene_file <- file.path(path_Group1, 'Genes_10X.txt')
  path_neftel_mn_group1 <- file.path(path_data, "neftel_mn_group1.rds")
  path_neftel_seurat_group1 <- file.path(path_data, "neftel_seurat_group1.rds")

  if (!file.exists(path_Group1)) {
    stop(paste("Please download Group1 from https://uppsala.box.com/s/aqf3fyeb43i5jy8v4s527u0l262virue and put Group1 folder in", path_Neftel2019))
  }

  if (!file.exists(path_neftel_mn_group1)) {
    mn_g1 <- Matrix::readMM(file = path_mtx_file)
    mn_g1 <- as.matrix(mn_g1)

    cells <- read.table(file = path_cell_file, sep = ' ', header = TRUE, stringsAsFactors = FALSE)
    genes <- read.table(file = path_gene_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
    genes <- genes[, 1]

    rownames(cells) <- cells[, 1]

    rownames(mn_g1) <- genes
    colnames(mn_g1) <- cells$cell_name

    mn_g1 <- mn_g1[!duplicated(genes),]
    genes <- genes[!duplicated(genes)]

    saveRDS(mn_g1, file = file.path(path_neftel_mn_group1))

  }else {
    mn_g1 <- readRDS(file = path_neftel_mn_group1)

    cells <- read.table(file = path_cell_file, sep = ' ', header = TRUE, stringsAsFactors = FALSE)
    genes <- read.table(file = path_gene_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
    genes <- genes[, 1]
    genes <- genes[!duplicated(genes)]

    rownames(cells) <- cells[, 1]
  }

  if (!file.exists(path_neftel_seurat_group1)) {
    Neftel_g1 <- Seurat::CreateSeuratObject(mn_g1, min.cells = 3, min.features = 500, meta.data = cells)
    Neftel_g1 <- Seurat::SCTransform(Neftel_g1)

    Neftel_g1 <- Seurat::RunPCA(Neftel_g1, verbose = FALSE)
    Neftel_g1 <- Seurat::RunUMAP(Neftel_g1, dims = 1:30, verbose = FALSE)

    saveRDS(Neftel_g1, file = path_neftel_seurat_group1)
  }else {
    Neftel_g1 <- readRDS(file = path_neftel_seurat_group1)
  }
}


# Group 2 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if (group_2_flag) {
  path_neftel_mn_group2 <- file.path(path_data, "neftel_mn_group2.rds")
  path_neftel_seurat_group2 <- file.path(path_data, "neftel_seurat_group2.rds")
  path_Group2 <- file.path(path_Neftel2019, "Group2")

  if (!file.exists(path_Group2)) {
    dir.create(path_Group2)
    stop(paste("Please download Group2 from https://uppsala.box.com/s/aqf3fyeb43i5jy8v4s527u0l262virue and put them in", path_Group2))
  }

  if (!file.exists(path_neftel_mn_group2)) {
    mn_g2 <- Matrix::readMM(file = file.path(path_Group2, 'exp_data_TPM.mtx'))
    mn_g2 <- as.matrix(mn_g2)

    cells <- read.table(file = file.path(path_Group2, 'Cells.txt'), sep = ' ', header = TRUE, stringsAsFactors = FALSE)
    genes <- read.table(file = file.path(path_Group2, 'Genes.txt'), sep = '\t', header = FALSE, stringsAsFactors = FALSE)
    genes <- genes[, 1]

    rownames(cells) <- cells[, 1]

    rownames(mn_g2) <- genes
    colnames(mn_g2) <- cells$cell_name

    saveRDS(mn_g2, file = path_neftel_mn_group2)
    base::rm(cells, genes, group_2_flag)

  }else {
    mn_g2 <- readRDS(file = path_neftel_mn_group2)
  }

  if (!file.exists(path_neftel_seurat_group2)) {
    cells <- read.table(file = file.path(path_Group2, 'Cells.txt'), sep = ' ', header = TRUE, stringsAsFactors = FALSE)
    rownames(cells) <- cells[, 1]

    Neftel_g2 <- Seurat::CreateSeuratObject(mn_g2, min.cells = 3, min.features = 500, meta.data = cells)
    Neftel_g2 <- Seurat::SCTransform(Neftel_g2)

    Neftel_g2 <- Seurat::RunPCA(Neftel_g2, verbose = FALSE)
    Neftel_g2 <- Seurat::RunUMAP(Neftel_g2, dims = 1:30, verbose = FALSE)

    saveRDS(Neftel_g2, file = path_neftel_seurat_group2)
    base::rm(cells)
  }else {
    Neftel_g2 <- readRDS(file = path_neftel_seurat_group2)
  }
}


