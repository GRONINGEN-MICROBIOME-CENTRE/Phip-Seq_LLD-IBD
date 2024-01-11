#!/usr/bin/env Rscript

# # setting working directory (cluster or local)
# path_cluster <- '/gpfs/projects/bsc83/Projects/scRNAseq/aripol1/CMV-sc/'
# path_local <- '/home/aripol1/Desktop/bsc/Projects/scRNAseq/aripol1/CMV-sc/'
# 
# if(file.exists(path_cluster)){
#   setwd(paste(path_cluster))
# }else if(file.exists(path_local)){
#   setwd(paste(path_local))
# }

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--dataset"), action="store", default=NA, type='character',
              help="v2, v3 or pilot3 filename"),
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell types in each of the cell_level argument."),
  make_option(c("--cell_level"), action="store", default="cell_type", type='character',
              help="Azimuth l1 (predicted.celltype.l1) or l2 (cell_type)."),
  make_option(c("--genes"), action="store", default=NULL, type='character',
              help="Create a new subset of cells."),
  make_option(c("--in_dir"), action="store", default='data/CMV_sc_2_DEA/so_split_by_celltype', type='character',
              help="Input main directory"),
  make_option(c("--out_dir"), action="store", default='output/CMV_sc_2.1_DEA_Pseudobulk', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(data.table))
shhh(library(Seurat))
shhh(library(lme4))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(RColorBrewer))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(ggrepel))
shhh(library(edgeR))
shhh(library(limma))
shhh(library(textTinyR))
shhh(library(pbapply))

#################### Define functions #################### 
# Create pseudobulk profiles
## main function
pseudobulk.func <- function(so, aggregates = 'assignment', vars_model = c('Age', 'Gender', 'date', 'lane'), so_assay = 'RNA'){
  # Create the aggregate_countMatrix (gene ~ aggregates_id)
  ## get the metadata
  metadata <- so@meta.data
  
  ## create aggregated metadata
  # aggregate_metadata <- unique(groups) #only one variable to aggregate
  # aggregate_metadata <- unique(metadata[, aggregates]) #more than one variable to aggregate
  # aggregate_metadata <- unique(metadata[, vars_model]) #we can be removing some individual (aggregates) that have the same vars_model values
  aggregate_metadata <- droplevels(unique(metadata[, c(aggregates,vars_model)]))
  rownames(aggregate_metadata) <- unique(aggregate_metadata[[aggregates]])
  
  ## get donor ids
  IDs <- metadata[,aggregates]
  unique_ID_list <- as.list(unique(IDs))
  
  ## grab the countmatrix
  countMatrix <- GetAssayData(so, assay = so_assay, slot = "counts")
  
  # Pseudobulk by aggregates
  print('Aggregating count matrix using lapply + textTinyR::sparse_Sums() ...')
  system.time(aggregate_countMatrix <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){sparse_Sums(countMatrix[,IDs == x, drop = FALSE], rowSums = TRUE)})))
  # cellcount <- pblapply(unique_ID_list, FUN = function(x){ncol(countMatrix[,IDs == x, drop = FALSE])})
  colnames(aggregate_countMatrix) <- unlist(unique_ID_list)
  rownames(aggregate_countMatrix) <- rownames(countMatrix)
  
  # Filter genes by minimum expression (number of counts) -as bulk-
  print('Normalizing like bulk (limma/DESeq2) using pseudobulk-sum...')
  isexpr = rowSums(cpm(aggregate_countMatrix)>0.1) >= 5
  isexpr_deseq2 = rowSums(aggregate_countMatrix>0) > 10
  aggregate_countMatrix_expr <- aggregate_countMatrix[isexpr,]
  aggregate_countMatrix_expr.deseq2 <- aggregate_countMatrix[isexpr_deseq2,]
  
  # Info
  n_genes <- nrow(countMatrix)
  n_genes_expr <- nrow(aggregate_countMatrix_expr)
  n_genes_expr.deseq2 <- nrow(aggregate_countMatrix_expr.deseq2)
  expressed_genes.prop_round <- round((n_genes_expr/n_genes),2)
  expressed_genes.deseq2.prop_round <- round((n_genes_expr.deseq2/n_genes),2)
  n_cells <- ncol(countMatrix)
  print(paste0('Calculating pseudo-bulk expression matrix for ', opt$cell_type, ':'))
  print(paste0('   ### n_cells: ', n_cells))
  print(paste0('   ### Limma voom --> n_genes (expressed): ', n_genes_expr, ' (proportion of ', expressed_genes.prop_round,')'))
  print(paste0('   ### DESeq2 --> n_genes (expressed): ', n_genes_expr.deseq2, ' (proportion of ', expressed_genes.deseq2.prop_round,')'))
  if(n_genes_expr<2){
    err_in <- paste0('Pseudobulk have not been performed for cell type: ', opt$cell_type, ' in the dataset: ', chem,'. voom() needs at least 2 genes to fit a mean-variance trend. We have ', n_genes_expr, ' genes expressed.')
    stop(err_in)
  }
  
  # Normalization (standard usage of limma/voom)
  geneExpr = DGEList( aggregate_countMatrix_expr )
  geneExpr = calcNormFactors( geneExpr )
  
  # Gene expression matrices (for the correlation analysis) --> use options 2,3 and 4; but save all of them
  ## 1. raw pseudo-bulk counts (sum counts from cells from the same cell type)
  counts <- geneExpr$counts
  
  ## 2. log2 pseudo-bulk counts (sum counts from cells from the same cell type)
  counts_log <- log2(counts+0.5)
  # counts_log <- log2(counts+1)
  
  ## 3. log2CPM normalized using normalization factors -voom- (previously calculated with calcNormFactors() from edgeR)
  # from source code of voom() --> https://rdrr.io/bioc/limma/src/R/voom.R --> log2-counts-per-million
  # y <- t(log2(t(counts+0.5)/(lib.size+1)*1e6))
  # y <- normalizeBetweenArrays(y,method=normalize.method) #from https://rdrr.io/bioc/limma/src/R/norm.R --> in normalizeBetweenArrays() it seems if object=matrix, quantile normalization
  # out$E <- y #out represents 'voom.out' object in this script
  voom.out <- voom(geneExpr) #Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observation-level weights. The data are then ready for linear modelling.
  counts_log_norm <- voom.out$E
  
  ## 4. log2 of weighted counts by number of cells per donor
  md.ss <- metadata[c('assignment','bare_barcode_lane')]
  md.ss %>%
    count(!!sym(aggregates)) %>%
    as.data.frame() -> cells_by_donor
  w <- cells_by_donor$n # weight
  counts_w.t <- t(counts)/w
  counts_weighted <- t(counts_w.t)
  counts_log_weighted <- log2(counts_weighted+0.5)
  # counts_log_weighted <- log2(counts_weighted+1)
  
  ## output
  pseudobulk_list <- list(counts = counts,
                          counts_log = counts_log,
                          counts_log_norm = counts_log_norm,
                          counts_weighted = counts_weighted,
                          counts_log_weighted = counts_log_weighted)
  out <- list(pseudobulk_LimmaVoom = pseudobulk_list,
              geneExpr_LimmaVoom = geneExpr,
              geneExpr_DESeq2 = aggregate_countMatrix_expr.deseq2, 
              metadata = aggregate_metadata)
  return(out)
}

#################### Set Variables and load Data #################### 
# Input/Output directories
in.dir <- paste0(opt$in_dir, '/')
out.dir <- paste0(opt$out_dir, '/')

# Dataset
# opt$dataset <- 'v2'
# opt$dataset <- 'v3'

# Cell type
# opt$cell_type <- 'CD8_TEM'

# Selection of new cells
# opt$genes <- 'data/genes_B3GAT1.tab' #B3GAT1
# opt$cell_type <- 'CD8_TEMpos'
# opt$cell_type <- 'CD8_TEMneg' 
if(!is.null(opt$genes)){
  genes_fn <- paste0(opt$genes)
  genes <- read.table(genes_fn)$V1
  genes_tag <- paste(genes, collapse='_')
  in.dir <- paste0(in.dir, genes_tag, '/')
  out.dir <- paste0(out.dir, genes_tag, '/')
}

# Input file
so.fn <- paste0(in.dir, '/', opt$dataset, '/',
                opt$cell_level, '/', opt$cell_type, '.so.rds')

# Output directory
out.dir <- paste0(out.dir, '/', opt$dataset, '/',
                  opt$cell_level, '/', opt$cell_type, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Report
print('############################')
print(paste0('Dataset: ', opt$dataset))
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Input file: ', so.fn))
print(paste0('Output directory: ', out.dir))
print('############################')

# Read seurat object
if(file.exists(so.fn)){
  print(paste0('Reading Seurat object: ', so.fn))
  system.time(pbmc <- readRDS(so.fn))
}else{
  err_in <- paste0('The file does not exist: ', so.fn)
  stop(err_in)
}

# # Downsampling (testing)
# set.seed(123)
# cells_to_subset <- sample(colnames(pbmc), 5000)
# pbmc <- subset(pbmc, cells = cells_to_subset)
# pbmc@meta.data <- droplevels(pbmc@meta.data)

#################### Convert sc matrix to pseudo-bulk matrix (by cell type - donor combination) #################### 
# Apply function
system.time(res <- pseudobulk.func(pbmc)) # complete

# Save pseudobulk matrices
fn <- paste0(out.dir, 'pseudobulk_matrices.rds')
print(paste0('Saving pseudobulk matrices in: ',fn))
saveRDS(res, fn)