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
  make_option(c("--cell_type"), action="store", default=NULL, type='character',
              help="Cell types in each of the cell_level argument."),
  make_option(c("--cell_level"), action="store", default="cell_type", type='character',
              help="Azimuth l1 (predicted.celltype.l1) or l2 (cell_type)."),
  make_option(c("--phenotype"), action="store", default="CMV_status", type='character',
              help="Gender/Age/CMV_status"),
  make_option(c("--ds_file"), action="store", default="data/Oelen2022_chemistries.tab", type='character',
              help="File with datasets to integrate."),
  make_option(c("--genes"), action="store", default=NULL, type='character',
              help="Create a new subset of cells."),
  make_option(c("--cmv_metadata"), action="store", default="data/Metadata_matchLLD_CMV.tsv", type='character',
              help="CMV metadata"),
  make_option(c("--debug"), action="store", default=FALSE, type='logical',
              help="Run dream() by gene. In case the execution stops."),
  make_option(c("--random"), action="store", default="date", type='character',
              help="Date, lane or any variables."),
  make_option(c("--covs"), action="store", default='data/CMV_sc_2_DEA/covariates.tab', type='character',
              help="Covariates file."),
  make_option(c("--in_dir"), action="store", default='output/CMV_sc_2.1_DEA_Pseudobulk', type='character',
              help="Main directory"),
  make_option(c("--out_dir"), action="store", default='output/CMV_sc_2.2_DEA_LimmaDream', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(data.table))
shhh(library(lme4))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(tibble))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(MAST))
shhh(library(Matrix))
shhh(library(edgeR))
shhh(library(limma))
shhh(library(variancePartition))
shhh(library(BiocParallel))

#################### Define functions #################### 
# Read pseudobulk input data (geneExpr_LimmaVoom and metadata)
read_data <- function(dataset){
  pseudobulk_fn <- paste0(in.dir, dataset, '/', 
                          opt$cell_level, '/', opt$cell_type, '/pseudobulk_matrices.rds')
  if(file.exists(pseudobulk_fn)){
    print(paste0('Reading pseudobulk object file in: ',pseudobulk_fn))
    system.time(pseudobulk_matrices <- readRDS(pseudobulk_fn))
  }else{
    err_in <- paste0('Pseudobulk have not been performed for cell type: ', opt$cell_type, ' in the dataset: ', dataset,'. voom() needs at least 2 genes to fit a mean-variance trend.')
    print(err_in)
    pseudobulk_matrices <- NULL
  }
  return(pseudobulk_matrices)
}

# Limma dream
## accessory function (used in 'dreamer' function): getContrast, dream, extract (from dream package)
# comb <- combinations[[1]]
# outdir <- out_dir
# form_model <- form
# vobjDream_obj <- vobjDream
# metadata_obj <- metadata
# contrast <- contrast_var
# comb_param <- comb_parameter
# debug <- debug
getContrast_by_comb <- function(comb = NULL, outdir, form_model, vobjDream_obj, metadata_obj, contrast, comb_param = NULL, debug){
  # get contrast
  if(is.null(comb)){
    print('Continous variable...')
    L = getContrast(vobjDream_obj, form_model, metadata_obj, contrast)
    comb <- contrast
  }else{
    print('Categorical variable...')
    metadata_obj[[contrast]] <- as.factor(metadata_obj[[contrast]])
    if(comb_param){
      # not tested --> try it by simulating data with another level in Gender
      print('The contrast variable has >2 levels...')
      L = getContrast(vobjDream_obj, form_model, metadata_obj, c(comb[[1]], comb[[2]]))
    }else{
      print('The contrast variable has <=2 levels...')
      comb_levels <- gsub(contrast, '', comb)
      metadata_obj[[contrast]] <- factor(metadata_obj[[contrast]], 
                                         levels=comb_levels)
      L = getContrast(vobjDream_obj, form_model, metadata_obj, comb[[2]])
    }
  }
  
  if(debug){
    print('Dream() alternative: subsetting the object by gene...')
    # Alternative: to skip potential errors
    ## Dream() for each gene -using try() to detect the error-
    dream_list <- list()
    for(j in 1:nrow(vobjDream_obj)){
      g <- rownames(vobjDream_obj)[j]
      print(paste0(j, ': ', g))
      fit = try(dream(vobjDream_obj[j,], form_model, metadata_obj, L))
      limma_res <- NA
      g_report <- TRUE
      if(class(fit)!='try-error'){
        fit = eBayes(fit)
        limma_res <- topTable(fit, coef=c('L1'), number=length(fit$F.p.value))
        limma_res <- limma_res %>%
          data.frame() %>%
          rownames_to_column(var="gene") %>%
          arrange(adj.P.Val)
      }else{
        fit <- NA
        g_report <- FALSE
      }
      names(g_report) <- g
      dream_list[[g]] <- list(fit = fit, 
                              top_table = limma_res,
                              g_report = g_report)
    }
    
    ## Get only the genes that do not give an error
    g_passed <- unlist(unname(lapply(dream_list, function(x) x$g_report)))
    g_out.names <- names(which(g_passed==FALSE))
    dream_list.passed <- dream_list[g_passed]
    g_out.tag <- paste(g_out.names, collapse = ', ')
    print(paste0('Genes with an error in dream(): ', g_out.tag))
    
    ## Get the toptable (data frame)
    limma_res <- do.call("rbind", lapply(dream_list.passed, function(x) x$top_table))
    limma_res$adj.P.Val <- p.adjust(limma_res$P.Value, method = "BH")
    # limma_res <- limma_res[order(limma_res$adj.P.Val),]
    limma_res <- limma_res[order(limma_res$P.Value),]
    
    ## Get the fit (list)
    fit <- lapply(dream_list.passed, function(x) x$fit)
    
  }else{
    print('Dream() default: using the whole object...')
    # Regular (but we can encounter an error --> example with gene 336 (RPL11): opt$phenotype <- 'Gender',  geneExpr <- geneExpr[501:1000,] as an example) --> Alternative
    # fit contrast
    fit = dream(vobjDream_obj, form_model, metadata_obj, L)
    fit = eBayes(fit)
    
    # grab the exact fit
    limma_res <- topTable(fit, coef=c('L1'), number=length(fit$F.p.value))
    limma_res <- limma_res %>%
      data.frame() %>%
      rownames_to_column(var="gene") %>%
      arrange(adj.P.Val)
    limma_res <- limma_res[order(limma_res$P.Value),]
  }
  
  # create list of variables
  vars_list <- list()
  vars_list[['form']] <- form_model
  vars_list[['vobjDream']] <- vobjDream_obj
  vars_list[['metadata']] <- metadata_obj
  vars_list[['L']] <- L
  vars_list[['fit']] <- fit #regular (fit dataframe for the whole object with all genes), alternative (fit list with only the genes without an error)
  
  # save/write the result
  comb_tag <- paste(comb,collapse='_')
  limma_res.fn <- paste0(outdir, comb_tag)
  limma_res_tsv.fn <- paste0(limma_res.fn, '.tsv')
  write.table(limma_res, limma_res_tsv.fn, sep = '\t', row.names = T)
  limma_res_rds.fn <- paste0(limma_res.fn, '.rds')
  saveRDS(limma_res, limma_res_rds.fn)
  vars.fn <- paste0(limma_res.fn, '.vars.rds')
  saveRDS(vars_list, vars.fn)
  
  return(limma_res)
}

## main function
# ge_dge = geneExpr
# covariates = covs.df
# contrast_var = opt$phenotype
# random_var = opt$random
# metadata = aggregate_metadata
# out_dir = out.dir
# debug = opt$debug
dreamer <- function(ge_dge, covariates, contrast_var, random_var, metadata, out_dir, debug){
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(4, "SOCK", progressbar=TRUE)
  register(param)
  
  # The variable to be tested must be a fixed effect
  ## defining the form
  fixed_var <- covariates[covariates$type=='fixed' & !covariates$covariate%in%contrast_var,]$covariate
  if(contrast_var%in%c('Age_cat', 'Age_cat_all')){
    contrast.var_out <- 'Age'
    fixed_var <- fixed_var[!fixed_var%in%contrast.var_out]
  }
  contrast_fixed.fmla <- paste(c(contrast_var,fixed_var),collapse='+')
  form_vars <- contrast_fixed.fmla
  contrast.var_out <- contrast_var
  if(!is.null(random_var)){
    print(random_var)
    random_var <- c(random_var, 'dataset')
    random_effects.fmla <- paste(paste0('(1|',random_var,')'),collapse='+')
    form_vars <- paste(c(contrast_fixed.fmla,random_effects.fmla), collapse='+')
  }
  form_vars <- paste0('~',form_vars)
  # form_vars <- "~Age+Gender+(1|date)+(1|dataset)" #testing
  form <- as.formula(form_vars)
  print(paste0('Fitting lmer: ',form_vars))
  print(paste0('Testing: ',contrast_var))
  
  # estimate weights using linear mixed model of dream
  vobjDream = voomWithDreamWeights(ge_dge, form, metadata, span = 0.1)
  vobjDream_max <- apply(vobjDream, 1, max)
  summary(vobjDream_max)
  head(sort(vobjDream_max,decreasing=T))
  
  # Fit the dream model on each gene
  # By default, uses the Satterthwaite approximation for the hypothesis test
  #fit = dream( vobjDream, form, metadata )
  
  ## use a contrast matrix
  combinations <- NULL
  if(is.factor(metadata[[contrast_var]])){
    contrast_levels <- levels(metadata[[contrast_var]])
    combinations <- split(combn(contrast_levels,2),  col(combn(contrast_levels,2)))
    combinations <- lapply(combinations, function(x) paste0(contrast_var,x))
    combinations <- unname(combinations)
    names(combinations) <- unlist(lapply(combinations, function(x) paste(x,collapse='_')))
    comb_parameter <- ifelse(nlevels(metadata[[contrast_var]])>2, TRUE, FALSE)
    # testing
    # combinations <- c(combinations, list(c("GenderF", "GenderM")))
  }
  
  ## apply getContrast_by_comb function
  if(!is.null(combinations)){
    getContrast.list <- lapply(combinations, function(i) getContrast_by_comb(comb = i,
                                                                             outdir = out_dir,
                                                                             form_model = form,
                                                                             vobjDream_obj = vobjDream,
                                                                             metadata_obj = metadata, 
                                                                             contrast = contrast_var, 
                                                                             comb_param = comb_parameter,
                                                                             debug = debug))
  }else{
    getContrast <- getContrast_by_comb(outdir = out_dir,
                                       form_model = form,
                                       vobjDream_obj = vobjDream,
                                       metadata_obj = metadata, 
                                       contrast = contrast_var,
                                       debug = debug)
    getContrast.list <- list()
    getContrast.list[[contrast_var]] <- getContrast
  }
  
  ## save
  suffix_fn <- paste(names(getContrast.list),collapse = '.')
  getContrast.list_fn <- paste0(out_dir, suffix_fn, '.rds')
  print(paste0('Saving DE results: ',getContrast.list_fn))
  saveRDS(getContrast.list, getContrast.list_fn)
  return(getContrast.list)
}

#################### Set Variables and load Data #################### 
# Input/Output directories
in.dir <- paste0(opt$in_dir, '/')
out.dir <- paste0(opt$out_dir, '/')

# Datasets
datasets_fn <- opt$ds_file
datasets <- read.table(datasets_fn)$V1
datasets_tag <- paste(datasets, collapse='_')

# cell_level
# opt$cell_level <- 'predicted.celltype.l1'
# opt$cell_level <- 'cell_type'

# cell_type
# opt$cell_type <- 'CD8_TEM'

# Selection of new cells 
opt$genes <- 'data/genes_B3GAT1.tab' #B3GAT1
opt$cell_type <- 'CD8_TEMpos'
# opt$cell_type <- 'CD8_TEMneg' 
if(!is.null(opt$genes)){
  genes_fn <- paste0(opt$genes)
  genes <- read.table(genes_fn)$V1
  genes_tag <- paste(genes, collapse='_')
  in.dir <- paste0(in.dir, genes_tag, '/')
  out.dir <- paste0(out.dir, genes_tag, '/')
}

# Debug
# opt$debug <- TRUE

# Random variable
# opt$random <- 'date'
# opt$random_var <- 'lane'

# Phenotype
# opt$phenotype <- 'Gender'
# opt$phenotype <- 'Age' 
# opt$phenotype <- 'CMV_status'

# Covariates
covs_fn <- opt$covs
print(paste0('Reading covariates file in: ',covs_fn))
covs.df <- read.table(covs_fn, header = T)

## Output directory
out.dir <- paste0(out.dir, '/', datasets_tag, '/', 
                  opt$cell_level, '/', opt$cell_type, '/',
                  opt$random, '/', opt$phenotype, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Report
print('############################')
print(paste0('Datasets: ', datasets_tag))
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Phenotype: ', opt$phenotype))
print(paste0('Input directory: ', in.dir))
print(paste0('Output directory: ', out.dir))
print('############################')
cat('\n')

#################### Data pre-processing #################### 
# Read pseudobulk input data (geneExpr_LimmaVoom and metadata)
## Apply function
system.time(pb_list <- sapply(datasets, function(i) read_data(i), simplify = FALSE))
pb_list <- Filter(Negate(is.null), pb_list)

## Join data from different datasets
types <- c('metadata', 'geneExpr_LimmaVoom')
types_datasets.list <- sapply(types, function(i) lapply(pb_list, function(x) x[[i]]), simplify = FALSE)

# Metadata (metadata)
## Join
aggregate_metadata <- do.call("rbind", types_datasets.list$metadata)
aggregate_metadata[,c('dataset', 'donor')] <- str_split_fixed(rownames(aggregate_metadata), '\\.', 2)[,c(1,2)]
rownames(aggregate_metadata) <- aggregate_metadata$donor
aggregate_metadata$dataset <- as.factor(aggregate_metadata$dataset)
aggregate_metadata$donor <- as.factor(aggregate_metadata$donor)
aggregate_metadata$Gender <- as.factor(aggregate_metadata$Gender)
aggregate_metadata$Age_cat <- ifelse(aggregate_metadata$Age<=40, 'Y',
                                     ifelse(aggregate_metadata$Age>=60, 'O', 'M')) 
aggregate_metadata$Age_cat_all <- ifelse(aggregate_metadata$Age<=40, 'Y', 'O')
aggregate_metadata$Age_cat <- as.factor(aggregate_metadata$Age_cat)
aggregate_metadata$Age_cat_all <- as.factor(aggregate_metadata$Age_cat_all)
if(opt$phenotype=='Age_cat'){
  aggregate_metadata <- droplevels(aggregate_metadata[!aggregate_metadata$Age_cat=='M',])
}
aggregate_metadata$LLDEEP_ID <- rownames(aggregate_metadata)

## Add CMV metadata
### Read data
cmv_fn <- opt$cmv_metadata
cmv_df <- read.delim(cmv_fn, check.names = FALSE)

### Subset by sc-donors
sc_donors <- unique(aggregate_metadata$LLDEEP_ID)
print(paste0('# of sc donors: ', length(sc_donors)))
cmv_df.sc <- droplevels(cmv_df[cmv_df$LLDEEP_ID%in%sc_donors,])
print(paste0('# of match donors: ', nrow(cmv_df.sc)))
cmv_df.sc <- cmv_df.sc[!is.na(cmv_df.sc$CMV_Baseline) | !is.na(cmv_df.sc$CMV_Followup),]
print(paste0('# of match donors (without NAs in CMV_Baseline and CMV_Followup): ', nrow(cmv_df.sc)))

### CMV status
cmv_df.sc$CMV_status <- ifelse(!is.na(cmv_df.sc$CMV_Followup), cmv_df.sc$CMV_Followup, cmv_df.sc$CMV_Baseline)
cmv_df.sc$Age_diff_SCvsLLD <- ifelse(!is.na(cmv_df.sc$CMV_Followup), cmv_df.sc$Age_diff_SCvsLLD2, cmv_df.sc$Age_diff_SCvsLLD1)
table(cmv_df.sc$CMV_status)
summary(cmv_df.sc$Age_diff_SCvsLLD)

### Add CMV metadata to seurat metadata
aggregate_metadata <- merge(aggregate_metadata, cmv_df.sc, by = 'LLDEEP_ID')
aggregate_metadata$CMV_status <- as.factor(as.character(aggregate_metadata$CMV_status))
rownames(aggregate_metadata) <- aggregate_metadata$LLDEEP_ID
table(aggregate_metadata$dataset, aggregate_metadata$CMV_status)

# Pseudobulk-gene expression profiles (geneExpr_LimmaVoom)
## Library size normalization: DGEList and calcNormFactors
counts_list <- lapply(types_datasets.list$geneExpr_LimmaVoom, function(x) x[['counts']])
rnames_common <- Reduce("intersect",lapply(counts_list, rownames))
counts_list <- lapply(counts_list, function(x) x[rnames_common,])
geneExpr <- do.call("cbind", counts_list)
geneExpr <- geneExpr[,rownames(aggregate_metadata)]
geneExpr = DGEList(geneExpr)
geneExpr = calcNormFactors(geneExpr)

## Make sure the donor ids from the metadata (aggregate_metadata) matches the ones in the expression data (geneExpr)
if(!identical(colnames(geneExpr), rownames(aggregate_metadata))){
  print('Matching donor ids from the metadata and gene expression...')
  idx <- match(colnames(geneExpr), rownames(aggregate_metadata))
  aggregate_metadata <- aggregate_metadata[idx,]
}

## Set contrast levels order
if(opt$phenotype%in%c('Gender', 'Age_cat', 'Age_cat_all', 'CMV_status')){
  if(opt$phenotype=='CMV_status'){Group_order <- c('0','1')}
  if(opt$phenotype=='Gender'){Group_order <- c('M','F')}
  if(opt$phenotype%in%c('Age_cat', 'Age_cat_all')){Group_order <- c('Y','O')}
  aggregate_metadata[[opt$phenotype]] <- factor(aggregate_metadata[[opt$phenotype]],
                                                levels = Group_order)
}
table(aggregate_metadata$dataset, aggregate_metadata$CMV_status)

contrast_LRT <- opt$phenotype
contrast_tag <- opt$phenotype
if(opt$phenotype%in%c('Gender','Age_cat', 'Age_cat_all', 'CMV_status')){
  contrast_LRT <- paste0(contrast_LRT, Group_order[2])
  contrast_tag <- paste(paste0(contrast_tag, Group_order[1]), paste0(contrast_tag, Group_order[2]), sep = '_')
}

#################### Limma dream ####################
# Apply function
print('Gene expression dimensions:')
dim(geneExpr)
# geneExpr.all <- geneExpr
# geneExpr <- geneExpr[1:100,] #testing
dreamer.res <- dreamer(ge_dge = geneExpr, 
                       covariates = covs.df, 
                       contrast_var = opt$phenotype, 
                       random_var = opt$random, 
                       metadata = aggregate_metadata, 
                       out_dir = out.dir, 
                       debug = opt$debug)