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
  make_option(c("--cell_level"), action="store", default=NULL, type='character',
              help="Azimuth l1 (predicted.celltype.l1) or l2 (cell_type)."),
  make_option(c("--interaction"), action="store", default=FALSE, type='logical',
              help="Interaction"),
  make_option(c("--genes"), action="store", default=NULL, type='character',
              help="Create a new subset of cells."),
  make_option(c("--cmv_metadata"), action="store", default="data/Metadata_matchLLD_CMV.tsv", type='character',
              help="CMV metadata"),
  make_option(c("--covs"), action="store", default='data/CMV_sc_1_CODA/covariates.tab', type='character',
              help="Covariates file."),
  make_option(c("--ds_file"), action="store", default="data/Oelen2022_chemistries.tab", type='character',
              help="File with datasets to integrate."),
  make_option(c("--in_dir"), action="store", default='data/CMV_sc_1_CODA/Azimuth_celltype_proportions', type='character',
              help="Main directory"),
  make_option(c("--out_dir"), action="store", default='output/CMV_sc_1_CODA', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(data.table))
shhh(library(tibble))
shhh(library(lmerTest))
shhh(library(broom))
shhh(library(broom.mixed))
shhh(library(ggplot2))
shhh(library(qvalue))

#################### Define functions #################### 
# Data description
## main function
bp_func <- function(df, fill_vec, width_var, out_dir, height_var = 3.5){
  title_var <- paste0(opt$dataset)
  
  p <- ggplot(df, aes(x=celltype, fill=CMV_status)) + 
    geom_bar(position="dodge", stat="count") + 
    theme_bw() +
    ylab(NULL) +
    ggtitle(paste0(title_var)) + 
    scale_fill_manual(values=fill_vec) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          strip.text = element_text(face="bold"),
          axis.text.x = element_text(angle = 90),
          plot.title = element_text(hjust=0.5, face="bold"),
          legend.title = element_text(hjust=0.5)) +
    facet_grid(Sex ~ .)
  
  p.fn <- paste0(out_dir,'Sex.CMV_status.bp.png')
  print(paste0('Saving barplot (sex - CMV status): ',p.fn))
  ggsave(p.fn, p, width = width_var, height = height_var)
}

# Adjust p-value
## accessory function (used in 'pval_adj.func()) --> https://github.com/StoreyLab/qvalue/blob/master/R/qvalue_trunc.R (https://support.bioconductor.org/p/124923/)
qvalue_truncp <- function(p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, pi0 = NULL, ...) {
  # Argument checks
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  } else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
    stop("'fdr.level' must be in (0, 1].")
  }
  p <- p / max(p)
  # Calculate pi0 estimate
  if (is.null(pi0)) {
    pi0s <- pi0est(p, ...)
  } else {
    if (pi0 > 0 && pi0 <= 1)  {
      pi0s = list()
      pi0s$pi0 = pi0
    } else {
      stop("pi0 is not (0,1]")
    }
  }
  
  # Calculate q-value estimates
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  if (pfdr) {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m / (i * (1 - (1 - p[o]) ^ m))))[ro]
  } else {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m /i ))[ro]
  }
  qvals_out[rm_na] <- qvals
  # Calculate local FDR estimates
  if (lfdr.out) {
    lfdr <- lfdr(p = p, pi0 = pi0s$pi0, ...)
    lfdr_out[rm_na] <- lfdr
  } else {
    lfdr_out <- NULL
  }
  
  # Return results
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level,
                   significant = (qvals <= fdr.level),
                   pi0.lambda = pi0s$pi0.lambda, lambda = pi0s$lambda,
                   pi0.smooth = pi0s$pi0.smooth)
  } else {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, pi0.lambda = pi0s$pi0.lambda,
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  class(retval) <- "qvalue"
  return(retval)
}
is.error <- function(x) inherits(x, "try-error")

## main function
pval_adj.func <- function(x, pval_col = 'p.value'){
  pval_raw <- x[[pval_col]]
  if(all(is.na(pval_raw))){
    x$pval_adj.holm <- NA
    x$pval_adj.fdr <- NA
    x$pval_adj.qvalue <- NA
  }else{
    x$pval_adj.holm <- p.adjust(pval_raw)
    x$pval_adj.fdr <- p.adjust(pval_raw, method = "fdr")
    qvals <- try(qvalue(pval_raw)$qvalues)
    if(is.error(qvals)){
      print('Default qvalue() is giving an erorr... Trying qvalue_truncp()')
      qvals <- try(qvalue_truncp(pval_raw)$qvalues)
      if(is.error(qvals)){
        print('Adapted qvalue_truncp() is giving an erorr... Trying a qvalue() with lambda=0.5')
        qvals <- try(qvalue(pval_raw, lambda=0.5)$qvalues)
        if(is.error(qvals)){
          print('qvalue() with lambda=0.5 is giving an error... Setting to NA')
          qvals <- rep(NA, nrow(x))
        }
      }
    }else{
      qvals <- qvalue(pval_raw)$qvalues
    }
    x$pval_adj.qvalue <- qvals
    
    # check
    pval_cnames <- c('pval_adj.holm', 'pval_adj.fdr', 'pval_adj.qvalue')
    pval_check <- sapply(pval_cnames, function(i) table(x[[i]]<=0.05), simplify = FALSE)
    print(pval_check)
    pval_cnames.df <- as.data.frame(combn(pval_cnames, 2))
    pval_cnames.list <- sapply(pval_cnames.df, function(x) as.character(x), simplify = FALSE)
    pval_cnames.cor.list <- lapply(pval_cnames.list, function(pval_vec){
      pval_tag <- paste(pval_vec, collapse = ' VS. ')
      print(pval_tag)
      cor_adj <- try(cor(x[[pval_vec[1]]], x[[pval_vec[2]]], use = "complete.obs"))
      if(is.error(cor_adj)){
        print('cor() is giving an erorr...')
        cor_adj <- NULL
      }
      print(cor_adj)
      return(cor_adj)
    })
  }
  cat('\n')
  return(x)
}

# Transform data (CLR)
## accessory functions
Geom_mean <- function(x){
  exp(mean(log(x)))
}
CLR <- function(D){
  log2(D / Geom_mean(D))
}

## main function
clr_func <- function(df){
  # Filter out cell types based on minimum nCells/donor: Long to wide --> rows = cell_type, cols = donor, value = n
  ## Check nCells/donor
  DF_n <- reshape2::dcast(df, donor ~ celltype, value.var = "n")
  rownames(DF_n) <- DF_n[,1]
  DF_n <- DF_n[,-1]
  DF_n <- t(DF_n)
  DF_n[is.na(DF_n)] <- 0
  
  ## Define filter: min 5 donors with at least 5 cells/donor
  nDonors_filt.byCT <- apply(DF_n, 1, function(x){sum(x>5)})
  CT_filt <- names(nDonors_filt.byCT[nDonors_filt.byCT>=5])
  
  # Long to wide --> rows = cell_type, cols = donor, value = freq
  DF_proportions <- reshape2::dcast(df, donor ~ celltype, value.var = "freq")
  rownames(DF_proportions) <- DF_proportions[,1]
  DF_proportions <- DF_proportions[,-1]
  DF_proportions <- t(DF_proportions)
  DF_proportions[is.na(DF_proportions)] <- 0
  
  # Add pseudocount
  pseudocount <- 1/3*min(DF_proportions[!DF_proportions==0])
  
  # CLR
  DF_proportions %>% 
    apply(2, function(x){CLR(x+pseudocount)}) -> DF_proportions.clr
  
  # # Check
  # apply(DF_proportions, 1, summary)
  # apply(DF_proportions.clr, 1, summary)
  
  ## Apply previously defined filter
  DF_proportions.clr.filt <- DF_proportions.clr[rownames(DF_proportions.clr)%in%CT_filt,]
  
  return(DF_proportions.clr.filt)
}

# Test per cell type (lmer/lm)
## main function
lm_by_ct <- function(cell_type, covs_df, df, md, interaction){
  print(cell_type)
  # Data
  vec_i <- df[rownames(df)==cell_type,]
  df_i <- as.data.frame(vec_i)
  colnames(df_i) <- 'freq'
  df_i$donor <- rownames(df_i)
  df_i <- merge(df_i, md, by = 'donor')
  rownames(df_i) <- df_i$donor
  df_i <- df_i[,-1]
  Group_order.CMV_status <- c('0','1')
  Group_order.Sex <- c('M','F')
  if('Sex'%in%colnames(df_i)){
    df_i[['Sex']] <- factor(df_i[['Sex']],
                            levels = Group_order.Sex)
  }
  if('CMV_status'%in%colnames(df_i)){
    df_i[['CMV_status']] <- factor(df_i[['CMV_status']],
                                   levels = Group_order.CMV_status)
  }
  
  # Formula
  covs_fixed <- covs_df[covs_df$type=='fixed',]$covariate
  covs_random <- covs_df[covs_df$type=='random',]$covariate
  if(length(covs_fixed)>0){
    fixed_fmla <- paste(covs_fixed,collapse='+')
  }
  random_fmla <- NULL
  if(length(covs_random)>0){
    random_fmla <- paste(paste0('(1|',covs_random,')'),collapse='+')
  }
  fmla <- paste(c(fixed_fmla, random_fmla), collapse = '+')
  if(interaction){
    interaction_fmla <- paste(vars_interaction, collapse = ':')
    fmla <- paste0(c(fmla, interaction_fmla), collapse = '+')
  }
  fmla <- paste0('freq ~ ',fmla)
  form <- as.formula(fmla)
  print(paste0('Fitting lmer: ',fmla))
  
  # lmer/lm
  if(!is.null(random_fmla)){
    # fit model
    print('lmer...')
    mod <-  lmerTest::lmer(form, data = df_i)
    
    # tidy model
    tidy_mod <- broom.mixed::tidy(mod, conf.int = TRUE, effects = "fixed")
    
  }else{
    # fit model
    print('lm...')
    mod <- lm(form, data = df_i)
    
    # tidy model
    tidy_mod <- broom::tidy(mod, conf.int = TRUE)
  }
  cat('\n')
  
  # tidy to dataframe
  cnames <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
  tidy_mod <- tidy_mod[,which(colnames(tidy_mod)%in%cnames)]
  tidy_mod.df <- as.data.frame(tidy_mod)
  tidy_mod.df$celltype <- cell_type
  
  return(tidy_mod.df)
}

# Dotplots: Output from the test per cell type (lmer/lm)
## wo/ facets
dp_func <- function(th_var, th = 0.05, cols_vec, alpha_vec, df, height_var, width_var, out_dir){
  print(th_var)
  df$signif <- ifelse(df[[th_var]]<=th, 'ss', 'ns')
  title_var <- paste0(opt$dataset, ' (', th_var, ' - ', as.character(th), ')')
  
  p <- ggplot(df, aes(x=estimate, y=celltype.label, color=direction, alpha=signif)) + 
    geom_point() +
    theme_bw() +
    ylab(NULL) +
    ggtitle(paste0(title_var)) + 
    geom_pointrange(aes(xmin=conf.low, xmax=conf.high), position=position_dodge(width=0.2), fatten = .25) +
    geom_vline(xintercept=0, linetype = "dashed") +
    scale_color_manual(values=cols_vec) +
    scale_alpha_manual(values=alpha_vec) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          strip.text = element_text(face="bold"),
          plot.title = element_text(hjust=0.5, face="bold"),
          legend.title = element_text(hjust=0.5)) +
    facet_grid(. ~ contrast,
               scales = 'free')
  
  suffix <- paste0(th_var, '_', as.character(th))
  p.fn <- paste0(out_dir, suffix,'.estimates_by_contrast.png')
  print(paste0('Saving dotplot + conf.int in: ',p.fn))
  ggsave(p.fn, p, width = width_var, height = height_var)
}

## w/ facets
dp_func.facets <- function(th_var, th = 0.05, cols_vec, alpha_vec, df, height_var, width_var, out_dir){
  print(th_var)
  df$signif <- ifelse(df[[th_var]]<=th, 'ss', 'ns')
  title_var <- paste0(opt$dataset, ' (', th_var, ' - ', as.character(th), ')')
  
  p <- ggplot(df, aes(x=estimate, y=celltype.label_short, color=direction, alpha=signif)) + 
    geom_point() +
    theme_bw() +
    ylab(NULL) +
    ggtitle(paste0(title_var)) + 
    geom_pointrange(aes(xmin=conf.low, xmax=conf.high), position=position_dodge(width=0.2), fatten = .25) +
    geom_vline(xintercept=0, linetype = "dashed") +
    scale_color_manual(values=cols_vec) +
    scale_alpha_manual(values=alpha_vec) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          strip.text = element_text(face="bold"),
          strip.text.y = element_text(angle=0),
          plot.title = element_text(hjust=0.5, face="bold"),
          legend.title = element_text(hjust=0.5)) +
    facet_grid(broad_celltype ~ contrast,
               scales = 'free')
  
  suffix <- paste0(th_var, '_', as.character(th))
  p.fn <- paste0(out_dir, suffix,'.estimates_by_contrast.facets.png')
  print(paste0('Saving dotplot + conf.int in: ',p.fn))
  ggsave(p.fn, p, width = width_var, height = height_var)
}

#################### Set Variables and load Data #################### 
# Interaction
# opt$interaction <- TRUE
# opt$covs <- 'data/CMV_sc_1_CODA/covariates.Sex_interaction.tab'
# opt$covs <- 'data/CMV_sc_1_CODA/covariates.Age_interaction.tab'

# Input directory
in.dir <- paste0(opt$in_dir, '/')

# Cell level
opt$cell_level <- 'cell_type'
# opt$cell_level <- 'predicted.celltype.l1'

# Datasets
datasets_fn <- opt$ds_file
datasets <- read.table(datasets_fn)$V1
datasets_tag <- paste(datasets, collapse='_')

# Covariates
covs_fn <- opt$covs
print(paste0('Reading covariates file in: ',covs_fn))
covs.df <- read.table(covs_fn, header = TRUE)

# Selection of new cells
opt$genes <- 'data/genes_B3GAT1.tab' #B3GAT1
# opt$genes <- 'data/genes_KLRC2.tab' #KLRC2
if(!is.null(opt$genes)){
  genes_fn <- paste0(opt$genes)
  genes <- read.table(genes_fn)$V1
  genes_tag <- paste(genes, collapse='_')
  in.dir <- paste0(in.dir, genes_tag, '/')
}

# Output directory
out.dir <- paste0(opt$out_dir, '/', datasets_tag, '/')
if(!is.null(opt$genes)){
  out.dir <- paste0(out.dir, genes_tag, '/')
}
out.dir <- paste0(out.dir, opt$cell_level, '/')
if(opt$interaction){
  vars_interaction <- covs.df[!is.na(covs.df$interaction) & covs.df$interaction=='interaction',]$covariate
  vars_interaction.tag <- paste(vars_interaction, collapse='.')
  out.dir <- paste0(out.dir, 'interaction/', vars_interaction.tag, '/')
  }
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Report
print('############################')
print(paste0('Datasets: ', datasets_tag))
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Interaction: ', as.character(opt$interaction)))
print(paste0('Input directory: ', in.dir))
print(paste0('Output directory: ', out.dir))
print('############################')
cat('\n')

# Proportions
## read data
in.dir <- paste0(in.dir, opt$cell_level, '/')
props_fn <- paste0(in.dir, 'proportions.rds')
print(paste0('Reading proportions file in: ',props_fn))
props.df <- readRDS(props_fn)
props.df <- droplevels(props.df[props.df$dataset%in%datasets,])

## change cell types names
cell_types <- unique(props.df$cell_type)
cell_type.vec <- gsub(' ', '_', cell_types)
names(cell_type.vec) <- cell_types
props.df$celltype <- unname(cell_type.vec[props.df$cell_type])
props.df <- droplevels(props.df)

## donor metadata
donor_md <- unique(props.df[,c('donor', 'Gender', 'Age', 'date', 'lane', 'dataset')])
colnames(donor_md)[colnames(donor_md)=='Gender'] <- 'Sex'
donor_md$lane <- as.factor(donor_md$lane)
donor_md$date <- as.factor(donor_md$date)

# Summarise info by celltype
props.df %>%
  group_by(celltype, .drop = FALSE) %>%
  summarise(n = sum(n)) -> count.df
count.df %>%
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) %>% as.data.frame() -> celltype.stats.df
celltype.stats.df$celltype.label <- paste0(celltype.stats.df$celltype,
                                           ' (n=', celltype.stats.df$n, ';freq=', as.character(round(celltype.stats.df$freq,2)), ')')

# Prepare variables
# Prepare variables
contrast_vec <- c('SexF', 'Age', 'CMV_status1')
names(contrast_vec) <- c('Sex', 'Age', 'CMV_status')

# Cell metadata
cell.md_fn <- paste0(in.dir, 'metadata.rds')
print(paste0('Reading proportions file in: ',cell.md_fn))
cell.md <- readRDS(cell.md_fn)
cell.md <- droplevels(cell.md[cell.md$dataset%in%datasets,])
cell_level.vec <- c('predicted.celltype.l1', 'cell_type') 
names(cell_level.vec) <- c('l1.pred', 'l2.pred')
cnames <- c('barcode', 'donor', 'date', names(cell_level.vec))
colnames(cell.md)[colnames(cell.md)=='Gender'] <- 'Sex'
cell.md <- cell.md[,cnames]
colnames(cell.md)[which(colnames(cell.md)%in%names(cell_level.vec))] <- unname(cell_level.vec)
cnames <- c('barcode', 'donor', opt$cell_level)
cell.md <- cell.md[,cnames]
colnames(cell.md)[which(colnames(cell.md)==opt$cell_level)] <- 'cell_type'

# CMV metadata
## read data
cmv_df <- read.delim(opt$cmv_metadata, check.names = FALSE)

## subset by sc-donors
sc_donors <- unique(props.df$donor)
print(paste0('# of sc donors: ', length(sc_donors)))
cmv_df.sc <- droplevels(cmv_df[cmv_df$LLDEEP_ID%in%sc_donors,])
print(paste0('# of match donors: ', nrow(cmv_df.sc)))
cmv_df.sc <- cmv_df.sc[!is.na(cmv_df.sc$CMV_Baseline) | !is.na(cmv_df.sc$CMV_Followup),]
print(paste0('# of match donors (without NAs in CMV_Baseline and CMV_Followup): ', nrow(cmv_df.sc)))

## CMV status
cmv_df.sc$CMV_status <- ifelse(!is.na(cmv_df.sc$CMV_Followup), cmv_df.sc$CMV_Followup, cmv_df.sc$CMV_Baseline)
cmv_df.sc$Age_diff_SCvsLLD <- ifelse(!is.na(cmv_df.sc$CMV_Followup), cmv_df.sc$Age_diff_SCvsLLD2, cmv_df.sc$Age_diff_SCvsLLD1)
cmv_df.sc$CMV_status <- as.factor(as.character(cmv_df.sc$CMV_status))
table(cmv_df.sc$CMV_status)
summary(cmv_df.sc$Age_diff_SCvsLLD)

## add CMV metadata to seurat metadata
props.df <- merge(props.df, cmv_df.sc, by.x = 'donor', by.y = 'LLDEEP_ID')
length(unique(props.df$donor))

# nDonors per cell type and CMV status --> to later filter datasets - cell type pairs in the pseudobulk-DEA meta-analysis
props.df %>%
  group_by(celltype, CMV_status, .drop = FALSE) %>%
  dplyr::count() -> donor_count.df

donor_count.df %>%
  group_by(celltype) %>%
  mutate(freq = n / sum(n)) %>% 
  as.data.frame() -> donor.stats.df

celltype.stats.df_link <- celltype.stats.df
colnames(celltype.stats.df_link)[c(2,3)] <- c('nCells', 'freqCells')
colnames(donor.stats.df)[c(3,4)] <- c('nDonors', 'freqDonors')
donor_celltype.stats.df <- merge(celltype.stats.df_link, donor.stats.df, by = 'celltype')
donor_celltype.stats.df <- donor_celltype.stats.df[order(-donor_celltype.stats.df$nCells),]
donor_celltype.stats.fn <- paste0(out.dir, 'donor_celltype.stats.txt')
write.table(donor_celltype.stats.df, donor_celltype.stats.fn, row.names = FALSE, quote = FALSE, sep = '\t')

# Pick donor metadata
donor_md <- unique(props.df[,c('donor', 'Sex.SC', 'Age.SC', 'CMV_status', 'dataset', 'date')])
colnames(donor_md)[colnames(donor_md)%in%c('Sex.SC', 'Age.SC')] <- c('Sex','Age')
donor_md$Sex <- ifelse(donor_md$Sex=='Male', 'M',
                       ifelse(donor_md$Sex=='Female', 'F', donor_md$Sex))
donor_md %>% mutate_if(is.character, as.factor) -> donor_md
donor_md$Sex <- factor(donor_md$Sex,
                       c('M','F'))
donor_md$CMV_status <- factor(donor_md$CMV_status,
                              c('0','1'))

#################### Data description: Barplot (CMV_status - Sex - cell type) ####################
# Data
celltype_md <- unique(props.df[,c('cell_type', 'celltype', 'donor', 'Sex.SC', 'Age.SC', 'CMV_status')])
colnames(celltype_md)[colnames(celltype_md)%in%c('Sex.SC', 'Age.SC')] <- c('Sex','Age')
celltype_md$Sex <- ifelse(celltype_md$Sex=='Male', 'M',
                          ifelse(celltype_md$Sex=='Female', 'F', celltype_md$Sex))
celltype_md <- merge(celltype_md, cell.md, by = c('donor','cell_type'))
celltypes.order <- names(sort(table(celltype_md$celltype),decreasing=T))
celltype_md$celltype <- factor(celltype_md$celltype,
                               celltypes.order)
celltype_md$Sex <- factor(celltype_md$Sex,
                          c('M','F'))
celltype_md$CMV_status <- factor(celltype_md$CMV_status,
                                 c('0','1'))

# Variables
CMV_status.vec <- c('#bfbfbf','#4d4d4d')
names(CMV_status.vec) <- c('0','1')
width.var <- ifelse(opt$cell_level=='predicted.celltype.l1', 6, 8)

# Apply function
bp_res <- bp_func(df = celltype_md, 
                  fill_vec = CMV_status.vec, 
                  width_var = width.var, 
                  out_dir = out.dir) 

#################### Transform data (CLR) ####################
# Apply function
props_clr.df <- clr_func(props.df)

#################### Save intermediate files ####################
props.fn <- paste0(out.dir, 'props.rds')
saveRDS(props.df, props.fn)
props_clr.fn <- paste0(out.dir, 'props_clr.rds')
saveRDS(props_clr.df, props_clr.fn)
donor_md.fn <- paste0(out.dir, 'donor_md.rds')
saveRDS(donor_md, donor_md.fn)

#################### Test per cell type (lmer/lm) ####################
# Prepare variables
celltypes.vec <- sort(apply(props_clr.df, 1, mean),decreasing=T)
celltypes <- names(celltypes.vec)

# Apply function
tidy_mod.list <- sapply(celltypes, function(i) lm_by_ct(cell_type = i,
                                                        covs_df = covs.df,
                                                        df = props_clr.df,
                                                        md = donor_md, 
                                                        interaction = opt$interaction), 
                        simplify = FALSE)
tidy_mod.df <- do.call("rbind", tidy_mod.list)

# Split by phenotype and compute FDR
tidy_mod.by_phe <- split(tidy_mod.df, tidy_mod.df$term)
tidy_mod.by_phe <- lapply(tidy_mod.by_phe, pval_adj.func)

# Get only interesting contrasts
if(opt$interaction){
  vars_interaction <- covs.df[!is.na(covs.df$interaction) & covs.df$interaction=='interaction',]$covariate
  vars_interaction.dict <- c('SexF', 'Age', 'CMV_status1')
  names(vars_interaction.dict) <- c('Sex', 'Age', 'CMV_status')
  vars_interaction.tag <- vars_interaction.dict[names(vars_interaction.dict)%in%vars_interaction]
  contrast_interaction.vec <- paste(unname(vars_interaction.tag), collapse=':')
  names(contrast_interaction.vec) <-  paste(names(vars_interaction.tag), collapse=':')
  contrast_interaction.vec_names <- names(contrast_interaction.vec)
  contrast.vec_names <- names(contrast_vec)
  contrast_vec <- c(unname(contrast_vec), unname(contrast_interaction.vec))
  names(contrast_vec) <- c(contrast.vec_names, contrast_interaction.vec_names) 
}
phe_stats.list <- tidy_mod.by_phe[names(tidy_mod.by_phe)%in%unname(contrast_vec)]
lapply(phe_stats.list, function(x) x[x$p.value<=0.05,]) #check
lapply(phe_stats.list, function(x) x[x$pval_adj.fdr<=0.05,]) #check
lapply(phe_stats.list, function(x) x[x$pval_adj.holm<=0.05,]) #check
phe_stats.list.fn <- paste0(out.dir, 'phenotype_stats.rds')
saveRDS(phe_stats.list, phe_stats.list.fn)

# List to DF
phe_stats.df <- do.call("rbind",phe_stats.list)
phe_stats.df$contrast <- str_split_fixed(rownames(phe_stats.df),'\\.',2)[,1]
phe_stats.df$direction <- ifelse(phe_stats.df$estimate>0, 'pos', 'neg')
contrast_vec <- c('Age', 'Sex', 'CMV_status', 'Sex:CMV_status', 'Age:CMV_status')
names(contrast_vec) <- c('Age', 'SexF', 'CMV_status1', 'SexF:CMV_status1', 'Age:CMV_status1')
phe_stats.df$contrast <- unname(contrast_vec[phe_stats.df$contrast])
phe_stats.df <- merge(phe_stats.df, celltype.stats.df[,c('celltype','celltype.label')], by = 'celltype')
base_levels <- c('CMV_status', 'Age', 'Sex')
all_levels <- unname(contrast_vec)[unname(contrast_vec)%in%unique(phe_stats.df$contrast)]
in_levels <- c(base_levels, setdiff(all_levels, base_levels))
phe_stats.df$contrast <- factor(phe_stats.df$contrast,
                                levels = in_levels)
celltype.label <- unique(celltype.stats.df$celltype.label)
phe_stats.df$celltype.label <- factor(phe_stats.df$celltype.label,
                                      rev(celltype.label))

#################### Dotplots: Output from the test per cell type (lmer/lm) ####################
# Dotplot with estimate+conf.int
## x-axis: estimate
## y-axis: cell types
## facets (x): phenotype
## color: positive/negative
## alpha: significance (p.value or fdr)

# wo/ facets
## Variables
direction.vec <- c('#cc0000','#003399')
names(direction.vec) <- c('pos','neg')
signif.vec <- c(0.2,1)
names(signif.vec) <- c('ns','ss')
height.var <- ifelse(opt$cell_level=='predicted.celltype.l1', 3.5, 6.5)
width.var <- ifelse(opt$cell_level=='predicted.celltype.l1', 8, 9.5)
th_var.vec <- c('p.value', 'pval_adj.fdr', 'pval_adj.holm')

## Apply function
dp_res <- sapply(th_var.vec, function(i) dp_func(th_var = i,
                                                 th = 0.05, 
                                                 cols_vec = direction.vec, 
                                                 alpha_vec = signif.vec, 
                                                 df = phe_stats.df, 
                                                 height_var = height.var, 
                                                 width_var = width.var, 
                                                 out_dir = out.dir), 
                 simplify = FALSE)

# w/ facets
## Variables
phe_stats.df$broad_celltype <- gsub('pos|neg', '', phe_stats.df$celltype)
props_broad.df <- props.df
props_broad.df$broad_celltype <- gsub('pos|neg', '', props_broad.df$celltype)
props_broad.df %>%
  group_by(broad_celltype, .drop = FALSE) %>%
  summarise(n = sum(n)) -> count_broad.df
count_broad.df %>%
  mutate(freq = n / sum(n)) %>% 
  arrange(desc(freq)) %>% as.data.frame() -> celltype.stats_broad.df
phe_stats.df$broad_celltype <- factor(phe_stats.df$broad_celltype,
                                      levels = unique(celltype.stats_broad.df$broad_celltype))
phe_stats.df$celltype.label_short <- gsub(".*(?=pos|neg)", "", phe_stats.df$celltype.label, perl = TRUE)

## Apply function
dp_res.facets <- sapply(th_var.vec, function(i) dp_func.facets(th_var = i,
                                                               th = 0.05, 
                                                               cols_vec = direction.vec, 
                                                               alpha_vec = signif.vec, 
                                                               df = phe_stats.df, 
                                                               height_var = height.var, 
                                                               width_var = width.var, 
                                                               out_dir = out.dir), 
                        simplify = FALSE)
