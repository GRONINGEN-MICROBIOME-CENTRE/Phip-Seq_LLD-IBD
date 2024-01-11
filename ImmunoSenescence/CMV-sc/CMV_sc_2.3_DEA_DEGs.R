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
  make_option(c("--cell_level"), action="store", default="cell_type", type='character',
              help="Azimuth l1 (predicted.celltype.l1) or l2 (cell_type)."),
  make_option(c("--phenotypes"), action="store", default="data/Age_Gender_CMV.phenotypes.tab", type='character',
              help="Gender/Age/Age_cat/Age_cat_all"),
  make_option(c("--ds_file"), action="store", default="data/Oelen2022_chemistries.tab", type='character',
              help="File with datasets to integrate."),
  make_option(c("--genes"), action="store", default=NULL, type='character',
              help="Create a new subset of cells."),
  make_option(c("--random"), action="store", default="date", type='character',
              help="Date, lane or any variables."),
  make_option(c("--props_dir"), action="store", default="data/CMV_sc_1_CODA/Azimuth_celltype_proportions", type='character',
              help="Proportions directory"),
  make_option(c("--in_dir"), action="store", default="output/CMV_sc_2.2_DEA_LimmaDream", type='character',
              help="Input main directory"),
  make_option(c("--out_dir"), action="store", default="output/CMV_sc_2.3_DEA_DEGs", type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(data.table))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(ggrepel))
shhh(library(RColorBrewer))
shhh(library(UpSetR))
shhh(library(grid))
shhh(library(openxlsx))
shhh(library(parallel))
shhh(library(collapse))
shhh(library(viridis))
shhh(library(rrvgo))

#################### Define functions #################### 
# Read input data
## main function
read_data <- function(phe, ct){
  print(phe)
  print(ct)
  contrast_LRT <- phe
  contrast_tag <- phe
  if(phe%in%c('Gender','Age_cat', 'Age_cat_all', 'CMV_status')){
    if(phe=='CMV_status'){Group_order <- c('0','1')}
    if(phe=='Gender'){Group_order <- c('M','F')}
    if(phe%in%c('Age_cat', 'Age_cat_all')){Group_order <- c('Y','O')}
    contrast_LRT <- paste0(contrast_LRT, Group_order[2])
    contrast_tag <- paste(paste0(contrast_tag, Group_order[1]), paste0(contrast_tag, Group_order[2]), sep = '_')
  }
  in_fn <- paste0(in.dir, ct, '/', opt$random, '/', phe, '/', contrast_tag, '.rds')
  if(file.exists(in_fn)){
    print(paste0('Reading DEA output in: ', in_fn))
    dea_df <- readRDS(in_fn)[[contrast_tag]]
    dea_df$cell_type <- ct
    dea_df$phenotype <- phe
  }else{
    dea_df <- NULL
  }
  return(dea_df)
}

# Upset plot function
## main function
upset.func <- function(phe, dir, degs_list){
  print(phe)
  print(dir)
  
  width.var <- 7
  if(phe=='CMV_status' & dir=='down'){width.var <- 12}
  
  # Pick DEGs
  degs_l <- degs_list[[phe]][[dir]]
  lapply(degs_l, length)
  
  # Upset plot
  degs.l <- degs_l[unname(unlist(lapply(degs_l, function(x) length(x)>0)))]
  if(length(degs.l)>1){
    upset_p <- upset(fromList(degs.l), nsets = length(degs.l), 
                     sets = names(degs.l), 
                     point.size = 3,
                     order.by = "freq", 
                     text.scale = 2)
    p_title <- paste0(phe, ' > ', dir, ' > ', opt$cell_level)
    upset.fn <- paste0(out.dir, phe, '.', dir, '.upset.pdf')
    print(paste0('Saving UpSet plot in: ',upset.fn))
    pdf(file=upset.fn, 
        onefile=FALSE,
        width = width.var)
    print(upset_p)
    grid.text(paste0(p_title),
              x = 0.65, y=0.95, 
              gp=gpar(fontsize=12, fontface="bold"))
    dev.off()
  }else{
    print(paste0('There are no DEGs in more than one cell type... Only in: ', 
                 names(degs.l), ' (', paste(unname(unlist(degs.l)),collapse=','), ')'))
    upset_p <- NULL
  }
  return(upset_p)
}

# Save files for FEA (webgestalt -website-): phe - direction - ct (all) combination
## main function
save_files_fea <- function(phe, dir, ct, degs_list, dea_df){
  print(phe)
  print(dir)
  print(ct)
  degs <- degs_list[[phe]][[dir]][[ct]]
  if(!is.null(degs)){
    # save DEGs
    degs_fn <- paste0(out.dir, phe, '.', dir, '.', ct, '.txt')
    print(degs_fn)
    write.table(as.data.frame(degs), degs_fn,
                row.names = FALSE, col.names = FALSE, 
                quote = FALSE, sep = '\t') 
    
    # save background
    bg_fn <- paste0(out.dir, phe, '.', ct, '.bg.txt')
    if(!file.exists(bg_fn)){
      print(bg_fn)
      dea <- dea_df[dea_df$phenotype==phe & dea_df$cell_type==ct,]
      genes <- dea$gene
      write.table(as.data.frame(genes), bg_fn,
                  row.names = FALSE, col.names = FALSE, 
                  quote = FALSE, sep = '\t') 
    }
  }
}

# Save files for FEA (webgestalt -website-): phe - direction - ct (intersect) combination
## main function
save_intersect_files_fea <- function(dir, phe = 'CMV_status', cts =  c('CD4_CTL', 'CD8_TEM'), degs_list, dea_df){
  cts_tag <- paste(cts, collapse = '.')
  
  print(phe)
  print(dir)
  print(cts_tag)
  degs.phe_dir <- degs_list[[phe]][[dir]]
  degs.phe_dir.cts <- degs.phe_dir[names(degs.phe_dir)%in%cts]
  degs <- Reduce("intersect", degs.phe_dir.cts)
  
  if(!is.null(degs)){
    # save DEGs
    degs_fn <- paste0(out.dir, phe, '.', dir, '.', cts_tag, '.txt')
    print(degs_fn)
    write.table(as.data.frame(degs), degs_fn,
                row.names = FALSE, col.names = FALSE, 
                quote = FALSE, sep = '\t') 
    
    # save background
    bg_fn <- paste0(out.dir, phe, '.', cts_tag, '.bg.txt')
    if(!file.exists(bg_fn)){
      print(bg_fn)
      dea <- dea_df[dea_df$phenotype==phe & dea_df$cell_type%in%cts,]
      dea.cts <- split(dea, dea$cell_type)
      genes.cts <- lapply(dea.cts, function(x) x$gene)
      genes <- Reduce("intersect", genes.cts)
      write.table(as.data.frame(genes), bg_fn,
                  row.names = FALSE, col.names = FALSE, 
                  quote = FALSE, sep = '\t') 
    }
  }
}

# Barplots
## all
bp.func <- function(var_y, df_plot){
  print(var_y)
  p <- ggplot(data=df_plot, aes(x=cell_type, y=.data[[var_y]])) + 
    geom_bar(stat="identity", position=position_dodge()) + 
    theme_bw() + 
    ylab(paste0(var_y)) + 
    xlab(NULL) + 
    ggtitle(paste0(opt$cell_level)) + 
    theme(axis.text.x = element_text(angle=90, hjust = 1),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.title = element_blank())+
    facet_grid(phenotype ~ .,
               scales = 'free')
  p.fn <- paste0(out.dir, var_y, '.barplot.png')
  ggsave(p.fn, p, width = width_var, height = height_var)
}

## logFC threshold
bp_th.func <- function(phe, var_y, df_plot, dir_vec, th_vec){
  print(phe)
  print(var_y)
  df <- droplevels(df_plot[df_plot$phenotype==phe,])
  ct.order.in <- ct.order[ct.order%in%unique(df$cell_type)]
  df$cell_type <- factor(df$cell_type,
                         levels = ct.order.in)
  df$th <- factor(df$th,
                  levels = names(th_vec))
  df$direction <- factor(df$direction,
                         levels = names(dir_vec))
  var_y.tag <- var_y
  if(var_y=='n'){var_y.tag <- 'Number of DEGs'}
  
  width.var <- 6
  if(phe!='CMV_status'){width.var <- 12}
  
  p <- ggplot(data=df, aes(x=th, y=.data[[var_y]])) + 
    geom_bar(aes(fill=direction, alpha=th), stat="identity", position=position_dodge()) + 
    theme_bw() + 
    ylab(paste0(var_y.tag)) + 
    xlab('abs(logFC) threshold') + 
    ggtitle(paste0(opt$cell_level, ' > ', phe)) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          strip.text = element_text(face="bold"),
          axis.text.x = element_text(face="bold"),
          legend.title = element_blank(),
          plot.title = element_text(hjust=0.5, face="bold"))+
    facet_grid(direction ~ cell_type,
               scales = 'free')+
    scale_fill_manual(values=dir_vec,
                      guide = "none") + 
    scale_alpha_manual(values=th_vec,
                       guide = "none")
  p.fn <- paste0(out.dir, phe, '.', var_y, '.barplot.png')
  ggsave(p.fn, p, width = width.var, height = 4)
}

################################## Set Variables and load Data ################################## 
# Input/Output directories
in.dir <- paste0(opt$in_dir, '/')
out.dir <- paste0(opt$out_dir, '/')
props.dir <- paste0(opt$props_dir, '/')

# phenotypes
phenotypes_fn <- opt$phenotypes
phenotypes <- read.table(phenotypes_fn)$V1
phenotypes_tag <- paste(phenotypes, collapse='_')

# cell_level
opt$cell_level <- 'cell_type'
width_var <- 6.5
height_var <-  4.5
height_var.upset <- 5
if(opt$cell_level=='cell_type'){
  width_var <- 9
  height_var <-  5
  height_var.upset <- 10
}

# Datasets
datasets_fn <- opt$ds_file
datasets <- read.table(datasets_fn)$V1
datasets_tag <- paste(datasets, collapse='_')

# Set common variables
cl_mean.order.fn <- 'data/cl_mean.order.rds'
cl_mean.order <- readRDS(cl_mean.order.fn)
ct.order <- cl_mean.order[[opt$cell_level]]

# Selection of new cells 
opt$genes <- 'data/genes_B3GAT1.tab' #B3GAT1
if(!is.null(opt$genes)){
  genes_fn <- paste0(opt$genes)
  genes <- read.table(genes_fn)$V1
  genes_tag <- paste(genes, collapse='_')
  in.dir <- paste0(in.dir, genes_tag, '/')
  out.dir <- paste0(out.dir, genes_tag, '/')
  props.dir <- paste0(props.dir, genes_tag, '/')
  
  ct_fn <- 'data/Azimuth_l2.cell_type.CMV.tab'
  cts <- read.table(ct_fn)$V1
  cts.idx <- unname(unlist(sapply(ct.order, function(i) grep(i, cts), simplify = FALSE)))
  ct.order <- cts[cts.idx]
}

# Input dir
in.dir <- paste0(in.dir, '/', datasets_tag, '/', opt$cell_level, '/') 

# Output directory
out.dir <- paste0(out.dir, '/',  datasets_tag, '/', opt$cell_level, '/', 
                  opt$random, '/', phenotypes_tag, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Report
print('############################')
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Phenotypes: ', phenotypes_tag))
print(paste0('Datasets: ', datasets_tag))
print(paste0('Input directory: ', in.dir))
print(paste0('Output directory: ', out.dir))
print('############################')
cat('\n')

#################### Data pre-processing #################### 
# Read input data and define DEGs
dea_list <- sapply(phenotypes, function(i) 
  sapply(ct.order, function(j) read_data(phe = i, ct = j), simplify = FALSE), simplify = FALSE)
dea_list <- lapply(dea_list, function(phe_l){
  phe_l.filt <- Filter(Negate(is.null), phe_l)
  dea.df <- do.call("rbind", phe_l.filt)
  return(dea.df)
})
dea.df <- do.call("rbind", dea_list)
dea.df$ss <- ifelse(dea.df$adj.P.Val<=0.05, 'ss', 'ns')
dea.df <- dea.df[!is.na(dea.df$ss),]

# # Check the N of genes tested in CD8_TEMneg and CD8_TEMpos (and with CD4_CTL) ###
# l <- split(dea.df, dea.df$cell_type)
# genes_tested.list <- lapply(l, rownames)
# common_genes.list <- sapply(names(genes_tested.list), function(x) 
#   sapply(names(genes_tested.list), function(y) 
#     Reduce("intersect", genes_tested.list[names(genes_tested.list)%in%c(x,y)]), 
#     simplify=FALSE), simplify=FALSE)
# ncommon_genes.list <- lapply(common_genes.list, function(x) lapply(x, length))
# ncommon_genes <- lapply(lapply(ncommon_genes.list, unlist), function(x) sort(x, decreasing=T))
# ncommon_genes$CD8_TEMneg #CD8_TEMneg = 12199 --> CD8_TEMpos = 2033
# ncommon_genes$CD4_CTLneg #CD4_CTLneg = 3436 --> CD4_CTLpos = 1023

# Filter out cell types based on minimum nCells/donor: Long to wide --> rows = cell_type, cols = donor, value = n
## read data
props.dir <- paste0(props.dir, '/', opt$cell_level, '/')
props_fn <- paste0(props.dir, '/proportions.rds')
print(paste0('Reading proportions file in: ',props_fn))
props.df <- readRDS(props_fn)
props.df <- droplevels(props.df[props.df$dataset%in%datasets,])
cell_types <- unique(props.df$cell_type)
cell_type.vec <- gsub(' ', '_', cell_types)
names(cell_type.vec) <- cell_types
props.df$celltype <- unname(cell_type.vec[props.df$cell_type])
props.df <- droplevels(props.df)

## check nCells/donor
DF_n <- reshape2::dcast(props.df, donor ~ celltype, value.var = "n")
rownames(DF_n) <- DF_n[,1]
DF_n <- DF_n[,-1]
DF_n <- t(DF_n)
DF_n[is.na(DF_n)] <- 0

## define filter: min 5 donors with at least 5 cells/donor
nDonors_filt.byCT <- apply(DF_n, 1, function(x){sum(x>5)})
CT_filt <- names(nDonors_filt.byCT[nDonors_filt.byCT>=5])

# pick DEGs
dea.df <- droplevels(dea.df[dea.df$cell_type%in%CT_filt,])
dea.df$direction <- ifelse(dea.df$logFC>0, 'up', 'down')
degs.df <- droplevels(dea.df[dea.df$ss=='ss',])
degs.df %>% 
  group_by(phenotype, direction, cell_type) %>% 
  dplyr::count() %>% as.data.frame() -> degs.count
degs_by_phe.list <- split(degs.df, degs.df$phenotype)
degs_phe_dir.list <- lapply(degs_by_phe.list, function(phe) split(phe, phe$direction))
degs_phe_dir_ct.list <- lapply(degs_phe_dir.list, function(phe) 
  lapply(phe, function(dir) split(dir, dir$cell_type)))
degs_phe_dir_ct.list <- lapply(degs_phe_dir_ct.list, function(phe)
  lapply(phe, function(dir) 
    lapply(dir, function(ct) ct$gene)))
dir_vec <- c('up', 'down')
cts_vec <- unique(degs.df$cell_type)

# Save files for FEA (webgestalt -website-)
## by phe - direction - ct (all) combination
save_degs_file.out <- lapply(phenotypes, function(i)
  lapply(dir_vec, function(j)
    lapply(cts_vec, function(k) save_files_fea(phe = i, 
                                               dir = j, 
                                               ct = k, 
                                               degs_list = degs_phe_dir_ct.list,
                                               dea_df = dea.df))))

## by phe - dir - cts (intersect) combination
save_intersect_file.out <- lapply(dir_vec, function(i) save_intersect_files_fea(dir = i,
                                                                                degs_list = degs_phe_dir_ct.list,
                                                                                dea_df = dea.df))

# Save DEA and DEGs
dea_fn <- paste0(out.dir, 'dea.rds')
degs_fn <- paste0(out.dir, 'degs.rds')
saveRDS(dea.df, dea_fn)
saveRDS(degs.df, degs_fn)

#################### Plots #################### 
phenotypes <- phenotypes[phenotypes%in%unique(dea.df$phenotype)]

# Upset Plot per phenotype-cell type
upset_phe_ct.res <- sapply(phenotypes, function(i) 
  sapply(dir_vec, function(j) upset.func(phe = i, 
                                         dir = j,
                                         degs_list = degs_phe_dir_ct.list), simplify = FALSE), simplify = FALSE)

# Barplots
## Prepare data (all)
### Counts and Proportions (of 'ss' category over dataset-cell_type-ss combination)
dea.df %>%
  group_by(phenotype, direction, cell_type, ss, .drop=FALSE) %>%
  dplyr::count() -> dea_n_ss.df
dea_n_ss.df %>%
  group_by(phenotype, direction, cell_type, .drop=FALSE) %>%
  mutate(proportion = n / sum(n)) %>% as.data.frame() -> dea_prop.df

### Set common variables
df <- dea_prop.df
df$ss <- as.character(df$ss)
df$ss <- factor(df$ss,
                levels = c('ns', 'ss'))
ct.order_in <- ct.order[ct.order%in%unique(df$cell_type)]
df$cell_type <- factor(df$cell_type,
                       levels = ct.order_in)
df$phenotype <- factor(df$phenotype,
                       levels = phenotypes)
df$n_log10 <- log10(df$n+1)
var_y.vec <- c('n','n_log10','proportion')

### Pick only DEGs
dea_ss.df <- droplevels(df[!is.na(df$ss) & df$ss=='ss',])
dea_ss.df <- dea_ss.df[order(dea_ss.df$phenotype, -dea_ss.df$n),]

## Apply function
bp_degs.res <- lapply(var_y.vec, function(i) bp.func(var_y = i,
                                                     df_plot = dea_ss.df))

## Prepare data (only significant, and logFC threshold)
### Thresholds
degs.df$th <- ifelse(abs(degs.df$logFC)>=2, '[2, ...)', 
                     ifelse(abs(degs.df$logFC)>=1.5, '[1.5, 2)', 
                            ifelse(abs(degs.df$logFC)>=1, '[1, 1.5)', 
                                   ifelse(abs(degs.df$logFC)>=0.5, '[0.5, 1)', '[0, 0.5)'))))
### Counts and Proportions (of 'ss' category over dataset-cell_type-ss combination)
degs.df %>%
  group_by(phenotype, direction, cell_type, th, .drop=FALSE) %>%
  dplyr::count() -> degs_n_ss.df
degs_n_ss.df %>%
  group_by(phenotype, direction, cell_type, .drop=FALSE) %>%
  mutate(proportion = n / sum(n)) %>% as.data.frame() -> degs_prop.df
degs_prop.df$n_log10 <- log10(degs_prop.df$n+1)

### Set common variables
# th_alpha <- c('#ced4da', '#adb5bd', '#6c757d', '#495057', '#212529')
th_alpha <- c(0.25, 0.45, 0.65, 0.85, 1)
names(th_alpha) <- c('[0, 0.5)', 
                     '[0.5, 1)',
                     '[1, 1.5)',
                     '[1.5, 2)',
                     '[2, ...)')
dir_cols <- c('#CE3C28', '#053398')
names(dir_cols) <- c('up', 'down')
var_y.vec <- c('n','n_log10','proportion')

## Apply function
bp_degs_th.res <- lapply(phenotypes, function(i) 
  lapply(var_y.vec, function(j) bp_th.func(phe = i,
                                           var_y = j,
                                           df_plot = degs_prop.df, 
                                           dir_vec = dir_cols, 
                                           th_vec = th_alpha)))

################# FEA from WebGestalt #################
# Read data
## function
# dir <- i
# ct <- j
# rr <- reduncancy_red
# phe <- phenotype
read_fea <- function(dir, ct, rr, phe){
  print(dir)
  print(ct)
  
  fea_dir <- paste0(out.dir, phe, '.', dir, '.', ct, '/')
  # read fea
  fea_fn <- list.files(fea_dir, full.names=TRUE, pattern='enrichment_results_wg_result.+.txt')
  fea_df <- read.delim(fea_fn, header=TRUE)
  fea_df$direction <- dir
  fea_df$cell_type <- ct
  cnames <- c('direction', 'cell_type', 'geneSet', 'description', 'size', 'overlap', 'expect', 
              'enrichmentRatio', 'pValue', 'FDR')
  fea_df <- fea_df[,cnames]
  
  if(!is.null(rr)){
    print('Redundancy reduction (webgestalt)...')
    # filter
    fit_pat <- paste0('enriched_geneset_', rr, '_wg_result.+.txt')
    filt_fn <- list.files(fea_dir, full.names=TRUE, pattern=fit_pat)
    filt_gs <- read.delim(filt_fn, header=TRUE)[,1]
    fea_df <- fea_df[fea_df$geneSet%in%filt_gs,]
  }
  return(fea_df)
}

# Dotplot
## function
# df <- fea_rrvgo.df
# height_var <- height.var
# tag <- reduncancy_red
# phe <- phenotype
dp_func <- function(df, height_var, tag, phe){
  my_breaks <- seq(floor(min(df$enrichmentRatio)), ceiling(max(df$enrichmentRatio)), length.out=5)
  my_limits <- c(floor(min(df$enrichmentRatio)), ceiling(max(df$enrichmentRatio)))
  p <- ggplot(df, aes(x=cell_type, y = description, color = enrichmentRatio, size = log10_fdr)) + 
    geom_point() + 
    theme_bw() + 
    ylab(NULL) + 
    xlab(NULL) + 
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, size = 13, face = "bold"),
          axis.text.y = element_text(size = 12),
          strip.text.x = element_text(face="bold", size = 13),
          strip.text.y = element_text(face="bold", size = 13, angle = 0),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust=0.5, face="bold"),
          panel.spacing.x = unit(0.75, "lines"))+
    facet_grid(parentTerm ~ direction,
               scales = 'free', space = 'free') + 
    scale_color_viridis_c(name = "Enrichment Ratio",
                          breaks = my_breaks, 
                          limits = my_limits) + 
    guides(size=guide_legend(title="-log10(FDR)"))
  if(is.null(tag)){
    tag <- 'all'
  }
  p.fn <- paste0(out.dir, phe, '.', tag, '.fea.png')
  ggsave(p.fn, p, width = 12, height = height_var)
}

# Main function
# reduncancy_red = NULL
# height.var <- 12
# reduncancy_red <- 'wsc_topsets' #or ap_clusters
# height.var <- 6.5
# dir_vec = c('up', 'down')
# ct_vec = c('CD4_CTL', 'CD8_TEM', 'CD4_CTL.CD8_TEM')
# phenotype = 'CMV_status'
dp_fea <- function(reduncancy_red = NULL, height.var = 12, dir_vec = c('up', 'down') , ct_vec = c('CD4_CTL', 'CD8_TEM', 'CD4_CTL.CD8_TEM'), phenotype = 'CMV_status'){
  print(paste0('Phenotype: ', phenotype))
  print(paste0('Redundancy reduction: ', reduncancy_red))
  
  # Read data
  print('Read data...')
  # i <- dir_vec[1]
  # j <- ct_vec[1]
  fea_list <- sapply(dir_vec, function(i)
    sapply(ct_vec, function(j) read_fea(dir = i,
                                        ct = j,
                                        rr = reduncancy_red,
                                        phe = phenotype), simplify = FALSE), simplify = FALSE)
  fea_df <- do.call("rbind", lapply(fea_list, function(dir) do.call("rbind", dir)))
  fea_df.ss <- droplevels(fea_df[fea_df$FDR<=0.05,])
  fea_df.ss$direction <- factor(fea_df.ss$direction,
                                levels = dir_vec)
  cts.in <- ct_vec[ct_vec%in%unique(fea_df.ss$cell_type)]
  fea_df.ss$cell_type <- factor(fea_df.ss$cell_type,
                                levels = cts.in)
  pseudocount <- 1/3*min(fea_df.ss$FDR[!fea_df.ss$FDR==0])
  fea_df.ss$log10_fdr <- -log10(fea_df.ss$FDR+pseudocount)
  colnames(fea_df.ss)[colnames(fea_df.ss)=='size'] <- 'size_webgestalt'
  
  # Reduce GO terms
  go_df <- fea_df.ss
  go_terms <- unique(go_df$geneSet)
  
  ## only if reduceSimMatrix(scores=go_vec)
  # go_vec <- go_df$log10_fdr 
  # names(go_vec) <- go_df$geneSet
  
  ## simMatrix
  simMatrix <- calculateSimMatrix(go_terms, 
                                  orgdb="org.Hs.eg.db", 
                                  ont="BP", 
                                  method="Rel")
  
  ## reduceSimMatrix (explore options --> keep default)
  reducedTerms.uniqueness <- reduceSimMatrix(simMatrix,
                                             scores = 'uniqueness',
                                             threshold=0.9,
                                             orgdb="org.Hs.eg.db") #default
  # reducedTerms.size <- reduceSimMatrix(simMatrix,
  #                                      scores = 'size',
  #                                      threshold=0.9,
  #                                      orgdb="org.Hs.eg.db") #alternative
  # reducedTerms.sign <- reduceSimMatrix(simMatrix,
  #                                      scores = go_vec,
  #                                      threshold=0.9,
  #                                      orgdb="org.Hs.eg.db") #it makes no sense cause we're merging different FEA tests
  reducedTerms <- reducedTerms.uniqueness
  
  # Check reduction
  ## parentTerms
  # reducedTerms <- reducedTerms.uniqueness
  # reducedTerms <- reducedTerms.size
  # reducedTerms <- reducedTerms.sign
  l <- lapply(split(reducedTerms, reducedTerms$parentTerm), function(x) unique(x$term))
  length(l)
  parentTerm_n <- sort(table(reducedTerms$parentTerm),decreasing=T)
  ll <- l[names(parentTerm_n)]
  
  ## plotting
  reducedTerms <- reducedTerms.uniqueness
  score_var <- 'uniqueness'
  # reducedTerms <- reducedTerms.size
  # score_var <- 'size'
  # reducedTerms <- reducedTerms.sign
  # score_var <- 'sign'
  
  rr.tag <- reduncancy_red
  if(is.null(reduncancy_red)){
    rr.tag <- 'all'
  }
  
  ### heatmap
  hm <- heatmapPlot(simMatrix,
                    reducedTerms,
                    annotateParent=TRUE,
                    annotationLabel="parentTerm",
                    fontsize=6)
  hm.fn <- paste0(out.dir, rr.tag, '.', score_var, '.hm.pdf')
  pdf(file=hm.fn, 
      onefile=FALSE)
  print(hm)
  dev.off()
  
  ### scatterplot
  sp <- scatterPlot(simMatrix, reducedTerms)
  sp.fn <- paste0(out.dir, rr.tag, '.', score_var, '.sp.pdf')
  pdf(file=sp.fn, 
      onefile=FALSE)
  print(sp)
  dev.off()
  
  # Merge webgestalt + rrvgo reduction
  fea_rrvgo.df <- merge(fea_df.ss, reducedTerms, by.x = 'geneSet', by.y = 'go')
  fea_rrvgo.df$parentTerm <- factor(fea_rrvgo.df$parentTerm,
                                    levels = names(parentTerm_n)) 
  fea_rrvgo.df$term <- factor(fea_rrvgo.df$term,
                              levels = unname(unlist(ll)))
  
  # Dotplot
  print('Dotplot...')
  dp_all <- dp_func(df = fea_rrvgo.df,
                    height_var = height.var,
                    tag = reduncancy_red,
                    phe = phenotype)
  
  # Save dataframe
  fea_rrvgo.tmp <- fea_rrvgo.df[,c(3,2,1,4:ncol(fea_rrvgo.df))]
  fea_rrvgo.list.tmp <- split(fea_rrvgo.tmp, fea_rrvgo.df[,c('cell_type', 'direction')])
  fea_rrvgo.list <- fea_rrvgo.list.tmp[which(lapply(fea_rrvgo.list.tmp, nrow) != 0)]
  fea_rrvgo.list <- lapply(fea_rrvgo.list, function(x){
    xx <- x[order(x$FDR),]
    return(xx)
  })
  fea_rrvgo.fn <- paste0(out.dir, 'fea_rrvgo.rds')
  print(paste0('Saving FEA: ', fea_rrvgo.fn))
  saveRDS(fea_rrvgo.list, fea_rrvgo.fn)
  return(NULL)
}

## apply
### all terms
dp_fea.all <- dp_fea()

### reduced
dp_fea.red <- dp_fea(reduncancy_red = 'wsc_topsets', height.var = 6.5)

################# Exploration only for l2 - B3GAT1 #################
degs.l <- split(degs.df, degs.df$cell_type)
degs <- degs.l$CD8_TEMpos$gene
degs.ll <- lapply(degs.l, function(x) x[x$gene%in%degs,])
degs.ll <- lapply(degs.ll, function(x){rownames(x)<-NULL;return(x[,c(1:7)])})
degs.ll$CD8_TEMneg <- degs.ll$CD8_TEMneg[c(1,3,4,2),]
degs.ll$CD4_CTLneg <- degs.ll$CD4_CTLneg[c(3,1,2),]
link_list <- degs.ll[c(3,2,1)]
# link_list <- link_list[c(1,2)]
link_list <- sapply(names(link_list), function(i){
  x <- link_list[[i]]
  x$log10_fdr <- -log10(x$adj.P.Val)
  cnames <- paste0(colnames(x)[-1] ,'.', i)
  colnames(x)[-1] <- cnames
  return(x)
}, simplify = FALSE)
df_to_plot <- Reduce(function(x, y) merge(x, y, by = "gene", all.x = TRUE, all.y = TRUE), link_list)

# CD8_TEMpos vs. CD8TEM_neg
max_val <- round(max(c(max(df_to_plot$logFC.CD8_TEMneg), max(df_to_plot$logFC.CD8_TEMpos))),1)+0.2
p <- ggplot(data=df_to_plot, aes(x=logFC.CD8_TEMneg, y=logFC.CD8_TEMpos, label = gene)) + 
  geom_point(color='blue') + geom_text_repel() +  
  theme_bw() +
  geom_abline(slope = 1) + 
  xlim(0, max_val) + ylim(0, max_val) 
p_fn <- paste0(out.dir, 'CD8_TEMposCD8TEM_neg.logFC.pdf')
ggsave(p_fn, p, width = 5, height = 4)

min_val <- round(-log10(0.05),1)
max_val <- round(max(c(max(df_to_plot$log10_fdr.CD8_TEMneg), max(df_to_plot$log10_fdr.CD8_TEMpos))),1)+0.2
p <- ggplot(data=df_to_plot, aes(x=log10_fdr.CD8_TEMneg, y=log10_fdr.CD8_TEMpos, label = gene)) + 
  geom_point(color='blue') + geom_text_repel() +  
  theme_bw() +
  geom_abline(slope = 1) + 
  xlim(0, max_val) + ylim(0, max_val) 
p_fn <- paste0(out.dir, 'CD8_TEMposCD8TEM_neg.log10_fdr.pdf')
ggsave(p_fn, p, width = 5, height = 4)

# CD8_TEMpos vs. CD4_CTLneg
df_to_plot <- droplevels(df_to_plot[!is.na(df_to_plot$logFC.CD4_CTLneg),])
max_val <- round(max(c(max(df_to_plot$logFC.CD4_CTLneg), max(df_to_plot$logFC.CD8_TEMpos))),1)+0.2
p <- ggplot(data=df_to_plot, aes(x=logFC.CD4_CTLneg, y=logFC.CD8_TEMpos, label = gene)) + 
  geom_point(color='blue') + geom_text_repel() +  
  theme_bw() +
  geom_abline(slope = 1) + 
  xlim(0, max_val) + ylim(0, max_val) 
p_fn <- paste0(out.dir, 'CD8_TEMposCD4_CTLneg.logFC.pdf')
ggsave(p_fn, p, width = 5, height = 4)

min_val <- round(-log10(0.05),1)
max_val <- round(max(c(max(df_to_plot$log10_fdr.CD4_CTLneg), max(df_to_plot$log10_fdr.CD8_TEMpos))),1)+0.2
p <- ggplot(data=df_to_plot, aes(x=log10_fdr.CD4_CTLneg, y=log10_fdr.CD8_TEMpos, label = gene)) + 
  geom_point(color='blue') + geom_text_repel() +  
  theme_bw() +
  geom_abline(slope = 1) + 
  xlim(0, max_val) + ylim(0, max_val) 
p_fn <- paste0(out.dir, 'CD8_TEMposCD4_CTLneg.log10_fdr.pdf')
ggsave(p_fn, p, width = 5, height = 4)