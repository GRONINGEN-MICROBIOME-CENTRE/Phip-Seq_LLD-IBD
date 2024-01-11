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
  make_option(c("--phenotype"), action="store", default="CMV_status", type='character',
              help="Phenotype of interest (CMV_status)."),
  make_option(c("--phenotypes"), action="store", default="data/Age_Gender_CMV.phenotypes.tab", type='character',
              help="Phenotypes tested in the DEA."),
  make_option(c("--ds_file"), action="store", default="data/Oelen2022_chemistries.tab", type='character',
              help="File with datasets to integrate."),
  make_option(c("--random"), action="store", default="date", type='character',
              help="Date, lane or any variables."),
  make_option(c("--in_dir"), action="store", default="output/CMV_sc_2.3_DEA_DEGs", type='character',
              help="Input main directory"),
  make_option(c("--out_dir"), action="store", default="output/CMV_sc_2.4_DEA_FEA", type='character',
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
# Dotplot
## accessory function: read output from Webgestalt
read_fea <- function(dir, ct, rr, phe, in_dir){
  print(dir)
  print(ct)
  
  # read fea
  fea_dir <- paste0(in_dir, phe, '.', dir, '.', ct, '/')
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

## accessory function: dotplot
dp_func <- function(df, height_var, tag, score, phe, out_dir){
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
  p.fn <- paste0(out_dir, phe, '.', tag, '.', score, '.fea.png')
  ggsave(p.fn, p, width = 12, height = height_var)
}

## main function
dp_fea <- function(reduncancy_red = NULL, height.var = 12, dir_vec = c('up', 'down') , ct_vec = c('CD4_CTL', 'CD8_TEM', 'CD4_CTL.CD8_TEM'), phenotype, in_dir, out_dir){
  print(paste0('Phenotype: ', phenotype))
  print(paste0('Redundancy reduction: ', reduncancy_red))
  
  # Read data
  print('Read data...')
  fea_list <- sapply(dir_vec, function(i)
    sapply(ct_vec, function(j) read_fea(dir = i,
                                        ct = j,
                                        rr = reduncancy_red,
                                        phe = phenotype,
                                        in_dir = in_dir), simplify = FALSE), simplify = FALSE)
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
  
  ## reduceSimMatrix (explore options --> keep default 'uniqueness')
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
  hm.fn <- paste0(out_dir,  phenotype, '.', rr.tag, '.', score_var, '.hm.pdf')
  pdf(file=hm.fn, 
      onefile=FALSE)
  print(hm)
  dev.off()
  
  ### scatterplot
  sp <- scatterPlot(simMatrix, reducedTerms)
  sp.fn <- paste0(out_dir, phenotype, '.', rr.tag, '.', score_var, '.sp.pdf')
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
                    score = score_var,
                    phe = phenotype,
                    out_dir = out_dir)
  
  # Save dataframe
  fea_rrvgo.tmp <- fea_rrvgo.df[,c(3,2,1,4:ncol(fea_rrvgo.df))]
  fea_rrvgo.list.tmp <- split(fea_rrvgo.tmp, fea_rrvgo.df[,c('cell_type', 'direction')])
  fea_rrvgo.list <- fea_rrvgo.list.tmp[which(lapply(fea_rrvgo.list.tmp, nrow) != 0)]
  fea_rrvgo.list <- lapply(fea_rrvgo.list, function(x){
    xx <- x[order(x$FDR),]
    return(xx)
  })
  fea_rrvgo.fn <- paste0(out_dir,  phenotype, '.', rr.tag, '.', score_var, '.fea_rrvgo.rds')
  print(paste0('Saving FEA: ', fea_rrvgo.fn))
  saveRDS(fea_rrvgo.list, fea_rrvgo.fn)
  return(NULL)
}

################################## Set Variables and load Data ################################## 
# Input/Output directories
in.dir <- paste0(opt$in_dir, '/')
out.dir <- paste0(opt$out_dir, '/')

# phenotypes
phenotypes_fn <- opt$phenotypes
phenotypes <- read.table(phenotypes_fn)$V1
phenotypes_tag <- paste(phenotypes, collapse='_')

# Datasets
datasets_fn <- opt$ds_file
datasets <- read.table(datasets_fn)$V1
datasets_tag <- paste(datasets, collapse='_')

# Input dir
in.dir <- paste0(in.dir, '/', datasets_tag, '/', opt$cell_level, '/', opt$random, '/', phenotypes_tag, '/')

# Output directory
out.dir <- paste0(out.dir, '/',  datasets_tag, '/', opt$cell_level, '/', 
                  opt$random, '/', phenotypes_tag, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Report
print('############################')
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Phenotype: ', opt$phenotype))
print(paste0('Datasets: ', datasets_tag))
print(paste0('Input directory: ', in.dir))
print(paste0('Output directory: ', out.dir))
print('############################')
cat('\n')

#################### Plots #################### 
# Apply function
## all terms
dp_fea.all <- dp_fea(phenotype = opt$phenotype, 
                     in_dir = in.dir, 
                     out_dir = out.dir)

## reduced terms
dp_fea.red <- dp_fea(reduncancy_red = 'wsc_topsets', height.var = 6.5)