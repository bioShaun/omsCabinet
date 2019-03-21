suppressMessages(library(omplotr))
suppressMessages(library(argparser))

options(stringsAsFactors = F)

p <- arg_parser('A rnaseq stats plot collection.')
p <- add_argument(p, '--exp_table',
                  help = 'Expression matrix for plot, required for [heatmap|PCA]',
                  default = NULL)
p <- add_argument(p, '--gene_list',
                  help = 'Sellected genes list for analysis.',
                  default = NULL)
p <- add_argument(p, '--sample_list',
                  help = 'Sellected sample list for analysis.',
                  default = NULL)
p <- add_argument(p, '--group_vs_sample',
                  help = 'Group vs Sample file. <Group_id>tab<Sample_id>',
                  default = NULL)
p <- add_argument(p, '--out_prefix',
                  help = 'Plot output prefix.')

p <- add_argument(p, '--heatmap',
                  help = 'Plot heatmap.',
                  flag = TRUE)
p <- add_argument(p, '--pca',
                  help = 'Plot PCA scatterplot.',
                  flag = TRUE)
argv <- parse_args(p)


exp_table <- argv$exp_table
gene_list <- argv$gene_list
sample_list <- argv$sample_list
group_vs_sample <- argv$group_vs_sample
heatmap <- argv$heatmap
pca <- argv$pca
out_prefix <- argv$out_prefix

# exp_table <- 'Gene.tpm.xls'
# gene_list <- 'photo.genes'
# sample_list <- 'sample.list'
# group_vs_sample <- 'group_sample'
# out_prefix <- 'test'
# heatmap <- TRUE
# pca <- TRUE


check_input <- function(file_path, err_message=NULL) {
  if ( is.null(file_path) || ! file.exists(file_path)) {
    if (is.null(err_message)) {
      err_message <- paste(file_path, 'Not exist!')
    }
    stop(err_message)
  }
}


valid_input <- function() {
  print('Input file is valid.')
}


select_exp_data <- function(exp_obj, item_file, select='row') {
  if (! is.na(item_file)) {
    check_input(item_file)
    item_df <- read.delim(item_file, header = F)
    item_num <- dim(item_df)[1]
    if (select == 'row') {
      exp_obj <- dplyr::filter(exp_obj, target_id %in% item_df$V1)
      item_name <- 'targets'
    } else if (select == 'column') {
      exp_obj <- dplyr::select(exp_obj, c('target_id', item_df$V1))
      item_name <- 'samples'
    } else {
      stop('Wrong select parameter [row, column].')
    }
    print(paste('Select', item_num, item_name))    
  }
  return(exp_obj)
}


load_exp_file <- function(exp_file, genes, samples) {
  exp_obj <- data.table::fread(exp_file, check.names=F)
  total_target <- dim(exp_obj)[1]
  total_sample <- dim(exp_obj)[2] - 1
  print(paste('Total', total_target, 'targets,',
              total_sample, 'samples.'))
  colnames(exp_obj)[1] <- 'target_id'
  exp_obj <- select_exp_data(exp_obj, genes, 'row')
  exp_obj <- select_exp_data(exp_obj, samples, 'column')
  valid_input()
  exp_df <- data.frame(exp_obj, check.names=F)
  rownames(exp_df) <- exp_df$target_id
  exp_df <- exp_df[, -1]
  return(exp_df)
}


label_sample <- function(exp_df, group_vs_sample) {
  if (! is.na(group_vs_sample)) {
    sample_inf <- read.delim(group_vs_sample, header=F)
    colnames(sample_inf) <- c("condition", "sample")
  } else {
    sample_inf <- data.frame(condition=colnames(exp_df), 
                             sample=colnames(exp_df))
  }
  return(sample_inf)
}


if ( heatmap | pca ){
  check_input(exp_table, 'Expression table is needed!')
  exp_df <- load_exp_file(exp_table, 
                          gene_list,
                          sample_list)
  sample_inf <- label_sample(exp_df, 
                             group_vs_sample)
  if (heatmap) {
    print('Heatmap plot.')
    heatmap_prefix <- paste(out_prefix, 'heatmap', sep = '.')
    om_heatmap(exp_df, sample_inf, out_prefix=heatmap_prefix)
  }
  if (pca) {
    print('PCA plot.')
    pca_prefix <- paste(out_prefix, 'pca', sep = '.')
    om_pca_plot(exp_df, sample_inf, out_prefix=pca_prefix)
  }
}
