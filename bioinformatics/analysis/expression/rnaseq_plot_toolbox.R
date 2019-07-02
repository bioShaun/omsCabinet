suppressMessages(library(omplotr))
suppressMessages(library(argparser))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))

options(stringsAsFactors = F)

p <- arg_parser('A rnaseq stats plot collection.')
p <- add_argument(p, '--exp_table',
                  help = 'Expression matrix for plot, required for [heatmap|PCA|cluster]',
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
p <- add_argument(p, '--group_mean_exp',
                  help = 'Plot using mean value of sample group.',
                  flag=TRUE)
p <- add_argument(p, '--test',
                  help = 'Analysis on first 1000 line of data.',
                  flag = TRUE)
p <- add_argument(p, '--cluster_cut_tree',
                  help = 'Cluster analysis cut tree height portion.',
                  default = 20)
p <- add_argument(p, '--out_prefix',
                  help = 'Plot output prefix.')

p <- add_argument(p, '--heatmap',
                  help = 'Plot heatmap.',
                  flag = TRUE)
p <- add_argument(p, '--heatmap_not_cluster_cols',
                  help = 'boolean values determining if columns should be clustered or hclust object.',
                  flag = TRUE)
p <- add_argument(p, '--heatmap_not_cluster_rows',
                  help = 'boolean values determining if rows should be clustered or hclust object.',
                  flag = TRUE)
p <- add_argument(p, '--heatmap_scale',
                  help = 'boolean values determining if rows should be clustered or hclust object.',
                  default = 'row')
p <- add_argument(p, '--pca',
                  help = 'Plot PCA scatterplot.',
                  flag = TRUE)
p <- add_argument(p, '--cluster',
                  help = 'Cluster analysis.',
                  flag = TRUE)
p <- add_argument(p, '--cluster_method',
                  help = 'Cluster method [hcluster or kmeans]',
                  default = 'hcluster')
p <- add_argument(p, '--kmeans_center',
                  help = 'kmeans center number.',
                  default = NULL)
argv <- parse_args(p)


MIN_CLUSTER_NUM = 10
MIN_CLUSTER_POR = 0.005
DIFF_HEATMAP_GENE = 40000

exp_table <- argv$exp_table
gene_list <- argv$gene_list
sample_list <- argv$sample_list
group_vs_sample <- argv$group_vs_sample
group_mean_exp <- argv$group_mean_exp
heatmap <- argv$heatmap
pca <- argv$pca
cluster <- argv$cluster
out_prefix <- argv$out_prefix
is_test <- argv$test
cluster_cut_tree <- argv$cluster_cut_tree
heatmap_cluster_cols <- !(argv$heatmap_not_cluster_cols)
heatmap_cluster_rows <- !(argv$heatmap_not_cluster_rows)
heatmap_scale <- argv$heatmap_scale
cluster_method <- argv$cluster_method
kmeans_center <- argv$kmeans_center



# out_prefix <- '203_4_vs_CKWT_3'


# cluster <- F
# heatmap <- TRUE
# pca <- F
# 
# heatmap_cluster_cols <- F
# heatmap_cluster_rows <- T
# heatmap_scale <- 'none'

# exp_table <- 'WT&58_expression.txt'
# gene_list <- NA
# sample_list <- NA
# group_vs_sample <- NA
# group_mean_exp <- FALSE
# cluster_method <- 'kmeans'
# kmeans_center <- 4
# is_test <- T
# out_prefix <- 'test'
# cluster <- T
# heatmap <- F
# pca <- F

if ( heatmap | pca | cluster){
  check_input(exp_table, 'Expression table is needed!')
  exp_df <- load_exp_file(exp_table, 
                          gene_list,
                          sample_list)
  exp_df <- test_data(exp_df, is_test)
  if (group_mean_exp) {
    print('Group sample expression by mean.')
    check_input(group_vs_sample, 'Calculate group mean expression need group information.')
    exp_df <- exp_by_group(exp_df, group_vs_sample)
  }
  sample_inf <- label_sample(exp_df, 
                             group_vs_sample,
                             sample_list,
                             group_mean_exp)
  if (heatmap) {
    print('Heatmap plot.')
    heatmap_prefix <- paste(out_prefix, 'heatmap', sep = '.')
    om_heatmap(exp_df, sample_inf, out_prefix=heatmap_prefix,
               scale=heatmap_scale,
               cluster_rows=heatmap_cluster_rows,
               cluster_cols=heatmap_cluster_cols)
  }
  if (pca) {
    print('PCA plot.')
    pca_prefix <- paste(out_prefix, 'pca', sep = '.')
    om_pca_plot(exp_df, sample_inf, out_prefix=pca_prefix)
  }
  if (cluster) {
    print('Cluster plot.')
    cluster_plot(exp_df, sample_inf, out_prefix, 
                 method=cluster_method, kmeans_center=kmeans_center,
                 cluster_cut_tree=cluster_cut_tree)
  }
}
