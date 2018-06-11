library(dplyr)
suppressMessages(library(argparser))

options(stringsAsFactors = F)
p <- arg_parser("perform correlation analysis")
p <- add_argument(p, '--exp_table', help = 'expression table.')
p <- add_argument(p, '--gene_type', help = 'gene type file.')
p <- add_argument(p, '--method', help = 'correlation analysis method.')
p <- add_argument(p, '--sample_order', help = 'output sample order.')
argv <- parse_args(p)

exp_df <- read.delim(argv$exp_table, row.names = 1,com = "", check.names = F)

# filter gene do not expressed in any sample
exp_df <- exp_df[rowSums(exp_df)>0, ]
write.table(exp_df, file='test.txt', sep='\t', quote=F)

# log transform
exp_df <- log10(exp_df + 1)

# read gene type file
gene_type_df <- read.delim(argv$gene_type)

# select and reorder samples in expression table
exp_samples <- colnames(exp_df)
gene_order_df <- read.delim(argv$sample_order, header=F)
out_samples <- gene_order_df$V1[gene_order_df$V1 %in% exp_samples]
exp_df <- exp_df[, out_samples]

# get all type of genes in gene type file
gene_types <- unique(gene_type_df$gene_biotype)

for (each_type in gene_types) {
    each_type_genes <- filter(gene_type_df, gene_biotype == each_type)$gene_id
    each_type_exp_df <- exp_df[each_type_genes,]
    each_type_exp <- paste(each_type, 'tpm.txt', sep='.')
    if (!file.exists(each_type_exp)) {
        write.table(each_type_exp_df, file=each_type_exp, sep='\t', quote=F)
    }
    sample_cor <- cor(each_type_exp_df, method=argv$method, use='pairwise.complete.obs')
    sample_cor_df <- as.data.frame(sample_cor)
    sample_dist_df <- 1 - sample_cor_df
    sample_cor_df <- cbind(Sample = rownames(sample_cor_df), sample_cor_df)
    sample_dist_df <- cbind(Sample = rownames(sample_dist_df), sample_dist_df)
    out_name <- paste(each_type, argv$method, 'correlation.txt', sep='.')
    dist_out_name <- paste(each_type, argv$method, 'dist.txt', sep='.')
    write.table(sample_cor_df, file=out_name, sep='\t', quote=F,
                row.names = F)
    write.table(sample_dist_df, file=dist_out_name, sep='\t', quote=F,
                row.names = F)
    
}
