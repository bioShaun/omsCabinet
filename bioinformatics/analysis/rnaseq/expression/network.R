suppressMessages(library(WGCNA, quietly=TRUE))
suppressMessages(library(argparser))
suppressMessages(library(tibble))
suppressMessages(library(ggdendro))
options(stringsAsFactors = FALSE);
enableWGCNAThreads()


p <- arg_parser("perform WGCNA analysis")
p <- add_argument(p, '--exp_table', help = 'expression table.')
p <- add_argument(p, '--out_dir', help = 'output directory.')
p <- add_argument(p, '--name', help = 'output name.')
p <- add_argument(p, '--blocksize', help = 'output name.')
p <- add_argument(p, '--power', help = 'WGCNA power', default = NULL)
argv <- parse_args(p)

## get parameter
exp_table <- argv$exp_table
out_dir <- argv$out_dir
out_name <- argv$name
sample_order <- argv$sample_order
blocksize <- as.numeric(argv$blocksize)
power_num <- argv$power

## mk out_dir
if (! dir.exists(out_dir)) {
    dir.create(out_dir)
}

# get sample order by cluster order
cal_dendro <- function(exp_df) {
  data.mat <- as.matrix(exp_df)
  sample_cor = cor(data.mat, method='pearson', use='pairwise.complete.obs')
  sample_dist = as.dist(1-sample_cor)
  sample_dist = sample_dist *100
  fit_hc <- hclust(sample_dist)
  dendr <- dendro_data(fit_hc,type="rectangle")
  text.df <- label(dendr)
  return(text.df)
}


## Choose a set of soft-thresholding powers
powers = c(1:30)

## Call the network topology analysis function
myData <- read.delim(exp_table, check.names = F, row.names=1)
mylogData <- log2(myData + 1)
samples <- colnames(myData)
cluster_df <- cal_dendro(mylogData)
sample_order_v <- cluster_df$label
out_samples <- sample_order_v[sample_order_v %in% samples]
mylogData <- mylogData[, out_samples]

datExpr <- as.data.frame(t(mylogData));
out_prefix <- paste(out_dir, out_name, sep='/')

if (is.null(power_num)) {
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
cex1 = 0.9;

pdf(paste(out_prefix, 'connectivity.pdf', sep='.'), width=8, height=8)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

scale_fit <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
scale_fit_cut <- rev(seq(0.6, 0.9, 0.01))

power_num <- which(scale_fit == max(scale_fit[2:30]))
scale_cut <- max(scale_fit[2:30])
for (each_cut in scale_fit_cut) {
    power_num <- which(scale_fit >= each_cut)[1]
    scale_cut <- each_cut
    if (! is.na(power_num)) {
        break
    }
}

print(scale_fit)
print(paste('power', power_num, scale_fit[power_num], sep=':'))

pdf(paste(out_prefix, 'scale.pdf', sep='.'), width=8, height=8)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=scale_cut,col="red")
dev.off()
} else {
    power_num <- as.numeric(power_num)
}

print(power_num)
print(blocksize)

## Constructing the gene network and identifying modules 
net = blockwiseModules(datExpr, power = power_num, maxBlockSize = blocksize,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = paste(out_prefix, "TOM", sep='.'),
                       verbose = 3)

## Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
pdf(file = paste(out_prefix, "module_colors.pdf", sep='.'), width = 12, height = 9);
## Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

probes <- names(datExpr)
gene_module_df <- data.frame(gene_id=probes, module=mergedColors)
write.table(gene_module_df, file=paste(out_prefix, 'gene.module.txt', sep='.'), row.names=F, quote=F, sep='\t')

samples <- row.names(datExpr)
sample_num <- length(samples)
sample_matrix <- matrix(rep(0, sample_num*sample_num), nrow=sample_num)
sample_df <- data.frame(sample_matrix)
row.names(sample_df) <- samples
colnames(sample_df) <- samples

for (each in samples) {
    sample_df[each, each] = 1
}

#write.table(sample_df, file=paste(out_prefix, 'sample.matrix.txt', sep='.'), quote=F, sep='\t')


nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
out_samples <- sample_order_v[sample_order_v %in% samples]
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
mes_df <- data.frame(MEs)
row.names(mes_df) <- row.names(datExpr)
mes_out_df <- rownames_to_column(mes_df, var = 'sample_id')
write.table(mes_out_df, file=paste(out_prefix, 'eigengenes.txt', sep='.'),
            quote=F, sep='\t', row.names=F)
moduleTraitCor = cor(MEs, sample_df, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#moduleTraitPvalue <- moduleTraitPvalue[, out_samples]
#moduleTraitCor <- moduleTraitCor[, out_samples]
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf(file = paste(out_prefix, "ME.sample.heatmap.pdf", sep='.'), width = 40, height = 12);
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

## output table
moduleTraitPvalue_df <- data.frame(moduleTraitPvalue)
moduleTraitCor_df <- data.frame(moduleTraitCor)
moduleTraitCor_out_df <- rownames_to_column(moduleTraitCor_df, var = 'Modules')
moduleTraitPvalue_out_df <- rownames_to_column(moduleTraitPvalue_df, var = 'Modules')
write.table(moduleTraitCor_out_df,
            file=paste(out_prefix, 'ME.sample.cor.matrix.txt', sep='.'),
            quote=F, sep='\t', row.names=F)
write.table(moduleTraitPvalue_out_df,
            file=paste(out_prefix, 'ME.sample.pval.matrix.txt', sep='.'),
            quote=F, sep='\t', row.names=F)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = paste(out_name,"wgcna.RData",sep='.'))
