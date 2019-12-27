library(igraph)
library(data.table)
suppressMessages(library(argparser))
library(dplyr)

# setwd('C:\\work\\博士课题\\结果\\Analysis_results\\embryo\\network')
# 
# cor_file <- 'filter.TUCP_TUCP.correlation.txt'
# module_file <- 'TUCP.gene.module.txt'

p <- arg_parser("exp network")
p <- add_argument(p, '--cor_file', help = 'pearson correlation file.')
p <- add_argument(p, '--module_file', help = 'gene module file.')
p <- add_argument(p, '--out_prefix', help = 'output prefix')
argv <- parse_args(p)

cor_file <- argv$cor_file
module_file <- argv$module_file
out_prefix <- argv$out_prefix

link <- fread(cor_file)
node <- fread(module_file)

link <- link[!duplicated(link$V3), ]
link <- filter(link, V1 %in% node$gene_id &
               V2 %in% node$gene_id)


net <- graph.data.frame(link, node, directed=F)
#net=delete.vertices(net,which(degree(net)<5))
net=delete.vertices(net,which(degree(net)<2))


V(net)$color <- V(net)$module
V(net)$size <- 2
V(net)$frame.color <- NA

pdf(paste(out_prefix, 'network.pdf', sep = '.'), width=8,
    height=8, onefile = FALSE)
plot(net,
     vertex.label=NA,
     edge.width=.5,
     edge.color='lightgrey')
dev.off()

png(paste(out_prefix, 'network.png', sep = '.'), width=8,
    height=8, res = 300, units = 'in')
plot(net,
     vertex.label=NA,
     edge.width=.5,
     edge.color='lightgrey')
dev.off()
