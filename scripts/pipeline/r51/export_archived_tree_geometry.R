suppressPackageStartupMessages(library(ape))
base <- '../code/result_raw/KO_panel_repair_20260926/word_revision/delivery/Source_Data/core65'
out <- 'NEE_revision_20260926/51_evidence_figures_20261001/host'
for (domain in c('bacteria','archaea')) {
 tree <- read.tree(file.path(base,paste0(domain,'_monophyletic_lineages.tree')))
 dep <- node.depth.edgelength(tree)
 nodes <- data.frame(node=seq_along(dep),depth=dep,label=c(tree$tip.label,rep('',tree$Nnode)))
 write.table(nodes,file.path(out,paste0(domain,'_tree_nodes.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 write.table(data.frame(parent=tree$edge[,1],child=tree$edge[,2]),file.path(out,paste0(domain,'_tree_edges.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
}
