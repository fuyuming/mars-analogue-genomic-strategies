#!/usr/bin/env Rscript
# Contract: depict lineage-level marker breadth and KO clustering relative to
# the source-preserving null; do not claim historical adaptive convergence.
# Reuse existing collapsed-tree geometry; current student layout has A/B/C only.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
base <- "code/result_raw/KO_panel_repair_20260926/word_revision/delivery"
source <- readLines("code/207.make_ISME_figure2_lineage_mosaic.R")
source <- source[seq_len(which(grepl("^# d: root-robust",source))-1L)]
code <- paste(source,collapse="\n")
old <- 'rec <- fread(paths[["recurrence"]])'
new <- 'rec <- fread(recurrence_null_path)[null_model == "within_source_prevalence_preserving"]; rec[, variable_state := variable_trait]'
stopifnot(grepl(old,code,fixed=TRUE))
code <- gsub(old,new,code,fixed=TRUE)
old <- 'lineage_breadth <- fread(paths[["lineage_breadth"]])'
new <- 'lineage_breadth <- fread(corrected_breadth); setnames(lineage_breadth,"mean_module_breadth","mean_marker_breadth")'
stopifnot(grepl(old,code,fixed=TRUE));code<-gsub(old,new,code,fixed=TRUE)
code<-gsub('"Cold / PQC"','"Cold / lipid"',code,fixed=TRUE)
code<-gsub('"Dormancy"','"Sporulation"',code,fixed=TRUE)
code<-gsub('"Trace gas"','"H2/CO/CF"',code,fixed=TRUE)
Sys.setenv(ISME_UPPERCASE_PANEL_TAGS="1",ISME_PUBLICATION_TEXT_GATE="1")
branches<-commandArgs(trailingOnly=TRUE);if(!length(branches))branches<-c("core65","restored74")
for(branch in branches) {
  out<-file.path(base,"Figures_corrected",paste0("Figure2_",branch))
  sd<-file.path(base,"Source_Data",branch)
  e<-new.env(parent=globalenv())
  e$corrected_breadth<-file.path(sd,"Figure2B_lineage_module_breadth.csv")
  e$commandArgs<-function(trailingOnly=FALSE)c(out,file.path(sd,"recurrence/KO_recurrence_null_results.csv"))
  set.seed(20260927)
  branch_code<-if(branch=="restored74")gsub('"Sporulation"','"Spore / Rpf"',code,fixed=TRUE)else code
  eval(parse(text=branch_code),envir=e)
  fig<-e$p_a/e$p_b/e$p_c+plot_layout(heights=c(.68,3.55,1.15))
  prefix<-file.path(out,paste0("Figure2_",branch,"_corrected"))
  svglite::svglite(paste0(prefix,".svg"),width=183/25.4,height=215/25.4);print(fig);dev.off()
  cairo_pdf(paste0(prefix,".pdf"),width=183/25.4,height=215/25.4,family="Arial");print(fig);dev.off()
  ragg::agg_png(paste0(prefix,".png"),width=183,height=215,units="mm",res=300);print(fig);dev.off()
  ragg::agg_tiff(paste0(prefix,".tiff"),width=183,height=215,units="mm",res=600,compression="lzw");print(fig);dev.off()
  fwrite(e$c_dt,file.path(sd,"Figure2C_plotted_KO_null_ratios.csv"))
  fwrite(rbindlist(list(e$bacteria_panel$tip_order[,domain:="Bacteria"],e$archaea_panel$tip_order[,domain:="Archaea"])),file.path(sd,"Figure2B_tree_tip_order.csv"))
  for(name in c("bacteria_monophyletic_lineages.tree","archaea_monophyletic_lineages.tree","monophyletic_lineage_summary.csv"))
    file.copy(file.path("code/result_raw/ISME_Figure2_monophyletic_lineage_source_20260717_v2",name),file.path(sd,name),overwrite=FALSE)
  fwrite(e$catalogue,file.path(sd,"Figure2A_catalogue_counts.csv"))
  writeLines(c("Tree geometry unchanged; branch-specific marker breadth and null ratios updated.",
               "Strict null scores use primary trees pruned to quality strata; independent strict-tree checks remain separate.",
               "Eight maintenance/energy modules shown, surface-retention module not shown in lineage mosaic.",
               paste0("plotted_KOs_across_domain_quality_groups=",nrow(e$c_dt))),file.path(out,"figure_contract_and_QA.txt"))
}
