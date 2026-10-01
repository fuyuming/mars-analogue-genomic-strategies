#!/usr/bin/env Rscript
# Contract: independent WGS supports occurrence of 58/65 core markers;
# A = strategy denominators, B = across-library recurrence, C = gate robustness.
# Quantitative grid, R-only; no inference from four libraries as four studies.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
base <- "code/result_raw/KO_panel_repair_20260926/word_revision/delivery"
out <- file.path(base,"Figures_corrected")
dir.create(out,showWarnings=FALSE)
cols <- c(cell_maintenance="#A56779",energy_acquisition="#497778",surface_retention="#809DA5")
theme_set(theme_classic(base_size=8,base_family="Arial") + theme(
  axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.3),
  plot.title=element_text(size=9,face="bold"),plot.tag=element_text(size=11,face="bold"),
  legend.title=element_blank(),legend.text=element_text(size=7),plot.margin=margin(7,10,7,7)))
for(branch in c("core65","restored74")) {
  inputs <- file.path(base,"Source_Data",branch)
  a <- fread(file.path(inputs,"Figure5A_strategy_detection.csv"))
  b <- fread(file.path(inputs,"Figure5B_detection_frequency.csv"))
  c <- fread(file.path(inputs,"Figure5C_threshold_breadth.csv"))
  stopifnot(sum(a$represented_KOs)==as.integer(sub(".*([0-9]{2})$","\\1",branch)),
            sum(b$KO_count)==sum(a$represented_KOs),
            sum(b[metagenomes_with_dual_mate_detection>0,KO_count])==sum(a$detected_KOs),nrow(c)==12)
  a[, strategy_axis:=factor(strategy_axis,levels=rev(names(cols)))]
  a[, count_label:=paste0(detected_KOs,"/",represented_KOs)]
  pa <- ggplot(a,aes(y=strategy_axis,x=detection_fraction,fill=strategy_axis))+
    geom_col(width=.5)+geom_text(aes(label=count_label),hjust=-.18,size=2.8)+
    scale_fill_manual(values=cols,guide="none")+
    scale_y_discrete(labels=c(cell_maintenance="Maintenance",energy_acquisition="Energy",surface_retention="Surface retention"))+
    scale_x_continuous(limits=c(0,1.3),breaks=c(0,.5,1),expand=expansion(mult=0))+
    labs(x="Fraction detected",y=NULL,title="Detection by strategy")
  pb <- ggplot(b,aes(x=metagenomes_with_dual_mate_detection,y=KO_count))+
    geom_col(fill="#497778",width=.64)+geom_text(aes(label=KO_count),vjust=-.4,size=2.8)+
    scale_x_continuous(breaks=0:4)+scale_y_continuous(expand=expansion(mult=c(0,.13)))+
    labs(x="Metagenomes with detection",y="Number of KOs",title="Across-library occurrence")
  c[, gate:=factor(gate,levels=c("broad","primary","stringent"))]
  labels <- c(SRR18183467="Last Chance Cave",SRR18183468="Indian Tunnel - powder",
              SRR18183469="Indian Tunnel - crystals",SRR18183470="Hidden Cave")
  samplecols <- c(SRR18183467="#6C628A",SRR18183468="#71948B",SRR18183469="#C48B70",SRR18183470="#A56779")
  pc <- ggplot(c,aes(x=gate,y=detected_KOs,group=sample_id,color=sample_id,shape=sample_id))+
    geom_line(linewidth=.55)+geom_point(size=2.2)+
    scale_color_manual(values=samplecols,labels=labels)+scale_shape_manual(values=c(15,16,17,18),labels=labels)+
    scale_x_discrete(labels=c(broad="Broad",primary="Primary",stringent="Stringent"))+
    scale_y_continuous(breaks=seq(48,64,2),limits=c(min(c$detected_KOs)-1,max(c$detected_KOs)+1))+
    labs(x="Alignment gate",y="Detected KOs",title="Alignment-stringency sensitivity")+
    theme(legend.position="right")
  fig <- ((pa|pb)/pc)+plot_layout(heights=c(1,1.1))+plot_annotation(tag_levels="A")
  stem <- file.path(out,paste0("Figure5_",branch,"_corrected"))
  svglite::svglite(paste0(stem,".svg"),width=183/25.4,height=148/25.4);print(fig);dev.off()
  cairo_pdf(paste0(stem,".pdf"),width=183/25.4,height=148/25.4,family="Arial");print(fig);dev.off()
  ragg::agg_png(paste0(stem,".png"),width=183,height=148,units="mm",res=300);print(fig);dev.off()
  ragg::agg_tiff(paste0(stem,".tiff"),width=183,height=148,units="mm",res=600,compression="lzw");print(fig);dev.off()
}
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"))
