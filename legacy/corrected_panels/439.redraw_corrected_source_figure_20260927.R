#!/usr/bin/env Rscript
# Contract: lineage-adjusted source associations coexist with source-specific
# marker portfolios; available within-source physical measurements do not
# resolve a universal driver. Primary MAG heatmap is descriptive, not causal.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
base<-"code/result_raw/KO_panel_repair_20260926/word_revision/delivery"
old<-readLines("code/371.make_ISME_figure4_physicochemical_portfolio_20260803.R")
old<-old[seq_len(which(grepl("^# Panel c:",old))-1L)]
code<-paste(old,collapse="\n")
code<-sub('out_dir <- "[^"]+"','out_dir <- corrected_output',code)
inject<-paste0('for (key in names(corrected_paths)) paths[[key]] <- corrected_paths[[key]]; ',
              'stopifnot(all(file.exists(unlist(paths))))')
code<-sub('stopifnot(all(file.exists(unlist(paths))))',inject,code,fixed=TRUE)
code<-sub('all(phys_q$min_q > 0.05)','all(phys_q$min_q >= 0 & phys_q$min_q <= 1)',code,fixed=TRUE)
module_order<-c("osmotic_desiccation_salt","cold_protein_quality","dna_repair_radiation","oxidative_redox",
                "dormancy_resuscitation","trace_gas_energy","light_energy","sulfur_chemolithotrophy","biofilm_eps_surface")
module_labels<-c("Osmotic / ion","Cold / lipid","DNA repair","Oxidative","Sporulation","H2/CO/CF","Light / ion","Sulfur","EPS / surface")
names(module_labels)<-module_order
strategy_cols<-c(cell_maintenance="#7465A6",energy_acquisition="#238878",surface_retention="#8B8B8B")
branches<-commandArgs(trailingOnly=TRUE);if(!length(branches))branches<-c("core65","restored74")
for(branch in branches) {
  module_labels["dormancy_resuscitation"]<-if(branch=="restored74")"Spore / Rpf" else "Sporulation"
  sd<-file.path(base,"Source_Data",branch)
  out<-file.path(base,"Figures_corrected",paste0("Figure3_",branch))
  e<-new.env(parent=globalenv());e$corrected_output<-out
  e$corrected_paths<-list(
    qaidam_group=file.path(sd,"physicochemistry_Qaidam/primary_group_models_10.tsv"),
    qaidam_module=file.path(sd,"physicochemistry_Qaidam/secondary_module_models_45.tsv"),
    tunnel_group=file.path(sd,"physicochemistry_Tunnel/primary_group_zone_models_2.tsv"),
    tunnel_module=file.path(sd,"physicochemistry_Tunnel/secondary_module_zone_models_9.tsv"),
    atacama_tlt=file.path(sd,"physicochemistry_Atacama/TLT_physicochemistry_exact_tests.tsv"),
    atacama_boulder=file.path(sd,"physicochemistry_Atacama/boulder_geochemistry_exact_tests.tsv"),
    atacama_micro=file.path(sd,"physicochemistry_Atacama/boulder_microhabitat_exact_tests.tsv"))
  eval(parse(text=code),envir=e)
  roots<-rbindlist(lapply(c("Bacteria","Archaea"),function(d)fread(file.path(sd,paste0("root_concordance_",d),"marker_module_root_sensitivity_consensus.csv"))))
  roots[,column:=factor(paste(domain,ifelse(stratum=="primary_ANI95_representatives","Primary","Strict"),sep="\n"),
                        levels=c("Bacteria\nPrimary","Archaea\nPrimary","Bacteria\nStrict","Archaea\nStrict"))]
  roots[,module_y:=factor(module,levels=rev(module_order),labels=rev(module_labels))]
  roots[,tile_fill:=ifelse(root_sensitivity_class=="BH_positive_all_roots",strategy_cols[strategy_axis],
                     ifelse(root_sensitivity_class=="BH_nonpositive_all_roots","#E4E4E4","white"))]
  roots[,symbol:=ifelse(root_sensitivity_class=="BH_positive_all_roots","●",ifelse(root_sensitivity_class=="BH_nonpositive_all_roots","–","×"))]
  pa<-ggplot(roots,aes(column,module_y))+geom_tile(aes(fill=tile_fill),colour="white",linewidth=.5)+
    geom_text(aes(label=symbol),size=2.3)+scale_fill_identity()+labs(x=NULL,y=NULL,title="Root-robust source associations")+
    e$theme_pub()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(size=5.5))
  formal<-fread(file.path(sd,"portfolio/portfolio_source_partition_summary.tsv"))[is.na(excluded_source)|excluded_source==""]
  labels<-c("Bacteria: primary","Bacteria: strict 90/5","Bacteria: class blocks","Bacteria: order blocks","Archaea: primary","Archaea: strict")
  stopifnot(nrow(formal)==6)
  formal[,layer:=factor(seq_len(.N),levels=6:1,labels=rev(labels))]
  pb<-ggplot(formal,aes(energy_minus_maintenance,layer))+geom_vline(xintercept=0,linetype=2,colour="grey60",linewidth=.35)+
    geom_errorbar(aes(xmin=difference_null_low_95,xmax=difference_null_high_95),orientation="y",width=0,colour="grey80",linewidth=1)+
    geom_point(aes(colour=domain),size=1.8)+scale_colour_manual(values=c(Bacteria="#315A7D",Archaea="#D07A3A"),guide="none")+
    labs(x=expression(Delta*" source partial "*R^2),y=NULL,title="Energy-minus-maintenance contrast")+e$theme_pub()+
    theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
  heat<-fread(file.path(sd,"Figure3C_primary_MAG_prevalence.tsv"))
  source_ids<-c("qaidam_basin_mmag_1773","atacama_halite_rainfall_prjna484015","atacama_salt_crust_prjna351262","alaska_permafrost_reference",
                "stordalen_mire_2019_hybrid_mags","mauna_loa_lava_tube_fishman_2023","australian_basalt_lava_tubes_bay_2025")
  short<-c("Qaidam","Atacama\nhalite","Atacama\ncrust","Alaska","Stordalen","Mauna\nLoa","Australia")
  counts<-heat[,.(n=first(eligible_MAGs)),by=dataset_id]
  heat[,source:=factor(dataset_id,levels=source_ids,labels=paste0(short,"\n(n=",counts$n[match(source_ids,counts$dataset_id)],")"))]
  order<-unique(heat[order(match(module,module_order),display_order),ko_label])
  heat[,marker:=factor(ko_label,levels=rev(order))]
  pc<-ggplot(heat,aes(source,marker,fill=prevalence))+geom_tile(colour="white",linewidth=.3)+
    scale_fill_gradient(low="#F4F6F5",high="#285F59",limits=c(0,1),breaks=c(0,.5,1),labels=c("0","50","100"),name="MAGs (%)")+
    scale_x_discrete(position="top")+labs(x=NULL,y=NULL,title="Selected markers in primary MAGs")+
    e$theme_pub()+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(size=5),axis.text.y=element_text(size=5.2),
                       legend.position="right",legend.key.height=grid::unit(4,"mm"))
  pd<-e$p_a+labs(title="Physical-measurement coverage")
  pe<-e$p_b+labs(title="Within-source tests")
  design<-"AACCC\nBBCCC\nDDEEE"
  fig<-pa+pb+pc+pd+pe+plot_layout(design=design,heights=c(1.25,1,1.15))+plot_annotation(tag_levels="A")
  stem<-file.path(out,paste0("Figure3_",branch,"_corrected"))
  svglite::svglite(paste0(stem,".svg"),width=183/25.4,height=183/25.4);print(fig);dev.off()
  cairo_pdf(paste0(stem,".pdf"),width=183/25.4,height=183/25.4,family="Arial");print(fig);dev.off()
  ragg::agg_png(paste0(stem,".png"),width=183,height=183,units="mm",res=300);print(fig);dev.off()
  ragg::agg_tiff(paste0(stem,".tiff"),width=183,height=183,units="mm",res=600,compression="lzw");print(fig);dev.off()
  fwrite(roots,file.path(sd,"Figure3A_module_root_concordance.tsv"),sep="\t")
  fwrite(formal,file.path(sd,"Figure3B_portfolio_contrast.tsv"),sep="\t")
  fwrite(e$coverage,file.path(sd,"Figure3D_measurement_coverage.tsv"),sep="\t")
  fwrite(e$phys_q,file.path(sd,"Figure3E_corrected_tests.tsv"),sep="\t")
}
