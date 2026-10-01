#!/usr/bin/env Rscript
# Post-model export; fail closed on incomplete/nonconverged diagnostics.
suppressPackageStartupMessages(library(data.table))
base<-"code/result_raw/KO_panel_repair_20260926"
delivery<-file.path(base,"word_revision/delivery")
for(branch in c("core65","restored74")) {
 input<-file.path(base,"mobile_repair_20260927",paste0(branch,"_models"))
 models<-fread(file.path(input,"genomad_strategy_context_mixed_models.csv"))
 desc<-fread(file.path(input,"genomad_strategy_context_descriptive.csv"))
 cells<-fread(file.path(input,"genomad_strategy_context_cells.csv.gz"))
 panel<-fread(file.path(base,"amended_inputs",paste0(branch,"_panel.tsv")))
 stopifnot(nrow(models)==4,all(models$fit_status=="ok"),all(models$singular==FALSE))
 out<-file.path(delivery,"Source_Data",branch,"mobile_context");dir.create(out,showWarnings=FALSE)
 for(f in list.files(input,full.names=TRUE))file.copy(f,file.path(out,basename(f)),overwrite=TRUE)
 models[,stratum_label:=fifelse(stratum=="primary_50_10_members","Primary 50/10","Strict 90/5")]
 models[,preset_label:=fifelse(preset=="conservative","Conservative","Default")]
 models[,model_label:=paste(stratum_label,preset_label,sep=" · ")]
 incidence<-desc[strategy_axis %chin% c("cell_maintenance","energy_acquisition")]
 incidence[,strategy_label:=fifelse(strategy_axis=="cell_maintenance","Maintenance","Energy")]
 incidence[,stratum_label:=fifelse(stratum=="primary_50_10_members","Primary 50/10","Strict 90/5")]
 incidence[,preset_label:=fifelse(preset=="conservative","Conservative","Default")]
 incidence[,mobile_cells_per_1000:=1000*mobile_cells/cells]
 primary<-cells[stratum=="primary_50_10_members" & preset=="conservative" & strategy_axis %chin% c("cell_maintenance","energy_acquisition")]
 datasets<-primary[,.(cells=.N,mobile_cells=sum(mobile_context),target_proteins=sum(target_proteins),
                      mobile_target_proteins=sum(mobile_target_proteins),mobile_cells_per_1000=1000*mean(mobile_context)),by=.(dataset_id,strategy_axis)]
 datasets[,strategy_label:=fifelse(strategy_axis=="cell_maintenance","Maintenance","Energy")]
 positive<-primary[strategy_axis=="energy_acquisition" & mobile_context==1,
                    .(mobile_cells=.N,mobile_target_proteins=sum(mobile_target_proteins),clusters=uniqueN(ani_cluster_95)),by=.(dataset_id,KO)]
 events<-merge(CJ(dataset_id=sort(unique(cells$dataset_id)),KO=sort(unique(positive$KO))),positive,by=c("dataset_id","KO"),all.x=TRUE)
 events[is.na(mobile_cells),`:=`(mobile_cells=0L,mobile_target_proteins=0L,clusters=0L)]
 events<-merge(events,panel[,.(KO,gene,module,submodule)],by="KO",all.x=TRUE)
 stopifnot(!anyNA(events$gene),sum(events$mobile_cells)==sum(primary[strategy_axis=="energy_acquisition",mobile_context]))
 # Full symbols remain in panel; plotting uses the first officially listed symbol.
 events[,gene:=trimws(sub(",.*$","",gene))]
 fwrite(models,file.path(out,"Figure3_source_data_models.csv"))
 fwrite(incidence,file.path(out,"Figure3_source_data_incidence.csv"))
 fwrite(datasets,file.path(out,"Figure3_source_data_dataset_concentration.csv"))
 fwrite(events,file.path(out,"Figure3_source_data_energy_KO_events.csv"))
 code<-paste(readLines("code/216.make_ISME_figure3_mobile_context.R"),collapse="\n")
 code<-sub('audit_dir <- file.path(root, "code/result_raw/ISME_Result3_3_mobile_context_audit_20260717_v3")','audit_dir <- corrected_audit',code,fixed=TRUE)
 code<-sub('outdir <- file.path(root, "code/result_raw/ISME_Figure3_mobile_context_20260717_v2")','outdir <- corrected_figure',code,fixed=TRUE)
 code<-sub('limits = c(0.45, 18)','limits = range(c(models$conf_low_95,models$conf_high_95,.5,16))*c(.85,1.15)',code,fixed=TRUE)
 code<-sub('limits = c(-0.6, 15.4)','limits = c(-.6,max(datasets$mobile_cells_per_1000)*1.18)',code,fixed=TRUE)
 e<-new.env(parent=globalenv());e$corrected_audit<-out;e$corrected_figure<-file.path(delivery,"Figures_corrected",paste0("FigureS4_",branch))
 eval(parse(text=code),envir=e)
 fwrite(models,file.path(out,"FigureS4B_models.csv"))
 fwrite(incidence,file.path(out,"FigureS4A_incidence.csv"))
 fwrite(datasets,file.path(out,"FigureS4C_source_concentration.csv"))
 fwrite(events,file.path(out,"FigureS4D_energy_events.csv"))
 writeLines(c("Corrected model fits complete; all four fit-status and singularity gates passed.",
              "Figures and Source Data exported. Word/SI insertion and visual QA remain required.",
              "Present-day mobile context does not demonstrate historical HGT."),file.path(out,"MODEL_EXPORT_COMPLETE.txt"))
}
