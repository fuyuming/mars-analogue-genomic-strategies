#!/usr/bin/env Rscript
# Contract: corrected bacterial maintenance source partitioning exceeds energy
# partitioning across quality, lineage blocks and source omissions; archaeal
# boundary results remain visible. Grey intervals are null ranges, not CIs.
suppressPackageStartupMessages(library(data.table))
base<-"code/result_raw/KO_panel_repair_20260926/word_revision/delivery"
code<-paste(readLines("code/370.make_portfolio_supplement_figure_20260803.R"),collapse="\n")
code<-sub('input_dir <- "[^"]+"','input_dir <- corrected_input',code)
code<-sub('output_dir <- "[^"]+"','output_dir <- corrected_output',code)
code<-sub('loo_path <- file.path(input_dir, "portfolio_source_partition_LOO.tsv")','loo_path <- summary_path',code,fixed=TRUE)
code<-sub('formal <- fread(summary_path)','all_rows <- fread(summary_path); formal <- all_rows[is.na(excluded_source) | excluded_source == ""]',code,fixed=TRUE)
code<-sub('loo <- fread(loo_path)','loo <- all_rows[!is.na(excluded_source) & excluded_source != ""]',code,fixed=TRUE)
code<-sub('label = "Δ = −0.0278; P = 0.0001"','label = sprintf("Delta = %.4f\\nP = %.4f", primary$energy_minus_maintenance, primary$difference_p_two_sided)',code,fixed=TRUE)
code<-sub('"Prespecified sensitivity layers"','"Corrected sensitivity layers"',code,fixed=TRUE)
code<-sub('limits = c(-0.115, 0.17)','limits = range(c(formal$energy_minus_maintenance,formal$difference_null_low_95,formal$difference_null_high_95,0)) + c(-.01,.01)',code,fixed=TRUE)
code<-sub('limits = c(-0.045, 0.006)','limits = range(c(loo$energy_minus_maintenance,loo$difference_null_low_95,loo$difference_null_high_95,0)) + c(-.005,.005)',code,fixed=TRUE)
for(branch in c("core65","restored74")) {
 e<-new.env(parent=globalenv())
 e$corrected_input<-file.path(base,"Source_Data",branch,"portfolio")
 e$corrected_output<-file.path(base,"Figures_corrected",paste0("FigureS5_",branch))
 eval(parse(text=code),envir=e)
 dest<-file.path(base,"Source_Data",branch,"FigureS5_panels");dir.create(dest,showWarnings=FALSE)
 for(f in c("panel_a_source_data.tsv","panel_b_source_data.tsv","panel_c_source_data.tsv"))
  file.copy(file.path(e$corrected_output,f),file.path(dest,f),overwrite=TRUE)
}
