#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 209.integrate_de_novo_strict_tree_Figure2_sources.R OUTPUT_DIR")
out_dir <- args[[1L]]
if (dir.exists(out_dir)) stop("Refusing to overwrite: ", out_dir)
dir.create(out_dir, recursive = TRUE)

old_dir <- "code/result_raw/ISME_Result2_evidence_audit_20260717_v1"
cmp_dir <- "code/result_raw/bacterial_tree_primary_strict_comparison_20260717_v3"
catalogue_file <- file.path(old_dir, "Figure2_source_data_a_catalogue_counts.csv")
recurrence_file <- file.path(old_dir, "Figure2_source_data_c_KO_recurrence.csv")
module_file <- file.path(old_dir, "Figure2_source_data_d_module_root_concordance.csv")
strict_rec_file <- file.path(cmp_dir, "strict_recurrence_topology_sensitivity_by_KO.csv")
strict_module_file <- file.path(cmp_dir, "strict_module_environment_class_sensitivity.csv")
inputs <- c(catalogue_file, recurrence_file, module_file, strict_rec_file, strict_module_file)
if (!all(file.exists(inputs))) stop("Missing required input")

catalogue <- fread(catalogue_file)
rec <- fread(recurrence_file)
modules <- fread(module_file)
strict_rec <- fread(strict_rec_file)
strict_modules <- fread(strict_module_file)
strict_name <- "strict_90_5_ANI95_representatives"

target <- rec[domain == "Bacteria" & stratum == strict_name]
if (nrow(target) != 71L || nrow(strict_rec) != 71L || !setequal(target$KO, strict_rec$KO)) {
  stop("Strict bacterial recurrence KO gate failed")
}
strict_map <- strict_rec[, .(
  KO,
  de_novo_score = strict_de_novo_parsimony_score,
  de_novo_norm = strict_de_novo_normalized_transitions
)]
rec[strict_map, on = "KO", `:=`(
  parsimony_score = fifelse(domain == "Bacteria" & stratum == strict_name,
                            i.de_novo_score, parsimony_score),
  normalized_minimum_transitions = fifelse(domain == "Bacteria" & stratum == strict_name,
                                           i.de_novo_norm, normalized_minimum_transitions)
)]
rec[domain == "Bacteria" & stratum == strict_name,
    homoplasy_excess := parsimony_score - minimum_possible_changes]
rec[domain == "Bacteria" & stratum == strict_name,
    consistency_index := fifelse(parsimony_score > 0,
                                 minimum_possible_changes / parsimony_score, NA_real_)]

module_target <- modules[stratum == strict_name]
if (nrow(module_target) != 9L || nrow(strict_modules) != 9L) {
  stop("Strict bacterial module gate failed")
}
modules[strict_modules, on = .(module, strategy_axis),
        bacteria_root_class := fifelse(stratum == strict_name,
                                       i.strict_de_novo_root_class, bacteria_root_class)]

fwrite(catalogue, file.path(out_dir, "Figure2_source_data_a_catalogue_counts.csv"))
fwrite(rec, file.path(out_dir, "Figure2_source_data_c_KO_recurrence.csv"))
fwrite(modules, file.path(out_dir, "Figure2_source_data_d_module_root_concordance.csv"))

summary <- rec[
  domain == "Bacteria" & stratum == strict_name & variable_state == TRUE &
    strategy_axis %chin% c("cell_maintenance", "energy_acquisition"),
  .(variable_KOs = .N,
    median_normalized_minimum_transitions = median(normalized_minimum_transitions)),
  by = strategy_axis
]
fwrite(summary, file.path(out_dir, "strict_de_novo_recurrence_summary.csv"))
writeLines(c(
  "Figure 2 de novo strict-tree source integration",
  "status=PASS",
  "strict_bacterial_tips=879",
  "strict_recurrence_rows_replaced=71",
  "strict_module_rows_replaced=9",
  "strict_source=independently inferred strict bacterial topology",
  "primary_and_archaeal_rows=unchanged"
), file.path(out_dir, "run_summary.txt"))
system2("sha256sum", inputs, stdout = file.path(out_dir, "input_sha256.txt"))
