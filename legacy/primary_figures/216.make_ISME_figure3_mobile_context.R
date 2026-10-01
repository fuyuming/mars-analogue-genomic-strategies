#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

root <- normalizePath(".", mustWork = TRUE)
audit_dir <- file.path(root, "code/result_raw/ISME_Result3_3_mobile_context_audit_20260717_v3")
outdir <- file.path(root, "code/result_raw/ISME_Figure3_mobile_context_20260717_v2")
files <- list(
  models = file.path(audit_dir, "Figure3_source_data_models.csv"),
  incidence = file.path(audit_dir, "Figure3_source_data_incidence.csv"),
  datasets = file.path(audit_dir, "Figure3_source_data_dataset_concentration.csv"),
  events = file.path(audit_dir, "Figure3_source_data_energy_KO_events.csv")
)
if (any(!file.exists(unlist(files)))) stop("Missing audited Figure 3 source data")
if (dir.exists(outdir)) stop("Refusing to overwrite existing output directory: ", outdir)
dir.create(outdir, recursive = TRUE)

models <- fread(files$models)
incidence <- fread(files$incidence)
datasets <- fread(files$datasets)
events <- fread(files$events)

stopifnot(
  nrow(models) == 4L,
  nrow(incidence) == 8L,
  nrow(datasets) == 14L,
  sum(events$mobile_cells) == 20L,
  uniqueN(events$dataset_id[events$mobile_cells > 0]) == 2L,
  all(models$fit_status == "ok"),
  all(models$singular == FALSE)
)

COL_MAINT <- "#6F66A8"
COL_ENERGY <- "#2A8C7B"
COL_DEFAULT <- "#A9B0B7"
COL_TEXT <- "#24262A"
COL_GRID <- "#D9DDE1"

theme_set(
  theme_classic(base_size = 7.2, base_family = "Arial") +
    theme(
      text = element_text(colour = COL_TEXT),
      axis.text = element_text(size = 6.5, colour = COL_TEXT),
      axis.title = element_text(size = 7),
      axis.line = element_line(linewidth = 0.35, colour = COL_TEXT),
      axis.ticks = element_line(linewidth = 0.35, colour = COL_TEXT),
      strip.background = element_rect(fill = "white", colour = COL_TEXT, linewidth = 0.35),
      strip.text = element_text(size = 6.5, face = "bold"),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6.2),
      plot.title = element_text(size = 8.2, face = "bold", hjust = 0),
      plot.margin = margin(4, 5, 4, 5)
    )
)

incidence[, stratum_label := factor(stratum_label, levels = c("Primary 50/10", "Strict 90/5"))]
incidence[, preset_label := factor(preset_label, levels = c("Conservative", "Default"))]
incidence[, strategy_label := factor(strategy_label, levels = c("Maintenance", "Energy"))]

p_a <- ggplot(incidence, aes(preset_label, mobile_cells_per_1000, group = strategy_label, colour = strategy_label)) +
  geom_line(linewidth = 0.55, alpha = 0.75) +
  geom_point(size = 2.1) +
  facet_wrap(~stratum_label, nrow = 1) +
  scale_colour_manual(values = c(Maintenance = COL_MAINT, Energy = COL_ENERGY)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(
    title = "a  Mobile-context incidence",
    x = NULL,
    y = "Positive cells per 1,000",
    colour = NULL
  ) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    axis.text.x = element_text(angle = 22, hjust = 1),
    panel.spacing.x = unit(4, "pt")
  )

models[, display_order := factor(model_label, levels = rev(c(
  "Primary 50/10 · Conservative",
  "Strict 90/5 · Conservative",
  "Primary 50/10 · Default",
  "Strict 90/5 · Default"
)))]
models[, preset_label := factor(preset_label, levels = c("Conservative", "Default"))]
models[, stat_label := sprintf("OR %.2f; LRT P=%.3f", odds_ratio_energy, strategy_lrt_p)]

p_b <- ggplot(models, aes(odds_ratio_energy, display_order, colour = preset_label)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.4, colour = "#737980") +
  geom_errorbarh(aes(xmin = conf_low_95, xmax = conf_high_95), height = 0, linewidth = 0.75) +
  geom_point(size = 2.4) +
  scale_x_log10(
    limits = c(0.45, 18),
    breaks = c(0.5, 1, 2, 4, 8, 16),
    labels = label_number(accuracy = 0.1)
  ) +
  scale_colour_manual(values = c(Conservative = COL_ENERGY, Default = COL_DEFAULT)) +
  labs(
    title = "b  Adjusted energy-to-maintenance contrast",
    x = "Odds ratio for mobile context (95% CI)",
    y = NULL,
    colour = "geNomAd preset"
  ) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(colour = COL_GRID, linewidth = 0.3)
  )

dataset_labels <- c(
  alaska_permafrost_reference = "Alaska permafrost",
  australian_basalt_lava_tubes_bay_2025 = "Australian lava tubes",
  stordalen_mire_2019_hybrid_mags = "Stordalen thaw gradient",
  qaidam_basin_mmag_1773 = "Qaidam hyperarid desert",
  mauna_loa_lava_tube_fishman_2023 = "Mauna Loa lava tube",
  atacama_halite_rainfall_prjna484015 = "Atacama halite",
  atacama_salt_crust_prjna351262 = "Atacama salt crust"
)
dataset_order <- names(dataset_labels)
datasets[, dataset_label := factor(dataset_labels[dataset_id], levels = rev(unname(dataset_labels)))]
datasets[, strategy_label := factor(strategy_label, levels = c("Maintenance", "Energy"))]
datasets[, count_label := fifelse(mobile_cells > 0, as.character(mobile_cells), "")]

p_c <- ggplot(datasets, aes(mobile_cells_per_1000, dataset_label, colour = strategy_label, shape = strategy_label)) +
  geom_vline(xintercept = 0, colour = COL_GRID, linewidth = 0.35) +
  geom_point(position = position_dodge(width = 0.48), size = 2.25, stroke = 0.55) +
  geom_text(
    aes(label = count_label),
    position = position_dodge(width = 0.48),
    hjust = -0.55,
    size = 2.05,
    colour = COL_TEXT
  ) +
  scale_colour_manual(values = c(Maintenance = COL_MAINT, Energy = COL_ENERGY)) +
  scale_shape_manual(values = c(Maintenance = 16, Energy = 17)) +
  scale_x_continuous(limits = c(-0.6, 15.4), breaks = c(0, 5, 10, 15), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "c  Dataset concentration (conservative primary set)",
    x = "Positive cells per 1,000",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(colour = COL_GRID, linewidth = 0.3)
  )

events[, dataset_label := factor(dataset_labels[dataset_id], levels = rev(unname(dataset_labels)))]
ko_labels <- unique(events[, .(KO, gene, module)])
ko_labels[, module_order := fifelse(module == "trace_gas_energy", 1L, 2L)]
setorder(ko_labels, module_order, gene, KO)
ko_levels <- ko_labels[, KO]
ko_display <- setNames(sprintf("%s\n%s", ko_labels$gene, ko_labels$KO), ko_labels$KO)
events[, KO_factor := factor(KO, levels = ko_levels)]

p_d <- ggplot(events, aes(KO_factor, dataset_label, fill = mobile_cells)) +
  geom_tile(colour = "#E2E5E8", linewidth = 0.35) +
  geom_text(aes(label = fifelse(mobile_cells > 0, as.character(mobile_cells), "")), size = 2.15, colour = COL_TEXT) +
  scale_fill_gradientn(
    colours = c("white", "#D4ECE7", "#7AC3B6", COL_ENERGY),
    values = rescale(c(0, 1, 3, 5)),
    breaks = 0:5,
    limits = c(0, 5)
  ) +
  scale_x_discrete(labels = ko_display) +
  labs(
    title = "d  Energy-marker composition",
    x = NULL,
    y = NULL,
    fill = "Positive\ncells"
  ) +
  theme(
    axis.text.x = element_text(size = 5.7, angle = 42, hjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.height = unit(18, "pt")
  )

figure <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(widths = c(0.88, 1.12), heights = c(0.92, 1.08), guides = "keep") &
  theme(plot.title.position = "plot")

base <- file.path(outdir, "Figure3_mobile_context")
width_in <- 183 / 25.4
height_in <- 148 / 25.4

svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in, bg = "white")
print(figure)
dev.off()

grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_in, height = height_in, family = "Arial", bg = "white")
print(figure)
dev.off()

ragg::agg_tiff(paste0(base, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, background = "white", compression = "lzw")
print(figure)
dev.off()

ragg::agg_png(paste0(base, ".png"), width = width_in, height = height_in, units = "in", res = 300, background = "white")
print(figure)
dev.off()

fwrite(incidence, file.path(outdir, "Figure3a_source_data.csv"))
fwrite(models, file.path(outdir, "Figure3b_source_data.csv"))
fwrite(datasets, file.path(outdir, "Figure3c_source_data.csv"))
fwrite(events, file.path(outdir, "Figure3d_source_data.csv"))

hashes <- tools::md5sum(unlist(files))
fwrite(data.table(path = names(hashes), md5 = unname(hashes)), file.path(outdir, "input_md5.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(outdir, "R_sessionInfo.txt"))
writeLines(c(
  "status=PASS",
  "width_mm=183",
  "height_mm=148",
  "backend=R_only",
  sprintf("model_rows=%d", nrow(models)),
  sprintf("dataset_rows=%d", nrow(datasets)),
  sprintf("energy_event_grid_rows=%d", nrow(events)),
  sprintf("primary_conservative_energy_positive=%d", sum(events$mobile_cells)),
  "historical_HGT_claim=not_supported"
), file.path(outdir, "run_summary.txt"))

message("PASS: ", outdir)
