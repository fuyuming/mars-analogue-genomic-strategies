#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(ggtext)
  library(patchwork)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 1:2) {
  stop("Usage: 207.make_ISME_figure2_lineage_mosaic.R OUTPUT_DIR [RECURRENCE_NULL_CSV]")
}
out_dir <- args[[1L]]
recurrence_null_path <- if (length(args) == 2L) args[[2L]] else NA_character_
uppercase_panel_tags <- identical(Sys.getenv("ISME_UPPERCASE_PANEL_TAGS", "0"), "1")
panel_tag <- function(value) if (uppercase_panel_tags) toupper(value) else value
publication_text_gate <- identical(Sys.getenv("ISME_PUBLICATION_TEXT_GATE", "0"), "1")
point_text_size <- function(value, minimum = 5.0) {
  if (publication_text_gate) max(value, minimum) else value
}
if (dir.exists(out_dir)) stop("Refusing to overwrite existing output directory: ", out_dir)
dir.create(out_dir, recursive = TRUE)
if (!is.na(recurrence_null_path) && !file.exists(recurrence_null_path)) {
  stop("Recurrence-null input does not exist: ", recurrence_null_path)
}

audit_dir <- "code/result_raw/ISME_Result2_evidence_strict_tree_update_20260717_v1"
lineage_dir <- "code/result_raw/ISME_Figure2_monophyletic_lineage_source_20260717_v2"
paths <- c(
  catalogue = file.path(audit_dir, "Figure2_source_data_a_catalogue_counts.csv"),
  recurrence = file.path(audit_dir, "Figure2_source_data_c_KO_recurrence.csv"),
  modules = file.path(audit_dir, "Figure2_source_data_d_module_root_concordance.csv"),
  lineage_summary = file.path(lineage_dir, "monophyletic_lineage_summary.csv"),
  lineage_breadth = file.path(lineage_dir, "monophyletic_lineage_module_breadth.csv"),
  bacteria_tree = file.path(lineage_dir, "bacteria_monophyletic_lineages.tree"),
  archaea_tree = file.path(lineage_dir, "archaea_monophyletic_lineages.tree")
)
if (!all(file.exists(paths))) stop("One or more required inputs are missing")

catalogue <- fread(paths[["catalogue"]])
rec <- fread(paths[["recurrence"]])
modules <- fread(paths[["modules"]])
lineage_summary <- fread(paths[["lineage_summary"]])
lineage_breadth <- fread(paths[["lineage_breadth"]])

strategy_cols <- c(
  cell_maintenance = "#7465A6",
  energy_acquisition = "#238878",
  surface_retention = "#8B8B8B"
)
strategy_labels <- c(
  cell_maintenance = "Cell maintenance",
  energy_acquisition = "Energy acquisition",
  surface_retention = "Surface retention"
)
module_order <- c(
  "osmotic_desiccation_salt", "cold_protein_quality", "dna_repair_radiation",
  "oxidative_redox", "dormancy_resuscitation", "trace_gas_energy",
  "light_energy", "sulfur_chemolithotrophy"
)
module_short <- c(
  osmotic_desiccation_salt = "Osmotic",
  cold_protein_quality = "Cold / PQC",
  dna_repair_radiation = "DNA repair",
  oxidative_redox = "Oxidative",
  dormancy_resuscitation = "Dormancy",
  trace_gas_energy = "Trace gas",
  light_energy = "Light",
  sulfur_chemolithotrophy = "Sulfur"
)
module_labels <- c(
  osmotic_desiccation_salt = "Osmotic/desiccation/salt",
  cold_protein_quality = "Cold/protein quality",
  dna_repair_radiation = "DNA repair/radiation",
  oxidative_redox = "Oxidative/redox",
  dormancy_resuscitation = "Dormancy/resuscitation",
  trace_gas_energy = "Trace-gas energy",
  light_energy = "Light energy",
  sulfur_chemolithotrophy = "Sulfur chemolithotrophy",
  biofilm_eps_surface = "Biofilm/EPS/surface"
)

theme_pub <- function(base_size = 6.4) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#222222"),
      plot.title = element_text(size = 7.2, face = "bold", margin = margin(b = 3)),
      plot.subtitle = element_text(size = 5.8, colour = "#555555", margin = margin(b = 4)),
      strip.text = element_text(size = 6.2, face = "bold"),
      legend.title = element_text(size = 5.8),
      legend.text = element_text(size = 5.4),
      plot.margin = margin(3, 3, 3, 3)
    )
}
theme_set(theme_pub())

# a: catalogue accounting, retained as a compact provenance panel.
a_dt <- catalogue[stage %chin% c(
  "All source genomes", "Primary MAGs (>=50% complete, <=10% contamination)",
  "ANI95 species representatives", "ANI95 representatives: Bacteria",
  "ANI95 representatives: Archaea", "ANI95 representatives: Unclassified"
)]
a_dt[, label := c("All genomes", "Primary 50/10", "ANI95 species", "Bacteria", "Archaea", "Unclassified")[
  match(stage, c(
    "All source genomes", "Primary MAGs (>=50% complete, <=10% contamination)",
    "ANI95 species representatives", "ANI95 representatives: Bacteria",
    "ANI95 representatives: Archaea", "ANI95 representatives: Unclassified"
  ))
]]
a_nodes <- data.table(
  label = c("All genomes", "Primary 50/10", "ANI95 species", "Bacteria", "Archaea", "Unclassified"),
  x = c(1, 2.25, 3.55, 4.85, 4.85, 4.85), y = c(2, 2, 2, 2.75, 2, 1.25)
)
a_nodes <- merge(a_nodes, a_dt[, .(label, count)], by = "label", all.x = TRUE, sort = FALSE)
a_edges <- data.table(
  x = c(1.34, 2.59, 3.89, 3.89, 3.89), xend = c(1.91, 3.21, 4.51, 4.51, 4.51),
  y = c(2, 2, 2, 2, 2), yend = c(2, 2, 2.75, 2, 1.25)
)
p_a <- ggplot() +
  geom_segment(data = a_edges, aes(x, y, xend = xend, yend = yend),
               linewidth = 0.32, colour = "#AAAAAA", arrow = arrow(length = unit(1.1, "mm"))) +
  geom_label(data = a_nodes, aes(x, y, label = paste0(label, "\n", format(count, big.mark = ","))),
             size = 1.75, family = "Arial", linewidth = 0.22,
             label.padding = unit(0.75, "mm"),
             fill = c("#F2F2F2", "#ECE8F5", "#E2F0ED", "#E8EEF7", "#F4E8DA", "white")) +
  coord_cartesian(xlim = c(0.55, 5.35), ylim = c(0.8, 3.15), clip = "off") +
  labs(title = paste(panel_tag("a"), " Genome catalogue")) +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 7.2, face = "bold"),
        plot.margin = margin(3, 3, 3, 3))

# b: lineage-collapsed trees aligned to module-breadth bubbles.
make_lineage_panel <- function(domain) {
  domain_value <- domain
  tree_path <- paths[[if (domain == "Bacteria") "bacteria_tree" else "archaea_tree"]]
  tree <- read.tree(tree_path)
  summary_dt <- lineage_summary[domain == domain_value]
  breadth_dt <- lineage_breadth[domain == domain_value & module %chin% module_order]
  if (!setequal(tree$tip.label, summary_dt$lineage_id)) stop(domain, " lineage-tree mismatch")
  if (nrow(breadth_dt) != length(tree$tip.label) * length(module_order)) {
    stop(domain, " lineage-module matrix is incomplete")
  }

  p_tree <- ggtree(tree, branch.length = "branch.length", linewidth = 0.24,
                   colour = "#50565A")
  tip_order <- as.data.table(p_tree$data)[isTip == TRUE, .(lineage_id = label, y)]
  summary_dt <- merge(summary_dt, tip_order, by = "lineage_id", all.x = TRUE)
  if (anyNA(summary_dt$y)) stop(domain, " lineage y-position mapping failed")
  setorder(summary_dt, y)
  summary_dt[, display_name := collapsed_name]
  duplicated_name <- summary_dt[, duplicated(display_name) | duplicated(display_name, fromLast = TRUE)]
  summary_dt[duplicated_name, display_name := paste0(display_name, " [", collapsed_rank, "]")]
  summary_dt[, row_label := ""]
  if (domain == "Bacteria") {
    candidate_ids <- summary_dt[genomes >= 20L][order(-genomes, lineage_id), lineage_id]
    selected_ids <- character()
    selected_y <- numeric()
    for (candidate_id in candidate_ids) {
      candidate_y <- summary_dt[lineage_id == candidate_id, y][1L]
      if (!length(selected_y) || all(abs(candidate_y - selected_y) >= 2)) {
        selected_ids <- c(selected_ids, candidate_id)
        selected_y <- c(selected_y, candidate_y)
      }
    }
    summary_dt[lineage_id %chin% selected_ids,
               row_label := paste0(display_name, "  (", genomes, ")")]
  } else {
    summary_dt[, row_label := paste0(display_name, "  (", genomes, ")")]
  }

  n_lineages <- nrow(summary_dt)
  p_tree <- p_tree +
    geom_tippoint(size = if (domain == "Bacteria") 0.45 else 0.8,
                  colour = if (domain == "Bacteria") "#5C6F82" else "#9A6B32") +
    geom_treescale(x = 0, y = 0.8, fontsize = 1.8, linesize = 0.3) +
    scale_y_continuous(limits = c(0.5, n_lineages + 0.5), expand = c(0, 0)) +
    labs(title = paste0(
      domain, " (", n_lineages, " lineages)"
    )) +
    theme_void(base_family = "Arial") +
    theme(plot.title = element_text(size = 6.5, face = "bold", hjust = 0),
          plot.margin = margin(2, 0, 2, 2))

  bubble <- merge(
    breadth_dt,
    summary_dt[, .(lineage_id, y, row_label)],
    by = "lineage_id", all.x = TRUE
  )
  bubble[, module_display := factor(module_short[module], levels = unname(module_short[module_order]))]
  bubble[, strategy_axis := factor(strategy_axis,
    levels = c("cell_maintenance", "energy_acquisition"),
    labels = c("Maintenance", "Energy"))]
  axis_labels <- c(
    paste0("<span style='color:", strategy_cols[["cell_maintenance"]], "'>",
           module_short[module_order[1:5]], "</span>"),
    paste0("<span style='color:", strategy_cols[["energy_acquisition"]], "'>",
           module_short[module_order[6:8]], "</span>")
  )
  names(axis_labels) <- unname(module_short[module_order])
  y_breaks <- summary_dt$y
  names(y_breaks) <- summary_dt$row_label
  label_size <- point_text_size(if (domain == "Bacteria") 4.5 else 5.0)
  point_max <- if (domain == "Bacteria") 3.5 else 5.2

  p_bubble <- ggplot(bubble, aes(module_display, y)) +
    annotate("rect", xmin = 0.5, xmax = 5.5, ymin = 0.5, ymax = n_lineages + 0.5,
             fill = "#F4F1F9", colour = NA) +
    annotate("rect", xmin = 5.5, xmax = 8.5, ymin = 0.5, ymax = n_lineages + 0.5,
             fill = "#EDF7F4", colour = NA) +
    geom_hline(yintercept = seq(1.5, n_lineages - 0.5, by = 1),
               linewidth = 0.12, colour = "#FFFFFF") +
    geom_vline(xintercept = 5.5, linewidth = 0.35, colour = "#FFFFFF") +
    geom_point(aes(size = mean_marker_breadth, fill = strategy_axis),
               shape = 21, colour = "#34383B", stroke = 0.18, alpha = 0.90) +
    scale_fill_manual(values = c(
      Maintenance = strategy_cols[["cell_maintenance"]],
      Energy = strategy_cols[["energy_acquisition"]]
    ), name = NULL) +
    scale_size_area(max_size = point_max, limits = c(0, 1),
                    breaks = c(0.1, 0.3, 0.5),
                    labels = c("0.1", "0.3", "0.5"),
                    name = "Mean marker breadth") +
    scale_x_discrete(position = "top",
                     labels = if (domain == "Bacteria") axis_labels else rep("", length(axis_labels)),
                     drop = FALSE) +
    scale_y_continuous(
      limits = c(0.5, n_lineages + 0.5), breaks = y_breaks,
      labels = names(y_breaks), expand = c(0, 0)
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_family = "Arial", base_size = 5.6) +
    theme(
      panel.grid = element_blank(), axis.ticks = element_blank(),
      axis.text.x.top = if (domain == "Bacteria") {
        element_markdown(size = 5.1, angle = 35, hjust = 0)
      } else {
        element_blank()
      },
      axis.text.y = element_text(size = label_size, colour = "#333333", lineheight = 0.9),
      legend.position = if (domain == "Bacteria") "bottom" else "none",
      legend.direction = "horizontal", legend.box = "horizontal",
      legend.title = element_text(size = 5.2),
      legend.text = element_text(size = point_text_size(4.9)),
      legend.key.width = unit(3.0, "mm"), legend.key.height = unit(2.4, "mm"),
      plot.margin = margin(2, 2, 2, 0)
    ) +
    guides(fill = guide_legend(order = 1, override.aes = list(size = 3)),
           size = guide_legend(order = 2))

  combined <- p_tree + p_bubble + plot_layout(widths = c(0.37, 0.63))
  list(plot = combined, tip_order = tip_order, source = bubble)
}

bacteria_panel <- make_lineage_panel("Bacteria")
archaea_panel <- make_lineage_panel("Archaea")
p_b <- wrap_elements(full =
  (bacteria_panel$plot / archaea_panel$plot) +
    plot_layout(heights = c(4.3, 1)) +
    plot_annotation(
      title = paste(panel_tag("b"), " Phylogenetic distribution of functional breadth"),
      theme = theme(
        plot.title = element_text(size = 7.5, face = "bold", family = "Arial")
      )
    )
)

# c: quantitative recurrence summary or optional source-aware null comparison.
if (is.na(recurrence_null_path)) {
  c_dt <- rec[variable_state == TRUE & strategy_axis %chin% c("cell_maintenance", "energy_acquisition")]
  c_dt[, recurrence_value := normalized_minimum_transitions]
  c_panel_title <- paste(panel_tag("c"), " Minimum functional transitions")
  c_y_label <- "Minimum transitions / tree tips"
  c_panel_mode <- "observed_normalized_minimum_transitions"
} else {
  c_dt <- fread(recurrence_null_path)[
    null_model == "within_source_prevalence_preserving" & variable_trait == TRUE &
      strategy_axis %chin% c("cell_maintenance", "energy_acquisition")
  ]
  expected_c_rows <- nrow(rec[
    variable_state == TRUE & strategy_axis %chin% c("cell_maintenance", "energy_acquisition")
  ])
  if (nrow(c_dt) != expected_c_rows) {
    stop("Null-panel rows do not match the frozen observed recurrence rows")
  }
  if (any(c_dt$observed_score <= 0 | c_dt$null_mean <= 0)) {
    stop("Invalid positive-transition requirement for log-ratio panel")
  }
  c_dt[, recurrence_value := log2(observed_score / null_mean)]
  c_panel_title <- paste(panel_tag("c"), " Phylogenetic constraint relative to source-aware null")
  c_y_label <- "log2 observed / null transitions"
  c_panel_mode <- "log2_observed_over_within_source_null_mean"
}
c_dt[, strategy_axis := factor(
  strategy_axis,
  levels = c("cell_maintenance", "energy_acquisition"),
  labels = c("Maintenance", "Energy")
)]
c_dt[, stratum_label := factor(
  stratum,
  levels = c("primary_ANI95_representatives", "strict_90_5_ANI95_representatives"),
  labels = c("Primary", "Strict 90/5")
)]
c_dt[, domain := factor(domain, levels = c("Bacteria", "Archaea"))]
p_c <- ggplot(c_dt, aes(strategy_axis, recurrence_value, fill = stratum_label)) +
  geom_hline(yintercept = if (is.na(recurrence_null_path)) NA_real_ else 0,
             linetype = 2, linewidth = 0.3, colour = "#777777") +
  geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.3,
               position = position_dodge(width = 0.7), colour = "#444444") +
  geom_point(aes(colour = strategy_axis, group = stratum_label), size = 0.68, alpha = 0.70,
             position = position_jitterdodge(jitter.width = 0.13, dodge.width = 0.7), stroke = 0) +
  facet_wrap(~domain, nrow = 1) +
  scale_colour_manual(values = unname(strategy_cols), guide = "none") +
  scale_fill_manual(values = c("Primary" = "#FFFFFF", "Strict 90/5" = "#CFCFCF")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(title = c_panel_title,
       x = NULL, y = c_y_label, fill = "Genome set") +
  theme_pub() +
  theme(legend.position = "top", panel.spacing = unit(4, "mm"),
        axis.text.x = element_text(size = 5.2))

# d: root-robust environmental-background associations.
d_long <- rbindlist(list(
  modules[, .(module, strategy_axis, stratum, domain = "Archaea", root_class = archaea_root_class)],
  modules[, .(module, strategy_axis, stratum, domain = "Bacteria", root_class = bacteria_root_class)]
))
d_long[, stratum_label := fifelse(stratum == "primary_ANI95_representatives", "Primary", "Strict 90/5")]
d_long[, column := factor(
  paste(domain, stratum_label, sep = "\n"),
  levels = c("Bacteria\nPrimary", "Archaea\nPrimary", "Bacteria\nStrict 90/5", "Archaea\nStrict 90/5")
)]
full_module_order <- c(
  "osmotic_desiccation_salt", "cold_protein_quality", "dna_repair_radiation",
  "oxidative_redox", "dormancy_resuscitation", "trace_gas_energy",
  "light_energy", "sulfur_chemolithotrophy", "biofilm_eps_surface"
)
d_long[, module := factor(module, levels = rev(full_module_order),
                          labels = rev(unname(module_labels[full_module_order])))]
d_long[, status := fifelse(root_class == "BH_positive_all_roots", "Root-robust positive",
                    fifelse(root_class == "BH_nonpositive_all_roots", "BH nonpositive", "Not eligible"))]
d_long[, tile_fill := fifelse(status == "Root-robust positive", strategy_cols[strategy_axis],
                       fifelse(status == "BH nonpositive", "#E3E3E3", "#FFFFFF"))]
d_long[, symbol := fifelse(status == "Root-robust positive", "●",
                    fifelse(status == "BH nonpositive", "–", "×"))]
p_d <- ggplot(d_long, aes(column, module)) +
  geom_tile(aes(fill = tile_fill), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = symbol), size = 2.7, family = "Arial", colour = "#333333") +
  scale_fill_identity() +
  labs(title = paste(panel_tag("d"), " Source-background associations"), x = NULL, y = NULL) +
  theme_minimal(base_size = 6.2, base_family = "Arial") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 5.4),
        axis.text.y = element_text(size = 5.6),
        plot.title = element_text(size = 7.2, face = "bold"),
        plot.margin = margin(3, 3, 3, 3))

design <- "
AAAACCCC
BBBBBBBB
BBBBBBBB
BBBBBBBB
DDDDDDDD
"
figure <- p_a + p_b + p_c + p_d +
  plot_layout(design = design, heights = c(0.58, 1.28, 1.28, 1.28, 0.88))

prefix <- file.path(out_dir, "Figure2_lineage_functional_mosaic")
width_in <- 183 / 25.4
height_in <- 205 / 25.4
svglite::svglite(paste0(prefix, ".svg"), width = width_in, height = height_in)
print(figure)
dev.off()
grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width_in, height = height_in, family = "Arial")
print(figure)
dev.off()
ragg::agg_png(paste0(prefix, ".png"), width = width_in, height = height_in, units = "in", res = 300)
print(figure)
dev.off()
ragg::agg_tiff(paste0(prefix, ".tiff"), width = width_in, height = height_in, units = "in", res = 600,
               compression = "lzw")
print(figure)
dev.off()

fwrite(rbindlist(list(
  bacteria_panel$source[, display_domain := "Bacteria"],
  archaea_panel$source[, display_domain := "Archaea"]
), fill = TRUE), file.path(out_dir, "Figure2b_source_data_lineage_module_breadth.csv.gz"))
fwrite(rbindlist(list(
  bacteria_panel$tip_order[, display_domain := "Bacteria"],
  archaea_panel$tip_order[, display_domain := "Archaea"]
), fill = TRUE), file.path(out_dir, "Figure2b_source_data_lineage_tree_order.csv"))
fwrite(c_dt, file.path(out_dir, "Figure2c_source_data_recurrence.csv"))

writeLines(c(
  "ISME Figure 2 lineage-functional mosaic redesign",
  "status=PASS",
  "backend=R-only",
  "dimensions_mm=183x205",
  paste0("publication_text_gate=", if (publication_text_gate) "PASS_minimum_5pt" else "not_requested"),
  "bacterial_primary_tips=2845",
  "archaeal_primary_tips=233",
  "bacterial_display_lineages=73",
  "archaeal_display_lineages=6",
  "lineage_module_bubbles=632",
  paste0("recurrence_points=", nrow(c_dt)),
  paste0("recurrence_panel_mode=", c_panel_mode),
  "tree_role=locate lineage-aggregated present-day marker breadth across accepted species-tree relationships",
  "bubble_area=mean proportion of prespecified module markers detected across all genomes in each lineage",
  "strict_tree_status=completed; bacterial strict recurrence and environment classes use the independently inferred 879-tip topology",
  "claim_boundary=lineage distributions and minimum transitions do not identify gain/loss direction, adaptive convergence or HGT"
), file.path(out_dir, "run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
