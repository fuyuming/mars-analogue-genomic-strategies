#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggtree)
  library(gggenes)
  library(ggnewscale)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) {
  stop(paste(
    "Usage: script HMM_HITS CONTEXT_STRATA NULL_FREQUENCY NEIGHBORHOODS",
    "GENE_TREE MAG_MANIFEST TREE_CONGRUENCE TREE_SENSITIVITY",
    "CONVENTIONAL_EVIDENCE OUT_DIR"
  ))
}
hit_file <- args[[1L]]
context_file <- args[[2L]]
null_file <- args[[3L]]
neighborhood_file <- args[[4L]]
gene_tree_file <- args[[5L]]
manifest_file <- args[[6L]]
congruence_file <- args[[7L]]
sensitivity_file <- args[[8L]]
conventional_file <- args[[9L]]
out_dir <- args[[10L]]
uppercase_panel_tags <- identical(Sys.getenv("ISME_UPPERCASE_PANEL_TAGS", "0"), "1")
panel_tag <- function(value) if (uppercase_panel_tags) toupper(value) else value
publication_text_gate <- identical(Sys.getenv("ISME_PUBLICATION_TEXT_GATE", "0"), "1")
point_text_size <- function(value, minimum = 5.0) {
  if (publication_text_gate) max(value, minimum) else value
}
geom_text_size <- function(value, minimum = 5.0) {
  if (publication_text_gate) max(value, minimum / ggplot2::.pt) else value
}
stopifnot(all(file.exists(args[1:9])))
if (dir.exists(out_dir) || file.exists(out_dir)) stop("Refusing to overwrite output directory")

hits <- fread(hit_file)
context <- fread(context_file)
null_frequency <- fread(null_file)
neighborhoods <- fread(neighborhood_file)
manifest <- fread(manifest_file)
congruence <- fread(congruence_file)
sensitivity <- fread(sensitivity_file)
conventional <- fread(conventional_file)
gene_tree <- ape::read.tree(gene_tree_file)

strict <- hits[strict == 1L]
stopifnot(
  nrow(strict) == 63L,
  uniqueN(strict$genome_id) == 61L,
  uniqueN(strict$dataset_id) == 4L,
  sum(strict$cxxch_count == 2L) == 54L,
  sum(strict$cxxch_count == 1L) == 7L,
  sum(strict$cxxch_count == 0L) == 2L,
  nrow(strict[length_aa >= 300L]) == 52L,
  all(strict[length_aa >= 300L]$cxxch_count == 2L),
  setequal(strict$protein_id, gene_tree$tip.label),
  nrow(context) == 63L,
  sum(context$target_context_positive) == 54L,
  congruence$n_representatives == 58L,
  abs(congruence$spearman_rho - 0.66649036) < 1e-7,
  congruence$permutations == 9999L,
  sensitivity[
    analysis == "Sensitivity: pairs containing at least one non-Qaidam representative",
    n_distance_pairs
  ] == 378L,
  conventional[dpfunc_id == "DPFQ008"]$integrated_ipr_count == 2L,
  conventional[dpfunc_id == "DPFQ008"]$transmembrane_segment_count == 2L
)

palette <- c(
  ink = "#29323A", branch = "#69737C", light = "#E7EAEC",
  lighter = "#F5F6F7", purple = "#70508F", purple_light = "#C9B7D8",
  teal = "#26877D", teal_light = "#A9D2CC", orange = "#D7863B",
  blue = "#397FA1", magenta = "#A95F8F", qaidam = "#8B969F",
  red = "#B9544A"
)
source_palette <- c(
  "Qaidam (56)" = palette[["qaidam"]],
  "Mauna Loa lava tube (3)" = palette[["orange"]],
  "Australian lava tube (2)" = palette[["magenta"]],
  "Alaskan permafrost (2)" = palette[["blue"]]
)
dataset_labels <- c(
  qaidam_basin_mmag_1773 = "Qaidam (56)",
  mauna_loa_lava_tube_fishman_2023 = "Mauna Loa lava tube (3)",
  australian_basalt_lava_tubes_bay_2025 = "Australian lava tube (2)",
  alaska_permafrost_reference = "Alaskan permafrost (2)"
)

theme_pub <- function(base_size = 6.2) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = palette[["ink"]]),
      axis.ticks = element_line(linewidth = 0.3, colour = palette[["ink"]]),
      axis.text = element_text(
        size = point_text_size(base_size * 0.8),
        colour = palette[["ink"]]
      ),
      axis.title = element_text(colour = palette[["ink"]]),
      plot.title = element_text(size = 7.3, face = "bold", colour = palette[["ink"]], margin = margin(b = 3)),
      plot.subtitle = element_text(size = 5.7, colour = palette[["branch"]], margin = margin(b = 3)),
      plot.tag = element_text(size = 8, face = "bold", family = "Helvetica"),
      panel.grid = element_blank(),
      legend.title = element_text(size = 5.7),
      legend.text = element_text(size = 5.3),
      plot.margin = margin(4, 5, 4, 5)
    )
}

# Panel a1: exact representative architecture. Domain and topology counts are
# shown as independent evidence cards rather than invented coordinate spans.
representative_id <- paste(
  "qaidam_basin_mmag_1773",
  "qaidam_basin_mmag_1773_CHGRbin21",
  "k141_121601_length_11112_cov_23.0792_5",
  sep = "|"
)
rep <- strict[protein_id == representative_id]
stopifnot(nrow(rep) == 1L, rep$length_aa == 359L, rep$cxxch_count == 2L)
motif_positions <- as.integer(strsplit(rep$cxxch_positions, ";", fixed = TRUE)[[1L]])
stopifnot(identical(motif_positions, c(134L, 263L)))

p_arch <- ggplot() +
  annotate("segment", x = 1, xend = 359, y = 0.23, yend = 0.23,
           linewidth = 5.5, lineend = "round", colour = palette[["light"]]) +
  annotate("segment", x = 1, xend = 359, y = 0.23, yend = 0.23,
           linewidth = 0.35, colour = palette[["branch"]]) +
  annotate("point", x = motif_positions, y = 0.23, shape = 23, size = 3.2,
           stroke = 0.45, colour = palette[["purple"]], fill = "white") +
  annotate("text", x = motif_positions, y = 0.42, label = c("CXXCH 134", "CXXCH 263"),
           size = geom_text_size(1.75), family = "Helvetica", colour = palette[["purple"]]) +
  annotate("label", x = 90, y = 0.79, label = "InterPro\ncytochrome c-like ×2",
           size = geom_text_size(1.55), family = "Helvetica", linewidth = 0.22,
           label.padding = grid::unit(0.65, "mm"), colour = palette[["ink"]],
           fill = "white") +
  annotate("label", x = 269, y = 0.79, label = "Phobius / TMHMM\nTM segments ×2",
           size = geom_text_size(1.55), family = "Helvetica", linewidth = 0.22,
           label.padding = grid::unit(0.65, "mm"), colour = palette[["ink"]],
           fill = "white") +
  scale_x_continuous(limits = c(-8, 370), breaks = c(1, 100, 200, 300, 359), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
  labs(
    tag = panel_tag("a"), title = "Conserved double-CXXCH architecture",
    subtitle = "Representative protein · 359 aa", x = "Amino-acid position", y = NULL
  ) +
  theme_pub(5.8) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.title.y = element_blank(), plot.margin = margin(4, 4, 1, 6)
  )

# Panel a2: motif completeness by observed protein length.
motif_summary <- strict[, .N, by = .(
  length_group = fifelse(length_aa >= 300L, "≥300 aa", "<300 aa"),
  motif_group = factor(cxxch_count, levels = c(2, 1, 0), labels = c("Two", "One", "None"))
)]
motif_summary[, total := sum(N), by = length_group]
motif_summary[, fraction := N / total]
motif_summary[, label := as.character(N)]
motif_summary[, length_group := factor(length_group, levels = c("<300 aa", "≥300 aa"))]
p_motif <- ggplot(motif_summary, aes(fraction, length_group, fill = motif_group)) +
  geom_col(width = 0.60, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = 1.85, family = "Helvetica", colour = "white") +
  scale_fill_manual(values = c("Two" = palette[["purple"]], "One" = palette[["purple_light"]], "None" = palette[["branch"]])) +
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = c(0, .5, 1), expand = c(0, 0)) +
  labs(x = "Proteins within length class", y = NULL, fill = "CXXCH motifs") +
  theme_pub(5.7) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal",
    legend.key.size = grid::unit(2.2, "mm"),
    plot.margin = margin(1, 4, 1, 6)
  ) +
  guides(fill = guide_legend(nrow = 1, title.position = "left"))

# Panel a3: one complete, traceable neighborhood illustrating the annotation
# classes that enter the context test.
neighborhood_target <- paste(
  "qaidam_basin_mmag_1773",
  "qaidam_basin_mmag_1773_CJ0912_30_35bin52",
  "k141_1489845_length_25643_cov_24.2576_10",
  sep = "|"
)
nb <- copy(neighborhoods[target_protein_id == neighborhood_target])
stopifnot(nrow(nb) == 11L, unique(nb$target_full_window) == 1L, unique(nb$target_context_positive) == 1L)
nb[, x_start := oriented_gene - 0.34]
nb[, x_end := oriented_gene + 0.34]
nb[oriented_strand == "-", c("x_start", "x_end") := .(x_end, x_start)]
nb[, gene_label := fcase(
  gene_class == "target", "DPFQ008",
  gene_class == "cyto_b_N", "cyt b N",
  gene_class == "cyto_b_C", "cyt b C",
  Preferred_name %chin% c("hemL", "hemB"), Preferred_name,
  Preferred_name == "-" & grepl("Uroporphyrinogen III synthase", Description, fixed = TRUE), "Hem4",
  default = ""
)]
nb[nzchar(gene_label), label_y := fifelse(oriented_gene %% 2L == 0L, 0.62, 1.38)]
gene_palette <- c(
  target = palette[["purple"]], cyto_b_N = "#176B7A", cyto_b_C = "#4A9AAC",
  heme = palette[["orange"]], other = "#8D969E", unresolved = "#D6DADD"
)
p_neighborhood <- ggplot(nb) +
  geom_gene_arrow(
    aes(xmin = x_start, xmax = x_end, y = 1, fill = gene_class),
    arrowhead_height = grid::unit(2.1, "mm"), arrowhead_width = grid::unit(0.9, "mm"),
    arrow_body_height = grid::unit(1.35, "mm"), colour = "white", linewidth = 0.18
  ) +
  geom_segment(data = nb[nzchar(gene_label)], aes(x = oriented_gene, xend = oriented_gene, y = 0.92, yend = label_y),
               linewidth = 0.22, colour = palette[["branch"]]) +
  geom_text(data = nb[nzchar(gene_label)], aes(x = oriented_gene, y = label_y, label = gene_label),
            size = geom_text_size(1.65), family = "Helvetica", colour = palette[["ink"]]) +
  scale_fill_manual(values = gene_palette) +
  scale_x_continuous(breaks = -5:5, limits = c(-5.55, 5.55), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.48, 1.52), expand = c(0, 0)) +
  labs(
    title = "Complete representative neighborhood",
    subtitle = "Qaidam Actinomycetota · target oriented left-to-right",
    x = "Relative gene position", y = NULL
  ) +
  theme_pub(5.7) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    legend.position = "none", plot.margin = margin(2, 4, 4, 6)
  )

p_a <- wrap_plots(
  list(p_arch, p_motif, p_neighborhood), ncol = 1,
  heights = c(1.05, 0.75, 1.15)
)

# Panel b: all 63 strict hits, with source, motif and matched-context tracks.
tip_meta <- merge(
  strict,
  context[, .(target_protein_id, target_context_positive)],
  by.x = "protein_id", by.y = "target_protein_id", all.x = TRUE
)
tip_meta <- merge(
  tip_meta,
  manifest[, .(genome_id, gtdb_taxonomy_r226, ani_cluster_95)],
  by = "genome_id", all.x = TRUE, suffixes = c("", ".manifest")
)
tip_meta[, source := factor(dataset_labels[dataset_id], levels = names(source_palette))]
tip_meta[, context_state := factor(
  fifelse(target_context_positive == 1L, "Detected", "Not detected"),
  levels = c("Detected", "Not detected")
)]
stopifnot(!anyNA(tip_meta$source), !anyNA(tip_meta$target_context_positive))

p_tree_base <- ggtree(gene_tree)
p_tree_base$layers[[1L]]$aes_params$colour <- palette[["branch"]]
p_tree_base$layers[[1L]]$aes_params$linewidth <- 0.27
tip_xy <- as.data.table(p_tree_base$data)[isTip == TRUE, .(protein_id = label, x, y)]
tip_tracks <- merge(tip_xy, tip_meta, by = "protein_id", all.x = TRUE)
stopifnot(nrow(tip_tracks) == 63L, !anyNA(tip_tracks$source))
xmax <- max(p_tree_base$data$x, na.rm = TRUE)
offset <- xmax * 0.08
track_x <- c(source = xmax + offset * 0.65, motif = xmax + offset * 2.20, context = xmax + offset * 4.00)
motif_dx <- offset * 0.22
ymax <- max(tip_tracks$y)
masked_sensitivity <- sensitivity[
  analysis == "Sensitivity: pairs containing at least one non-Qaidam representative"
]
stopifnot(nrow(masked_sensitivity) == 1L)

p_b <- p_tree_base +
  geom_point(
    data = tip_tracks, aes(x = track_x[["source"]], y = y, fill = source),
    inherit.aes = FALSE, shape = 21, size = 1.55, stroke = 0.22, colour = "white"
  ) +
  scale_fill_manual(
    values = source_palette, name = "Source catalogue",
    guide = guide_legend(nrow = 2, byrow = TRUE, title.position = "left", override.aes = list(size = 1.7))
  ) +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = tip_tracks[cxxch_count == 2L],
    aes(x = track_x[["motif"]] - motif_dx, y = y),
    inherit.aes = FALSE, shape = 21, size = 0.85, stroke = 0.18,
    colour = palette[["purple"]], fill = palette[["purple"]]
  ) +
  geom_point(
    data = tip_tracks[cxxch_count == 2L],
    aes(x = track_x[["motif"]] + motif_dx, y = y),
    inherit.aes = FALSE, shape = 21, size = 0.85, stroke = 0.18,
    colour = palette[["purple"]], fill = palette[["purple"]]
  ) +
  geom_point(
    data = tip_tracks[cxxch_count == 1L],
    aes(x = track_x[["motif"]], y = y),
    inherit.aes = FALSE, shape = 21, size = 0.95, stroke = 0.18,
    colour = palette[["purple"]], fill = palette[["purple"]]
  ) +
  geom_point(
    data = tip_tracks[cxxch_count == 0L],
    aes(x = track_x[["motif"]], y = y),
    inherit.aes = FALSE, shape = 21, size = 0.95, stroke = 0.28,
    colour = palette[["purple"]], fill = "white"
  ) +
  geom_point(
    data = tip_tracks, aes(x = track_x[["context"]], y = y, fill = context_state),
    inherit.aes = FALSE, shape = 22, size = 1.35, stroke = 0.28, colour = palette[["teal"]]
  ) +
  scale_fill_manual(
    values = c("Detected" = palette[["teal"]], "Not detected" = "white"),
    name = "Cytochrome/haem context", guide = "none"
  ) +
  annotate("text", x = unname(track_x), y = ymax + 2.0,
           label = c("Source", "CXXCH", "Context"), size = geom_text_size(1.7),
           family = "Helvetica", fontface = "bold", colour = palette[["ink"]]) +
  annotate(
    "text", x = track_x[["context"]] / 2, y = -3.2, hjust = 0.5,
    label = sprintf(
      paste0(
        "Gene–host distances · %d ANI95 representatives\n",
        "all pairs: r = %.3f, P = 1.0 × 10⁻⁴\n",
        "non-Qaidam pairs: r = %.3f, P = 1.0 × 10⁻⁴"
      ),
      congruence$n_representatives, congruence$spearman_rho,
      masked_sensitivity$spearman_rho
    ),
    size = geom_text_size(1.62), family = "Helvetica", colour = palette[["ink"]]
  ) +
  geom_treescale(
    x = 0, y = -1.8, width = 0.5, offset = 0.15,
    label = "substitutions/site", linesize = 0.35, fontsize = geom_text_size(1.55),
    color = palette[["ink"]], family = "Helvetica"
  ) +
  coord_cartesian(xlim = c(0, track_x[["context"]] + offset * 0.85), ylim = c(-4.2, ymax + 3.4), clip = "off") +
  labs(
    tag = panel_tag("b"), title = "A cross-environment but lineage-structured family",
    subtitle = "63 strict profile-HMM hits · IQ-TREE2 LG+F+I+G4",
    x = NULL, y = NULL
  ) +
  theme_tree() +
  theme(
    text = element_text(size = 5.8, family = "Helvetica", colour = palette[["ink"]]),
    plot.title = element_text(size = 7.3, face = "bold", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 5.7, colour = palette[["branch"]], margin = margin(b = 3)),
    plot.tag = element_text(size = 8, face = "bold", family = "Helvetica"),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(),
    legend.position = "top", legend.box = "horizontal",
    legend.text = element_text(size = point_text_size(4.64)),
    legend.key.size = grid::unit(2.2, "mm"),
    legend.spacing.y = grid::unit(0.2, "mm"),
    plot.margin = margin(4, 8, 5, 5)
  )

# Panel c: the MAG-level matched null is the independent-unit primary test.
null_mag <- copy(null_frequency[level == "MAG"])
null_mag[, proportion := permutation_frequency / permutations]
perm <- 100000L
observed_mag <- 53L
stopifnot(sum(null_mag$permutation_frequency) == perm)
p_c <- ggplot(null_mag, aes(positive_count, proportion)) +
  annotate("rect", xmin = 1, xmax = 8, ymin = -Inf, ymax = Inf,
           fill = palette[["teal_light"]], alpha = 0.28) +
  geom_col(width = 0.82, fill = palette[["blue"]], colour = "white", linewidth = 0.18) +
  geom_vline(xintercept = 3.75206, linetype = 2, linewidth = 0.35, colour = palette[["branch"]]) +
  geom_vline(xintercept = observed_mag, linewidth = 0.65, colour = palette[["red"]]) +
  annotate("point", x = observed_mag, y = max(null_mag$proportion) * 0.94,
           shape = 21, size = 2.3, stroke = 0.4, fill = "white", colour = palette[["red"]]) +
  annotate("text", x = observed_mag, y = max(null_mag$proportion) * 0.79,
           label = "53/61", size = 2.0, family = "Helvetica", fontface = "bold",
           colour = palette[["red"]]) +
  annotate("text", x = 4.5, y = max(null_mag$proportion) * 0.94,
           label = "null 95%: 1–8", size = geom_text_size(1.65),
           family = "Helvetica", colour = palette[["ink"]]) +
  scale_x_continuous(limits = c(-0.5, 61.5), breaks = c(0, 10, 20, 30, 40, 50, 60), expand = c(0, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.07))) +
  labs(
    tag = panel_tag("c"), title = "Non-random genomic context",
    subtitle = "100,000 within-MAG matched permutations · empirical P = 9.9999 × 10⁻⁶",
    x = "Carrier MAGs with cytochrome/haem context", y = "Permutations"
  ) +
  theme_pub(5.8) +
  theme(plot.margin = margin(4, 4, 4, 6))

left <- wrap_plots(list(p_a, p_c), ncol = 1, heights = c(1.85, 1.0))
figure <- wrap_plots(list(left, p_b), ncol = 2, widths = c(1.02, 1.30))

dir.create(out_dir, recursive = TRUE)
source_dir <- file.path(out_dir, "Source_Data")
dir.create(source_dir)
fwrite(strict[, .(protein_id, dataset_id, genome_id, length_aa, cxxch_count,
                  cxxch_positions, domain_ievalue, target_coverage)],
       file.path(source_dir, "Figure3a_strict_hit_architecture.csv"))
fwrite(nb, file.path(source_dir, "Figure3a_representative_neighborhood.csv"))
fwrite(tip_tracks[, .(protein_id, dataset_id, source, genome_id, ani_cluster_95,
                      length_aa, cxxch_count, target_context_positive,
                      gtdb_taxonomy_r226, x, y)],
       file.path(source_dir, "Figure3b_gene_tree_tip_tracks.csv"))
fwrite(null_mag, file.path(source_dir, "Figure3c_matched_null_frequency.csv"))
fwrite(congruence, file.path(source_dir, "Figure3_statistics.csv"))
fwrite(sensitivity, file.path(source_dir, "Figure3_tree_sensitivity.csv"))

width_mm <- 183
height_mm <- 126
base <- file.path(out_dir, "Figure3_DPFQ008_lineage_context")
svglite::svglite(paste0(base, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4)
print(figure)
dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4, family = "Helvetica")
print(figure)
dev.off()
ragg::agg_tiff(paste0(base, ".tiff"), width = width_mm / 25.4, height = height_mm / 25.4,
               units = "in", res = 600, compression = "lzw")
print(figure)
dev.off()
ragg::agg_png(paste0(base, ".png"), width = width_mm / 25.4, height = height_mm / 25.4,
              units = "in", res = 300)
print(figure)
dev.off()

legend <- paste0(
  "**Fig. 3 | AI-assisted prioritization identifies a lineage-structured double-CXXCH cytochrome-c-like family.** ",
  paste0("**", panel_tag("a"), ",** The 359-aa DPFQ008 representative carries two CXXCH motifs, two integrated InterPro cytochrome-c-like domains and two transmembrane segments predicted independently by Phobius/TMHMM. "),
  "All 52 strict profile-HMM hits of at least 300 aa retained two CXXCH motifs; shorter proteins accounted for all reduced-motif calls. A complete target-oriented Qaidam neighborhood illustrates the cytochrome-b N/C and haem-biosynthesis annotations used in the context definition. ",
  paste0("**", panel_tag("b"), ",** Maximum-likelihood gene tree for all 63 strict profile-HMM hits, coloured by source catalogue and aligned to CXXCH count (filled dots; an open circle denotes no motif) and cytochrome/haem-context detection. The 63 proteins came from 61 MAGs because two MAGs each encoded two hits; the 61 carrier MAGs mapped to 58 ANI95 clusters because three clusters each contained two carrier MAGs. Hits occurred in Qaidam, Mauna Loa and Australian lava tubes, and Alaskan permafrost. After selecting one strongest hit per ANI95 species representative, gene-tree and host-species-tree patristic distances were correlated (Mantel Spearman r = 0.666, 9,999 label permutations, P = 0.0001, n = 58 ANI95 representatives). The association remained after excluding all Qaidam–Qaidam distance pairs (r = 0.305, 378 retained pairs, 9,999 label permutations, P = 0.0001), supporting lineage-structured sequence evolution without establishing exclusively vertical inheritance or horizontal transfer. "),
  paste0("**", panel_tag("c"), ",** Null distribution from 100,000 within-MAG permutations. Each of the 63 targets was matched to proteins from the same MAG by length (±20%) and complete versus contig-truncated ±5-gene window; controls overlapping target windows were excluded. Cytochrome/haem context occurred in 53 of 61 carrier MAGs versus a null mean of 3.75 (95% interval 1–8; one-sided empirical P = 9.9999 × 10⁻⁶). The protein-level sensitivity result was 54/63 with the same empirical P. These data support a conserved, non-random electron-transfer-associated genomic context but do not establish substrate, donor–acceptor pairing, expression, flux, activity or a new pathway.")
)
writeLines(legend, file.path(out_dir, "Figure3_legend.md"))

contract <- c(
  "Core conclusion: AI-assisted triage recovered a cross-environment double-CXXCH cytochrome-c-like family whose sequence evolution follows host lineage structure and whose cytochrome/haem genomic context is non-random.",
  "Figure archetype: asymmetric mixed-modality figure.",
  "Target journal/output: ISME Journal main figure; SVG/PDF/TIFF/PNG.",
  "Backend: R only.",
  "Final size: 183 x 126 mm.",
  "Panel a: exact architecture, motif-length audit and one complete representative neighborhood.",
  "Panel b: 63-tip gene tree with source, motif and context tracks; gene-host tree congruence statistic.",
  "Panel c: MAG-level within-MAG matched permutation null.",
  "Independent unit: ANI95 representative for tree congruence; carrier MAG for primary context enrichment.",
  "Reviewer risk: custom HMM circularity, fragmentary proteins, small non-Qaidam counts, Mantel interpretation and unmeasured biochemical activity.",
  "Claim boundary: no substrate, flux, expression, activity, HGT or new-pathway claim."
)
writeLines(contract, file.path(out_dir, "figure_contract.md"))
writeLines(c(
  "Figure 3 QA notes",
  "- Semantic gates: 63 strict hits; 61 MAGs; four catalogues; 54/63 double CXXCH; 52/52 >=300-aa hits double CXXCH.",
  "- Primary context unit: MAG; observed 53/61; 100,000 matched permutations; empirical P=9.9999e-06.",
  "- Tree congruence unit: one deterministic protein per 58 ANI95 representatives; Mantel Spearman r=0.66649036; 9,999 permutations; P=0.0001.",
  "- Qaidam-dominance sensitivity: after excluding all Qaidam-Qaidam pairs, 378 pairs remain; Spearman r=0.30537817; 9,999 label permutations; P=0.0001.",
  "- Denominators: two MAGs carry two strict proteins each; three ANI95 clusters contain two carrier MAGs each.",
  "- Gene-tree tracks are descriptive for all 63 proteins; inferential tree comparison removes within-ANI95 duplicate weighting.",
  "- Representative neighborhood has a complete +/-5-gene window and is selected by a fixed traceable target identifier.",
  "- No generative imagery or raster scientific content is used; all panels are deterministic R vectors.",
  "- Final pixel/font/overlap and editable-text QA must be completed after rendering."
), file.path(out_dir, "QA_notes.md"))
writeLines(c(
  "ISME Figure 3 DPFQ008 lineage-context redesign",
  "strict_HMM_hits=63", "carrier_MAGs=61", "catalogues=4",
  "strict_double_CXXCH=54/63", "length_ge_300_double_CXXCH=52/52",
  "context_positive_MAGs=53/61", "context_empirical_p=9.9999e-06",
  "tree_ANI95_representatives=58", "gene_host_mantel_r=0.66649036",
  "gene_host_mantel_p=0.0001", "non_Qaidam_pair_mask_r=0.30537817",
  "non_Qaidam_pair_mask_p=0.0001", "backend=R", "size_mm=183x126",
  paste0("panel_tags=", if (uppercase_panel_tags) "uppercase" else "lowercase"),
  paste0("publication_text_gate=", if (publication_text_gate) "PASS_minimum_5pt" else "not_requested"),
  "status=PASS"
), file.path(out_dir, "run_summary.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
cat(readLines(file.path(out_dir, "run_summary.txt")), sep = "\n")
