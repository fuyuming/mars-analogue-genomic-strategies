#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: script PRIORITY13_AUDIT NEIGHBORHOODS OUT_DIR")
}
priority_file <- args[[1L]]
neighborhood_file <- args[[2L]]
out_dir <- args[[3L]]
if (dir.exists(out_dir) || file.exists(out_dir)) stop("Refusing to overwrite output directory")
if (!file.exists(priority_file) || !file.exists(neighborhood_file)) stop("Input absent")

priority <- fread(priority_file)
neighbors <- fread(neighborhood_file)
stopifnot(
  uniqueN(priority$dpfunc_id) == 13L,
  uniqueN(neighbors[family == "DPFQ008"]$target_protein_id) == 20L,
  uniqueN(neighbors[family == "DPFQ008"]$genome_id) == 19L,
  uniqueN(neighbors[family == "DPFQ008"]$ani_cluster_95) == 19L
)

taxonomy_rank <- function(x, prefix) {
  vapply(strsplit(x, ";", fixed = TRUE), function(parts) {
    hit <- parts[startsWith(parts, prefix)]
    if (!length(hit) || nchar(hit[[1L]]) == nchar(prefix)) "unclassified" else
      substring(hit[[1L]], nchar(prefix) + 1L)
  }, character(1))
}

palette <- c(
  unresolved = "#D7DADF",
  other = "#8B929B",
  heme = "#D69A45",
  cyto_b_C = "#3D91A3",
  cyto_b_N = "#176B7A",
  target = "#75579A",
  dpfq003 = "#A9B2BA",
  dpfq008 = "#75579A",
  gemmatimonadota = "#D69A45",
  actinomycetota = "#59636E"
)

theme_pub <- function(base_size = 6.4) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "#30343A"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#30343A"),
      plot.title = element_text(size = 7.2, face = "bold", margin = margin(b = 3)),
      plot.subtitle = element_text(size = 5.8, colour = "#59636E", margin = margin(b = 3)),
      legend.title = element_text(size = 5.8),
      legend.text = element_text(size = 5.5),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

# Panel a: frozen evidence hierarchy, not an accuracy matrix.
evidence <- rbindlist(list(
  priority[, .(
    dpfunc_id,
    evidence_source = "AI agreement",
    evidence_state = fifelse(
      startsWith(evidence_tier, "A_"), "stronger candidate", "weak candidate"
    )
  )],
  priority[, .(
    dpfunc_id,
    evidence_source = "Domain / topology",
    evidence_state = fcase(
      annotation_status == "supported", "independent support",
      annotation_status == "suggestive", "weak candidate",
      annotation_status == "contradicted", "conflicting",
      default = "none"
    )
  )],
  priority[, .(
    dpfunc_id,
    evidence_source = "DPFunc recovery",
    evidence_state = fifelse(exact_prior_go_recovered == 1L, "weak candidate", "none")
  )],
  priority[, .(
    dpfunc_id,
    evidence_source = "Gene context",
    evidence_state = fcase(
      dpfunc_id == "DPFQ008", "independent support",
      dpfunc_id == "DPFQ003", "weak candidate",
      default = "not evaluated"
    )
  )]
))
evidence[, dpfunc_id := factor(dpfunc_id, levels = rev(sprintf("DPFQ%03d", 1:13)))]
evidence[, evidence_source := factor(
  evidence_source,
  levels = c("AI agreement", "Domain / topology", "DPFunc recovery", "Gene context")
)]
evidence_levels <- c(
  "independent support" = "#318A78",
  "stronger candidate" = "#6FA9BE",
  "weak candidate" = "#B7D2DB",
  "none" = "#E5E7E9",
  "conflicting" = "#C77B62",
  "not evaluated" = "#FFFFFF"
)
p_a <- ggplot(evidence, aes(evidence_source, dpfunc_id, fill = evidence_state)) +
  annotate("rect", xmin = 0.5, xmax = 4.5,
           ymin = which(levels(evidence$dpfunc_id) == "DPFQ008") - 0.48,
           ymax = which(levels(evidence$dpfunc_id) == "DPFQ008") + 0.48,
           fill = "#F2EDF7", colour = NA) +
  geom_tile(width = 0.78, height = 0.76, colour = "white", linewidth = 0.35) +
  scale_fill_manual(values = evidence_levels, drop = FALSE) +
  scale_x_discrete(labels = c("AI\nagreement", "Domain /\ntopology", "DPFunc\nrecovery", "Gene\ncontext")) +
  labs(title = "Orthogonal evidence", x = NULL, y = NULL, fill = NULL) +
  coord_cartesian(clip = "off") +
  theme_pub(5.9) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 5.2, lineheight = 0.9),
    axis.text.y = element_text(size = 5.2),
    legend.position = "bottom", legend.key.size = grid::unit(2.3, "mm"),
    legend.spacing.x = grid::unit(1.1, "mm")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# Panel b: orient every target protein in the positive direction and retain all members.
b <- copy(neighbors[family == "DPFQ008"])
b[, phylum := taxonomy_rank(gtdb_taxonomy_r226, "p__")]
targets <- b[relative_gene == 0L, .(
  target_strand = strand[[1L]], phylum = phylum[[1L]],
  ani_cluster_95 = ani_cluster_95[[1L]], genome_id = genome_id[[1L]]
), by = target_protein_id]
targets[, orient := fifelse(target_strand == "+", 1L, -1L)]
targets[, ani_short := sub("ANI95_", "", ani_cluster_95, fixed = TRUE)]
targets[, copy_index := seq_len(.N), by = .(genome_id, ani_cluster_95)]
targets[, copy_count := .N, by = .(genome_id, ani_cluster_95)]
targets[, row_label := paste0(
  ani_short,
  fifelse(copy_count > 1L, paste0(c("a", "b")[copy_index]), "")
)]
targets[phylum == "Gemmatimonadota", row_label := paste0(row_label, " · Gemm.")]
targets[, context_positive := target_protein_id %in% unique(
  b[relative_gene != 0L & grepl(
    "cytochrome|CXXCH|heme|porphyr",
    paste(Description, Preferred_name, PFAMs), ignore.case = TRUE
  )]$target_protein_id
)]
targets[, non_actinomycetota := phylum != "Actinomycetota"]
setorder(targets, non_actinomycetota, -context_positive, ani_short, row_label)
targets[, row_order := .I]
b <- merge(b, targets[, .(
  target_protein_id, orient, row_label, row_order, phylum, context_positive
)], by = "target_protein_id", suffixes = c("", ".target"))
b[, oriented_gene := relative_gene * orient]
b[, oriented_strand := fifelse(orient == 1L, strand, fifelse(strand == "+", "-", "+"))]
b[, gene_class := fcase(
  relative_gene == 0L, "target",
  grepl("Cytochrom_B_N_2", PFAMs, fixed = TRUE), "cyto_b_N",
  grepl("Cytochrom_B_C", PFAMs, fixed = TRUE), "cyto_b_C",
  grepl("heme|porphyr|tetrapyr|hem[A-Z]", paste(Description, Preferred_name, PFAMs), ignore.case = TRUE), "heme",
  eggnog_annotated == 1L, "other",
  default = "unresolved"
)]
b[, x_start := oriented_gene - 0.32]
b[, x_end := oriented_gene + 0.32]
b[oriented_strand == "-", c("x_start", "x_end") := .(x_end, x_start)]
b[, row_factor := factor(row_label, levels = rev(targets$row_label))]

window <- b[, .(
  min_x = min(oriented_gene), max_x = max(oriented_gene),
  truncated_left = min(oriented_gene) > -5L,
  truncated_right = max(oriented_gene) < 5L,
  phylum = phylum[[1L]], context_positive = context_positive[[1L]],
  row_factor = row_factor[[1L]]
), by = target_protein_id]
stopifnot(
  sum(targets$context_positive) == 17L,
  sum(!window$truncated_left & !window$truncated_right) == 6L,
  uniqueN(targets[context_positive == TRUE]$genome_id) == 16L,
  uniqueN(targets[context_positive == TRUE]$ani_cluster_95) == 16L
)

class_tiles <- unique(b[, .(row_factor, phylum)])
class_tiles[, class_colour := fifelse(
  phylum == "Actinomycetota", "actinomycetota", "gemmatimonadota"
)]
trunc_marks <- rbindlist(list(
  window[truncated_left == TRUE, .(row_factor, x = -5.45, label = "<")],
  window[truncated_right == TRUE, .(row_factor, x = 5.45, label = ">")]
))
gene_labels <- c(
  target = "DPFQ008", cyto_b_N = "Cytochrome b N", cyto_b_C = "Cytochrome b C",
  heme = "Heme/tetrapyrrole", other = "Other annotated", unresolved = "Unresolved"
)
p_b <- ggplot(b, aes(y = row_factor)) +
  geom_tile(
    data = class_tiles,
    aes(x = -5.82, fill = class_colour), width = 0.20, height = 0.72,
    inherit.aes = TRUE, show.legend = FALSE
  ) +
  geom_segment(
    aes(x = x_start, xend = x_end, yend = row_factor, colour = gene_class),
    linewidth = 2.25, lineend = "butt",
    arrow = grid::arrow(length = grid::unit(0.9, "mm"), type = "closed")
  ) +
  geom_text(
    data = trunc_marks, aes(x = x, y = row_factor, label = label),
    inherit.aes = FALSE, colour = "#8B929B", size = 1.8, family = "Helvetica"
  ) +
  scale_colour_manual(values = palette[names(gene_labels)], labels = gene_labels) +
  scale_fill_manual(values = palette[c("actinomycetota", "gemmatimonadota")]) +
  scale_x_continuous(breaks = -5:5, limits = c(-6.05, 5.8), expand = c(0, 0)) +
  labs(
    title = "DPFQ008 genomic neighborhoods",
    x = "Relative gene position", y = "ANI95 cluster", colour = NULL
  ) +
  theme_pub(6.1) +
  theme(
    axis.text.y = element_text(size = 4.9),
    axis.text.x = element_text(size = 5.3),
    legend.position = "bottom", legend.key.width = grid::unit(4.2, "mm"),
    legend.key.height = grid::unit(2.2, "mm"),
    legend.spacing.x = grid::unit(1.0, "mm"),
    plot.margin = margin(4, 4, 4, 2)
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

# Panel c: exact descriptive breadth and lineage/source boundaries.
c_data <- data.table(
  metric = factor(
    c("Protein context", "MAG context", "ANI95 context", "Qaidam carriers", "CADDZG01 carriers"),
    levels = rev(c("Protein context", "MAG context", "ANI95 context", "Qaidam carriers", "CADDZG01 carriers"))
  ),
  family = c("DPFQ003", "DPFQ003", "DPFQ003", "DPFQ008", "DPFQ008"),
  numerator = c(3L, 3L, 2L, 19L, 18L),
  denominator = c(8L, 8L, 5L, 19L, 19L)
)
c_data <- rbind(
  c_data,
  data.table(
    metric = factor(
      c("Protein context", "MAG context", "ANI95 context"),
      levels = levels(c_data$metric)
    ),
    family = "DPFQ008", numerator = c(17L, 16L, 16L), denominator = c(20L, 19L, 19L)
  )
)
c_data[, fraction := numerator / denominator]
c_data[, label := paste0(numerator, "/", denominator)]
c_data[, family := factor(family, levels = c("DPFQ003", "DPFQ008"))]
p_c <- ggplot(c_data, aes(fraction, metric, colour = family, shape = family)) +
  geom_segment(
    aes(x = 0, xend = fraction, yend = metric),
    linewidth = 0.45, colour = "#D7DADF", position = position_dodge(width = 0.42)
  ) +
  geom_point(size = 2.15, stroke = 0.35, position = position_dodge(width = 0.42)) +
  geom_text(
    aes(label = label), hjust = -0.18, size = 1.8, family = "Helvetica",
    position = position_dodge(width = 0.42), show.legend = FALSE
  ) +
  scale_colour_manual(values = setNames(
    unname(palette[c("dpfq003", "dpfq008")]), c("DPFQ003", "DPFQ008")
  )) +
  scale_shape_manual(values = c(DPFQ003 = 1, DPFQ008 = 16)) +
  scale_x_continuous(
    limits = c(0, 1.14), breaks = c(0, 0.5, 1),
    labels = scales::percent_format(accuracy = 1), expand = c(0, 0)
  ) +
  labs(title = "Context breadth and boundaries", x = "Observed fraction", y = NULL, colour = NULL, shape = NULL) +
  theme_pub(5.9) +
  theme(
    axis.text.y = element_text(size = 5.1), axis.text.x = element_text(size = 5.1),
    legend.position = "bottom", legend.key.width = grid::unit(3, "mm"),
    plot.margin = margin(4, 8, 4, 5)
  )

design <- "
BBBAA
BBBCC
"
figure <- p_a + p_b + p_c +
  plot_layout(design = design, guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "Helvetica"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source_dir <- file.path(out_dir, "Source_Data")
dir.create(source_dir, showWarnings = FALSE)
fwrite(evidence, file.path(source_dir, "Figure3a_candidate_evidence_matrix.csv"))
fwrite(b[, .(
  target_protein_id, genome_id, ani_cluster_95, phylum, row_label,
  relative_gene, oriented_gene, strand, oriented_strand, gene_class,
  neighbor_protein_id, Description, Preferred_name, KEGG_ko, PFAMs
)], file.path(source_dir, "Figure3b_DPFQ008_neighborhoods.csv"))
fwrite(c_data, file.path(source_dir, "Figure3c_context_breadth_boundaries.csv"))

base <- file.path(out_dir, "Figure3_AI_neighborhood")
width_in <- 183 / 25.4
height_in <- 165 / 25.4
svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in)
print(figure)
dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_in, height = height_in, family = "Helvetica")
print(figure)
dev.off()
ragg::agg_tiff(paste0(base, ".tiff"), width = width_in, height = height_in, units = "in", res = 600)
print(figure)
dev.off()
ragg::agg_png(paste0(base, ".png"), width = width_in, height = height_in, units = "in", res = 300)
print(figure)
dev.off()

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
summary_lines <- c(
  "ISME Figure 3 AI-neighborhood draft",
  "backend=R-only",
  "figure_width_mm=183",
  "figure_height_mm=165",
  "priority_candidates=13",
  "DPFQ008_family_members=20",
  "DPFQ008_carrier_MAGs=19",
  "DPFQ008_carrier_ANI95=19",
  "DPFQ008_context_positive_members=17",
  "DPFQ008_context_positive_MAGs=16",
  "DPFQ008_context_positive_ANI95=16",
  "DPFQ008_complete_windows=6",
  "claim_gate=broad electron-transfer context; not a new pathway or cross-environment convergence",
  "completed=1",
  "failed=0"
)
writeLines(summary_lines, file.path(out_dir, "run_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
