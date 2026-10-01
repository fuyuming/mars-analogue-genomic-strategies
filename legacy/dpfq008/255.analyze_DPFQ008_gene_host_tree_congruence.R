#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(vegan)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: script HMM_HITS MAG_MANIFEST GENE_TREE SPECIES_TREE OUT_DIR")
}
hit_file <- args[[1L]]
manifest_file <- args[[2L]]
gene_tree_file <- args[[3L]]
species_tree_file <- args[[4L]]
out_dir <- args[[5L]]
if (dir.exists(out_dir) || file.exists(out_dir)) stop("Refusing to overwrite output directory")
stopifnot(all(file.exists(args[1:4])))

hits <- fread(hit_file)
manifest <- fread(manifest_file)
gene_tree <- read.tree(gene_tree_file)
species_tree <- read.tree(species_tree_file)

x <- merge(
  hits[strict == 1L],
  manifest[, .(
    genome_id, ani_cluster_95, representative_genome_id, gtdb_taxonomy_r226
  )],
  by = "genome_id", all.x = TRUE
)
stopifnot(
  nrow(x) == 63L,
  uniqueN(x$genome_id) == 61L,
  uniqueN(x$ani_cluster_95) == 58L,
  !anyNA(x$representative_genome_id)
)

# One deterministic protein per ANI95 independent unit: strongest domain i-E,
# then greatest target coverage, then stable lexical protein identifier.
setorder(x, representative_genome_id, domain_ievalue, -target_coverage, protein_id)
best <- x[, .SD[1L], by = representative_genome_id]
stopifnot(
  nrow(best) == 58L,
  all(best$protein_id %in% gene_tree$tip.label),
  all(best$representative_genome_id %in% species_tree$tip.label)
)

gene_pruned <- keep.tip(gene_tree, best$protein_id)
tip_map <- setNames(best$representative_genome_id, best$protein_id)
gene_pruned$tip.label <- unname(tip_map[gene_pruned$tip.label])
species_pruned <- keep.tip(species_tree, best$representative_genome_id)

ids <- sort(best$representative_genome_id)
gene_dist <- cophenetic.phylo(gene_pruned)[ids, ids]
species_dist <- cophenetic.phylo(species_pruned)[ids, ids]
upper <- upper.tri(gene_dist)
rho <- cor(gene_dist[upper], species_dist[upper], method = "spearman")
set.seed(20260722)
mantel_fit <- mantel(
  as.dist(gene_dist), as.dist(species_dist),
  method = "spearman", permutations = 9999
)

# Sensitivity analyses retain the ANI95 representative as the independent unit.
# They do not replace the prespecified 58-representative primary analysis.
mantel_subset <- function(id_subset, seed) {
  stopifnot(length(id_subset) >= 4L)
  set.seed(seed)
  fit <- mantel(
    as.dist(gene_dist[id_subset, id_subset]),
    as.dist(species_dist[id_subset, id_subset]),
    method = "spearman", permutations = 9999
  )
  data.table(
    n_representatives = length(id_subset),
    n_distance_pairs = choose(length(id_subset), 2L),
    spearman_rho = unname(fit$statistic),
    permutations = fit$permutations,
    empirical_p = unname(fit$signif)
  )
}

dataset_by_id <- setNames(best$dataset_id, best$representative_genome_id)
qaidam_id <- "qaidam_basin_mmag_1773"
qaidam_ids <- ids[dataset_by_id[ids] == qaidam_id]
non_qaidam_ids <- ids[dataset_by_id[ids] != qaidam_id]
stopifnot(length(qaidam_ids) == 51L, length(non_qaidam_ids) == 7L)

qaidam_fit <- mantel_subset(qaidam_ids, 20260723)
non_qaidam_fit <- mantel_subset(non_qaidam_ids, 20260724)

# Retain every pair containing at least one non-Qaidam representative. This
# preserves all 58 units while removing Qaidam-Qaidam pairs from the statistic.
pair_mask <- upper & (
  outer(ids %in% non_qaidam_ids, ids %in% non_qaidam_ids, FUN = "|")
)
masked_rho <- cor(gene_dist[pair_mask], species_dist[pair_mask], method = "spearman")
set.seed(20260725)
masked_null <- replicate(9999L, {
  perm <- sample.int(length(ids))
  permuted_gene <- gene_dist[perm, perm]
  cor(permuted_gene[pair_mask], species_dist[pair_mask], method = "spearman")
})
masked_p <- (1 + sum(abs(masked_null) >= abs(masked_rho))) / (length(masked_null) + 1)

nearest <- function(x) {
  diag(x) <- Inf
  colnames(x)[max.col(-x, ties.method = "first")]
}
gene_nearest <- nearest(gene_dist)
species_nearest <- nearest(species_dist)

tip_data <- best[, .(
  representative_genome_id, protein_id, genome_id, ani_cluster_95,
  dataset_id, gtdb_taxonomy_r226, length_aa, cxxch_count,
  domain_ievalue, target_coverage
)]
setorder(tip_data, representative_genome_id)
dataset_lookup <- setNames(tip_data$dataset_id, tip_data$representative_genome_id)
nearest_data <- data.table(
  representative_genome_id = ids,
  gene_nearest_representative = gene_nearest,
  species_nearest_representative = species_nearest,
  same_nearest_neighbor = gene_nearest == species_nearest,
  dataset_id = dataset_lookup[ids],
  gene_nearest_dataset = dataset_lookup[gene_nearest]
)
nearest_data[, gene_nearest_same_dataset := dataset_id == gene_nearest_dataset]

summary <- data.table(
  independent_unit = "ANI95 species representative",
  n_representatives = length(ids),
  gene_protein_selection = "lowest domain i-E; then highest target coverage; then lexical protein ID",
  distance_measure = "patristic distance",
  association = "Spearman correlation with Mantel label permutation",
  spearman_rho = rho,
  permutations = 9999L,
  empirical_p = unname(mantel_fit$signif),
  identical_nearest_neighbors = sum(nearest_data$same_nearest_neighbor),
  nearest_neighbor_same_dataset = sum(nearest_data$gene_nearest_same_dataset),
  claim_boundary = paste(
    "Tree-distance association supports lineage-structured sequence evolution;",
    "it does not by itself establish vertical inheritance, horizontal transfer, causality or biochemical function"
  )
)

sensitivity <- rbindlist(list(
  data.table(
    analysis = "Primary: all ANI95 representatives",
    n_representatives = length(ids),
    n_distance_pairs = choose(length(ids), 2L),
    spearman_rho = rho,
    permutations = 9999L,
    empirical_p = unname(mantel_fit$signif),
    interpretation = "Prespecified primary gene-host patristic-distance association"
  ),
  cbind(
    data.table(
      analysis = "Sensitivity: Qaidam ANI95 representatives only",
      interpretation = "Within-dominant-source association"
    ),
    qaidam_fit
  ),
  cbind(
    data.table(
      analysis = "Sensitivity: non-Qaidam ANI95 representatives only",
      interpretation = "Small cross-environment subset; underpowered and descriptive"
    ),
    non_qaidam_fit
  ),
  data.table(
    analysis = "Sensitivity: pairs containing at least one non-Qaidam representative",
    n_representatives = length(ids),
    n_distance_pairs = sum(pair_mask),
    spearman_rho = masked_rho,
    permutations = 9999L,
    empirical_p = masked_p,
    interpretation = "Label permutation after excluding all Qaidam-Qaidam distance pairs"
  )
), use.names = TRUE)

dir.create(out_dir, recursive = TRUE)
fwrite(tip_data, file.path(out_dir, "DPFQ008_tree_tip_metadata.tsv"), sep = "\t")
fwrite(nearest_data, file.path(out_dir, "DPFQ008_nearest_neighbor_audit.tsv"), sep = "\t")
fwrite(summary, file.path(out_dir, "DPFQ008_gene_host_tree_congruence.tsv"), sep = "\t")
fwrite(sensitivity, file.path(out_dir, "DPFQ008_Qaidam_dominance_sensitivity.tsv"), sep = "\t")
write.tree(gene_pruned, file.path(out_dir, "DPFQ008_gene_tree_58_ANI95.treefile"))
write.tree(species_pruned, file.path(out_dir, "DPFQ008_host_tree_58_ANI95.treefile"))
writeLines(c(
  "DPFQ008 gene-host tree congruence",
  "strict_HMM_proteins=63",
  "carrier_MAGs=61",
  "independent_ANI95_representatives=58",
  sprintf("spearman_mantel_r=%.8f", rho),
  sprintf("mantel_permutations=%d", 9999L),
  sprintf("mantel_empirical_p=%.8g", unname(mantel_fit$signif)),
  sprintf("identical_nearest_neighbors=%d/58", sum(nearest_data$same_nearest_neighbor)),
  sprintf("gene_nearest_same_dataset=%d/58", sum(nearest_data$gene_nearest_same_dataset)),
  sprintf(
    "non_Qaidam_only_r=%.8f;non_Qaidam_only_p=%.8g;n=%d",
    non_qaidam_fit$spearman_rho, non_qaidam_fit$empirical_p,
    non_qaidam_fit$n_representatives
  ),
  sprintf(
    "non_Qaidam_pair_mask_r=%.8f;non_Qaidam_pair_mask_p=%.8g;pairs=%d",
    masked_rho, masked_p, sum(pair_mask)
  ),
  paste0("claim_boundary=", summary$claim_boundary),
  "status=PASS"
), file.path(out_dir, "run_summary.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
cat(readLines(file.path(out_dir, "run_summary.txt")), sep = "\n")
