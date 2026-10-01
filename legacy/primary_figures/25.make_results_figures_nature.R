#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

Sys.setenv(XDG_CACHE_HOME = "/private/tmp")

outdir <- "result_raw/manuscript_update_20260710"
figdir <- file.path(outdir, "figures")
datadir <- file.path(outdir, "plot_data")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)

figure_index <- tibble()

theme_nature <- function(base_size = 7, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(linewidth = 0.28, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.25, colour = "#222222"),
      axis.ticks.length = unit(1.5, "mm"),
      axis.text = element_text(colour = "#222222"),
      axis.title = element_text(colour = "#222222"),
      panel.grid.major.y = element_line(linewidth = 0.18, colour = "#E6E6E6"),
      panel.grid.major.x = element_blank(),
      legend.key.size = unit(3.5, "mm"),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.5),
      plot.title = element_text(face = "bold", size = base_size + 1.2, hjust = 0, colour = "#111111"),
      plot.subtitle = element_text(size = base_size, hjust = 0, colour = "#555555"),
      plot.caption = element_text(size = base_size - 1, colour = "#666666", hjust = 0),
      plot.margin = margin(5.5, 6.5, 5.5, 6.5),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size, colour = "#222222")
    )
}

theme_panel_void <- function(base_size = 7, base_family = "Helvetica") {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(face = "bold", size = base_size + 1.2, hjust = 0, colour = "#111111"),
      plot.subtitle = element_text(size = base_size, hjust = 0, colour = "#555555"),
      plot.caption = element_text(size = base_size - 1, colour = "#666666", hjust = 0),
      plot.margin = margin(5.5, 6.5, 5.5, 6.5)
    )
}

save_plot <- function(plot, name, width, height) {
  pdf_path <- file.path(figdir, paste0(name, ".pdf"))
  png_path <- file.path(figdir, paste0(name, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 600)
  figure_index <<- bind_rows(
    figure_index,
    tibble(
      figure_id = name,
      pdf_path = pdf_path,
      png_path = png_path,
      width_in = width,
      height_in = height
    )
  )
}

dataset_labels <- c(
  alaska_permafrost_reference = "Alaska\npermafrost",
  atacama_halite_rainfall_prjna484015 = "Atacama\nhalite",
  atacama_salt_crust_prjna351262 = "Atacama\nsalt crust",
  qaidam_basin_mmag_1773 = "Qaidam\ncold desert",
  stordalen_mire_2019_hybrid_mags = "Stordalen\nthaw gradient"
)

dataset_cols <- c(
  alaska_permafrost_reference = "#436B9A",
  atacama_halite_rainfall_prjna484015 = "#D99A2B",
  atacama_salt_crust_prjna351262 = "#B95B3F",
  qaidam_basin_mmag_1773 = "#4F9A94",
  stordalen_mire_2019_hybrid_mags = "#5A7F3C"
)

habitat_cols <- c(
  "1" = "#6F8FAF",
  "2" = "#C98F3A",
  "3" = "#8E8E8E",
  "4" = "#9B6A3C",
  "5" = "#4F7F5A"
)

module_cols <- c(
  light_energy = "#D8A22E",
  dna_repair_radiation = "#6D5FA7",
  trace_gas_energy = "#16866F",
  sulfur_chemolithotrophy = "#6C8E2F",
  osmotic_desiccation_salt = "#C76F2A",
  oxidative_redox = "#B2466B",
  biofilm_eps_surface = "#8E6B3D",
  dormancy_resuscitation = "#5F5F5F",
  cold_protein_quality = "#3E75A6"
)

module_labels <- c(
  light_energy = "Light energy",
  dna_repair_radiation = "DNA repair / radiation",
  trace_gas_energy = "Trace-gas energy",
  sulfur_chemolithotrophy = "Sulfur chemolithotrophy",
  osmotic_desiccation_salt = "Osmotic / salt",
  oxidative_redox = "Oxidative redox",
  biofilm_eps_surface = "Biofilm / EPS",
  dormancy_resuscitation = "Dormancy / resuscitation",
  cold_protein_quality = "Cold / protein quality"
)

strategy_cols <- c(
  "Light, repair,\nDNA protection" = "#D8A22E",
  "Trace-gas and\nchemolithotrophy" = "#16866F",
  "Osmotic and\nion homeostasis" = "#C76F2A",
  "Surface colonization\nand EPS" = "#8E6B3D",
  "Oxidative defense,\ndormancy-resuscitation" = "#B2466B"
)

stress_workbook <- "../environmental_stress_methods_formal_table_publish(2).xlsx"
input_manifest <- tibble(
  input_role = c(
    "Updated environmental stress scoring workbook",
    "Five-group environmental cluster assignment",
    "Study-aware 16S PERMANOVA",
    "Leave-one-study-out random forest summary",
    "Driver genera panel",
    "Driver random-forest metrics",
    "MicFunPred stress category table",
    "MAG catalogue Prodigal QC",
    "Targeted KOfam unique protein summary",
    "Targeted KOfam MAG-KO presence summary",
    "Targeted KO prevalence by dataset",
    "mDeepFRI-CNN gene support summary",
    "MAG seed dataset table",
    "GTDB-Tk r226 genome taxonomy",
    "Taxonomy-adjusted focal-KO models"
  ),
  path = c(
    stress_workbook,
    "dataset/env_cluster_k5.tsv",
    "result_raw/study_aware_ai_validation_20260708/study_aware_permanova.csv",
    "result_raw/study_aware_ai_validation_20260708/loso_rf_summary_overall.csv",
    "result_raw/driver_rf_5habitats_g123_boost_20260416/selected_5_drivers_per_group.csv",
    "result_raw/driver_rf_5habitats_g123_boost_20260416/baseline_vs_optimized_metrics.csv",
    "result_raw/micfunpred_kegg_relative_stress_20260424/category_diff_test/stress_category_dataset_level.tsv",
    "result_raw/mars_analog_mag_catalogue_manifest_20260709/prodigal_dataset_qc.csv",
    "result_raw/targeted_kofam_20260709/targeted_kofam_unique_summary_by_dataset.csv",
    "result_raw/targeted_kofam_20260709/targeted_kofam_mag_ko_presence_summary_by_dataset.csv",
    "result_raw/targeted_kofam_module_summary_20260709/ko_prevalence_by_dataset.csv",
    "result_raw/mdeepfri_candidates_20260709/mdeepfri_gene_support_summary.csv",
    "dataset/mars_analog_mag_seed_datasets.csv",
    "result_raw/gtdbtk_r226_taxonomy_control_20260710/gtdb_r226_taxonomy_by_genome.csv",
    "result_raw/gtdbtk_r226_taxonomy_control_models_20260710/key_ko_dataset_effect_after_gtdb_class_control.csv"
  )
)
write_csv(input_manifest, file.path(datadir, "figure_input_manifest.csv"))

env <- read_excel(stress_workbook, sheet = "Group_Index") %>%
  rename(dataset = Group) %>%
  rename(
    R = Radiation,
    D = Aridity,
    T = Temperature,
    S = Salinity,
    `ΔT` = `Thermal amplitude`
  )
stress_evidence <- read_excel(stress_workbook, sheet = "Evidence_Long_Revised")
stress_methods_notes <- read_excel(stress_workbook, sheet = "Original_Method_Notes")
env_cluster <- read_tsv("dataset/env_cluster_k5.tsv", show_col_types = FALSE)
permanova <- read_csv("result_raw/study_aware_ai_validation_20260708/study_aware_permanova.csv", show_col_types = FALSE)
loso <- read_csv("result_raw/study_aware_ai_validation_20260708/loso_rf_summary_overall.csv", show_col_types = FALSE)
drivers <- read_csv("result_raw/driver_rf_5habitats_g123_boost_20260416/selected_5_drivers_per_group.csv", show_col_types = FALSE)
driver_metrics <- read_csv("result_raw/driver_rf_5habitats_g123_boost_20260416/baseline_vs_optimized_metrics.csv", show_col_types = FALSE)
micfun <- read_tsv("result_raw/micfunpred_kegg_relative_stress_20260424/category_diff_test/stress_category_dataset_level.tsv", show_col_types = FALSE)
mag_qc <- read_csv("result_raw/mars_analog_mag_catalogue_manifest_20260709/prodigal_dataset_qc.csv", show_col_types = FALSE)
unique_summary <- read_csv("result_raw/targeted_kofam_20260709/targeted_kofam_unique_summary_by_dataset.csv", show_col_types = FALSE)
mag_ko_presence <- read_csv("result_raw/targeted_kofam_20260709/targeted_kofam_mag_ko_presence_summary_by_dataset.csv", show_col_types = FALSE)
ko_prev <- read_csv("result_raw/targeted_kofam_module_summary_20260709/ko_prevalence_by_dataset.csv", show_col_types = FALSE)
mdeep <- read_csv("result_raw/mdeepfri_candidates_20260709/mdeepfri_gene_support_summary.csv", show_col_types = FALSE)
seed <- read_csv("dataset/mars_analog_mag_seed_datasets.csv", show_col_types = FALSE)
gtdb_taxonomy <- read_csv("result_raw/gtdbtk_r226_taxonomy_control_20260710/gtdb_r226_taxonomy_by_genome.csv", show_col_types = FALSE)
taxonomy_models <- read_csv("result_raw/gtdbtk_r226_taxonomy_control_models_20260710/key_ko_dataset_effect_after_gtdb_class_control.csv", show_col_types = FALSE)

# -------------------------
# Data summaries
# -------------------------
study_overview <- tibble(
  metric = c(
    "16S dataset units",
    "Environmental stress dimensions",
    "Environmental habitat groups",
    "Public MAG/genome datasets used",
    "MAG/genome FASTA files",
    "Predicted proteins",
    "Targeted KOfam KOs",
    "Raw target-KO HMM hits",
    "Unique target proteins",
    "Unique MAG-KO presence records",
    "mDeepFRI-CNN candidate proteins"
  ),
  value = c(
    nrow(env),
    6,
    n_distinct(env_cluster$Group),
    nrow(mag_qc),
    sum(mag_qc$n_manifest),
    sum(mag_qc$total_proteins),
    71,
    sum(unique_summary$raw_hits),
    sum(unique_summary$unique_target_proteins),
    sum(mag_ko_presence$mag_ko_presence_records),
    248
  )
)
write_csv(study_overview, file.path(datadir, "study_overview_summary.csv"))
write_csv(env, file.path(datadir, "environmental_stress_group_index_publish2.csv"))
write_csv(stress_evidence, file.path(datadir, "environmental_stress_evidence_long_publish2.csv"))
write_csv(stress_methods_notes, file.path(datadir, "environmental_stress_method_notes_publish2.csv"))

mag_plot_data <- mag_qc %>%
  left_join(unique_summary, by = "dataset_id") %>%
  left_join(mag_ko_presence, by = "dataset_id") %>%
  mutate(
    dataset_label = recode(dataset_id, !!!dataset_labels),
    dataset_label = factor(dataset_label, levels = dataset_labels[dataset_id]),
    target_protein_fraction = unique_target_proteins / total_proteins
  )
write_csv(mag_plot_data, file.path(datadir, "mag_catalogue_unique_kofam_summary.csv"))

taxonomy_control_status <- gtdb_taxonomy %>%
  group_by(dataset_id) %>%
  summarise(
    classified_genomes = n_distinct(gtdb_id),
    bacterial_genomes = n_distinct(gtdb_id[domain == "Bacteria"]),
    archaeal_genomes = n_distinct(gtdb_id[domain == "Archaea"]),
    n_phyla = n_distinct(phylum),
    n_classes = n_distinct(class),
    .groups = "drop"
  ) %>%
  mutate(
    taxonomy_status = "GTDB-Tk r226 complete",
    ko_control_status = "Genome-level binomial model: GTDB class + dataset"
  )
write_csv(taxonomy_control_status, file.path(datadir, "taxonomy_control_status_table.csv"))

# -------------------------
# Figure 1: study inclusion / framework summary
# -------------------------
stress_long <- env %>%
  select(dataset, R, D, T, pH, S, `ΔT`) %>%
  pivot_longer(-dataset, names_to = "stress_axis", values_to = "score") %>%
  left_join(env_cluster, by = "dataset")
write_csv(stress_long, file.path(datadir, "environment_stress_scores_long.csv"))

infer_dataset_geography <- function(dataset) {
  case_when(
    str_detect(dataset, "^QB1|^QB2") ~ "Qaidam Basin, China",
    str_detect(dataset, "^ANT2") ~ "McMurdo Dry Valleys, Antarctica",
    str_detect(dataset, "^ASK1") ~ "CRREL permafrost tunnel, Alaska, USA",
    str_detect(dataset, "^ADS1") ~ "Ojos del Salado active layer, Andes",
    str_detect(dataset, "^LVA1") ~ "Craters of the Moon lava tubes, Idaho, USA",
    str_detect(dataset, "^LVA2") ~ "Lanzarote lava tubes, Canary Islands, Spain",
    str_detect(dataset, "^ATA1") ~ "Atacama lithic microhabitats, Chile",
    str_detect(dataset, "^ATA2") ~ "Atacama desert pavement/playa, Chile",
    str_detect(dataset, "^ATA3_YU") ~ "Yungay, Atacama Desert, Chile",
    str_detect(dataset, "^ATA3_ME") ~ "Maria Elena, Atacama Desert, Chile",
    str_detect(dataset, "^ATA3_LB") ~ "Lomas Bayas, Atacama Desert, Chile",
    str_detect(dataset, "^ATA3_AL") ~ "Atacama arid transect, Chile",
    str_detect(dataset, "^ATA3_CS") ~ "Atacama coastward transect, Chile",
    str_detect(dataset, "^ATA3_RS") ~ "Atacama regolith/salt transect, Chile",
    TRUE ~ "Unassigned"
  )
}

site_coordinates <- tribble(
  ~site_region, ~longitude, ~latitude, ~coordinate_basis,
  "Qaidam Basin, China", 94.6, 37.4, "regional centroid inferred from Qaidam Basin study groups",
  "McMurdo Dry Valleys, Antarctica", 162.9, -77.5, "regional centroid inferred from Antarctic Dry Valleys study groups",
  "CRREL permafrost tunnel, Alaska, USA", -147.7, 64.9, "site centroid inferred from CRREL permafrost tunnel location",
  "Ojos del Salado active layer, Andes", -68.55, -27.1, "regional centroid inferred from Ojos del Salado active-layer study",
  "Craters of the Moon lava tubes, Idaho, USA", -113.52, 43.42, "regional centroid inferred from CRMO lava-tube groups",
  "Lanzarote lava tubes, Canary Islands, Spain", -13.45, 29.17, "study-level coordinates reported for PRJNA816077",
  "Atacama lithic microhabitats, Chile", -69.45, -23.5, "regional centroid inferred from Atacama lithic microhabitat study",
  "Atacama desert pavement/playa, Chile", -69.9, -24.1, "regional centroid inferred from Atacama pavement/playa study",
  "Yungay, Atacama Desert, Chile", -69.93, -24.08, "locality centroid inferred from Yungay transect label",
  "Maria Elena, Atacama Desert, Chile", -69.67, -22.35, "locality centroid inferred from Maria Elena transect label",
  "Lomas Bayas, Atacama Desert, Chile", -69.32, -23.45, "locality centroid inferred from Lomas Bayas transect label",
  "Atacama arid transect, Chile", -69.6, -23.2, "regional centroid inferred from Atacama transect code",
  "Atacama coastward transect, Chile", -70.05, -24.4, "regional centroid inferred from Atacama transect code",
  "Atacama regolith/salt transect, Chile", -69.75, -24.6, "regional centroid inferred from Atacama transect code"
)

dataset_geography <- env %>%
  transmute(dataset, site_region = infer_dataset_geography(dataset)) %>%
  left_join(site_coordinates, by = "site_region") %>%
  left_join(env_cluster, by = "dataset") %>%
  mutate(coordinate_precision = "study/site-region centroid; not sample-level GPS")
write_csv(dataset_geography, file.path(datadir, "dataset_geography_inferred.csv"))

map_points <- dataset_geography %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  group_by(site_region, longitude, latitude) %>%
  summarise(
    n_dataset = n(),
    dominant_group = names(sort(table(Group), decreasing = TRUE))[1],
    group_mix = paste(sort(unique(paste0("G", Group))), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(dominant_group = factor(dominant_group, levels = names(habitat_cols)))
write_csv(map_points, file.path(datadir, "sampling_landscape_map_points_inferred.csv"))

metric_cards <- study_overview %>%
  mutate(
    display_value = case_when(
      metric == "Predicted proteins" ~ paste0(round(value / 1e6, 1), "M"),
      metric %in% c("Raw target-KO HMM hits", "Unique target proteins", "Unique MAG-KO presence records") ~ label_number(scale_cut = cut_short_scale())(value),
      TRUE ~ as.character(value)
    ),
    metric = recode(
      metric,
      "16S dataset units" = "16S datasets",
      "Environmental stress dimensions" = "Stress axes",
      "Environmental habitat groups" = "Habitat groups",
      "Public MAG/genome datasets used" = "MAG datasets",
      "MAG/genome FASTA files" = "MAG/genomes",
      "Predicted proteins" = "Predicted proteins",
      "Targeted KOfam KOs" = "Targeted KOs",
      "Raw target-KO HMM hits" = "Raw KOfam hits",
      "Unique target proteins" = "Unique proteins",
      "Unique MAG-KO presence records" = "MAG-KO records",
      "mDeepFRI-CNN candidate proteins" = "mDeepFRI proteins"
    ),
    card = row_number(),
    col = if_else(card <= 6, 1, 2),
    row = if_else(card <= 6, 7 - card, 12 - card)
  )

p1a <- ggplot(metric_cards) +
  geom_rect(aes(xmin = col - 0.47, xmax = col + 0.47, ymin = row - 0.38, ymax = row + 0.38),
            fill = "#F5F5F2", colour = "white", linewidth = 0.8) +
  geom_text(aes(x = col - 0.35, y = row + 0.11, label = display_value),
            hjust = 0, size = 3.0, fontface = "bold", colour = "#111111") +
  geom_text(aes(x = col - 0.35, y = row - 0.16, label = metric),
            hjust = 0, size = 2.0, colour = "#555555") +
  coord_cartesian(xlim = c(0.45, 2.55), ylim = c(0.45, 6.55), clip = "off") +
  labs(title = "Integrated data layers") +
  theme_panel_void(7)

world_map <- map_data("world") %>%
  filter(region != "Antarctica" | lat > -85)

p1map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(long, lat, group = group),
    fill = "#F0F0EC",
    colour = "white",
    linewidth = 0.15
  ) +
  geom_point(
    data = map_points,
    aes(longitude, latitude, size = n_dataset, fill = dominant_group),
    shape = 21,
    colour = "#222222",
    stroke = 0.25,
    alpha = 0.92
  ) +
  scale_fill_manual(values = habitat_cols, name = "Dominant\nhabitat") +
  scale_size_continuous(range = c(2.2, 6.5), breaks = c(1, 3, 6, 10), name = "Dataset\nunits") +
  coord_quickmap(xlim = c(-170, 180), ylim = c(-82, 78), expand = FALSE) +
  labs(title = "Sampling landscape") +
  theme_panel_void(7) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6),
    plot.margin = margin(5.5, 1, 5.5, 1)
  )

p1b <- stress_long %>%
  group_by(Group, stress_axis) %>%
  summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(Group = paste0("G", Group), stress_axis = factor(stress_axis, levels = c("R", "D", "T", "pH", "S", "ΔT"))) %>%
  ggplot(aes(stress_axis, Group, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "#F7F7F7", high = "#333333", name = "Mean\nscore") +
  labs(x = "Stress axis", y = "Habitat group", title = "Six-dimensional stress coverage") +
  theme_nature() +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

fig1 <- p1a + p1map + p1b + plot_layout(widths = c(0.95, 1.25, 1))
save_plot(fig1, "Fig1_study_inclusion_stress_framework", 9.2, 3.6)

# -------------------------
# Figure 2: study-aware 16S evidence
# -------------------------
perm_plot <- permanova %>%
  filter(term %in% c("habitat", "study")) %>%
  mutate(
    contrast = case_when(
      model == "habitat_only" ~ "Habitat only",
      model == "study_plus_habitat" & term == "study" ~ "Study term",
      model == "study_plus_habitat" & term == "habitat" ~ "Habitat after study",
      model == "habitat_stratified_permutation_by_study" ~ "Study-stratified habitat",
      TRUE ~ model
    ),
    contrast = factor(contrast, levels = c("Habitat only", "Study term", "Habitat after study", "Study-stratified habitat"))
  )
write_csv(perm_plot, file.path(datadir, "study_aware_permanova_plot_data.csv"))

p2a <- ggplot(perm_plot, aes(contrast, R2, fill = contrast)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = paste0("P=", p_label)), vjust = -0.4, size = 2.2) +
  scale_fill_manual(values = c("#999999", "#4C78A8", "#E45756", "#B279A2"), guide = "none") +
  labs(x = NULL, y = expression(PERMANOVA~R^2), title = "Study-sensitive 16S separation") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

loso_plot <- tibble(
  metric = c("Accuracy", "Balanced accuracy", "Macro-F1"),
  value = c(loso$accuracy[1], loso$balanced_accuracy[1], loso$macro_f1[1])
)
write_csv(loso_plot, file.path(datadir, "loso_rf_overall_plot_data.csv"))

p2b <- ggplot(loso_plot, aes(metric, value)) +
  geom_col(width = 0.6, fill = "#333333") +
  geom_hline(yintercept = 0.2, linetype = 2, linewidth = 0.25, colour = "grey45") +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(x = NULL, y = "Score", title = "Leave-one-study-out transfer") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

fig2 <- p2a + p2b + plot_layout(widths = c(1.35, 1))
save_plot(fig2, "Fig2_study_aware_16S_validation", 7.2, 3.2)

# -------------------------
# Figure 3: MAG catalogue and unique target protein burden
# -------------------------
p3a <- mag_plot_data %>%
  ggplot(aes(dataset_label, n_manifest, fill = dataset_id)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = dataset_cols, guide = "none") +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "MAG/genome count", title = "Genome-resolved catalogue") +
  theme_nature()

p3b <- mag_plot_data %>%
  ggplot(aes(dataset_label, total_proteins / 1e6, fill = dataset_id)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = dataset_cols, guide = "none") +
  labs(x = NULL, y = "Predicted proteins (million)", title = "Uniform Prodigal prediction") +
  theme_nature()

p3c <- mag_plot_data %>%
  ggplot(aes(dataset_label, target_protein_fraction, fill = dataset_id)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = dataset_cols, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Unique target proteins / all proteins", title = "Target-KO protein burden") +
  theme_nature()

fig3 <- p3a + p3b + p3c + plot_layout(ncol = 3)
save_plot(fig3, "Fig3_MAG_catalogue_KOfam_unique_summary", 7.4, 2.8)

# -------------------------
# Figure 4: variable KO prevalence heatmap
# -------------------------
ko_cols <- names(ko_prev)[!(names(ko_prev) %in% c("ko", "gene", "module", "submodule", "mean_prevalence", "sd_prevalence"))]
top_ko <- ko_prev %>%
  arrange(desc(sd_prevalence)) %>%
  slice_head(n = 30) %>%
  mutate(ko_label = paste0(gene, " (", ko, ")")) %>%
  select(ko_label, module, all_of(ko_cols)) %>%
  pivot_longer(all_of(ko_cols), names_to = "dataset_label_full", values_to = "prevalence") %>%
  mutate(
    dataset_id = names(dataset_labels)[match(dataset_label_full, c("Alaska permafrost", "Atacama halite rainfall", "Atacama salt crust", "Qaidam cold hyperarid desert", "Stordalen thaw gradient"))],
    dataset_short = factor(recode(dataset_id, !!!dataset_labels), levels = dataset_labels),
    ko_label = factor(ko_label, levels = rev(unique(ko_label))),
    module = factor(module, levels = names(module_cols))
  )
write_csv(top_ko, file.path(datadir, "top_variable_ko_prevalence_long.csv"))

ko_order <- top_ko %>%
  distinct(ko_label, module) %>%
  mutate(module = factor(module, levels = names(module_cols))) %>%
  left_join(
    top_ko %>% group_by(ko_label) %>% summarise(max_prevalence = max(prevalence, na.rm = TRUE), .groups = "drop"),
    by = "ko_label"
  ) %>%
  arrange(module, desc(max_prevalence), ko_label) %>%
  pull(ko_label)

top_ko <- top_ko %>%
  mutate(
    ko_label = factor(ko_label, levels = rev(ko_order)),
    dataset_short = factor(dataset_short, levels = dataset_labels)
  )

ko_module_annotation <- top_ko %>%
  distinct(ko_label, module) %>%
  mutate(
    ko_label = factor(ko_label, levels = levels(top_ko$ko_label)),
    module_label = module_labels[as.character(module)]
  )

module_positions <- ko_module_annotation %>%
  mutate(y = as.numeric(ko_label)) %>%
  group_by(module, module_label) %>%
  summarise(
    y = mean(range(y)),
    ymin = min(y) - 0.48,
    ymax = max(y) + 0.48,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(short_label = case_when(
    module == "light_energy" ~ "Light\nenergy",
    module == "dna_repair_radiation" ~ "DNA repair /\nradiation",
    module == "trace_gas_energy" ~ "Trace gas\nenergy",
    module == "sulfur_chemolithotrophy" ~ "Sulfur\nchemolith.",
    module == "osmotic_desiccation_salt" ~ "Osmotic /\nsalt",
    module == "oxidative_redox" ~ "Oxidative\nredox",
    module == "biofilm_eps_surface" ~ "Biofilm /\nEPS",
    module == "dormancy_resuscitation" ~ "Dormancy",
    module == "cold_protein_quality" ~ "Cold /\nprotein",
    TRUE ~ module_label
  ))

dataset_annotation <- tibble(
  dataset_short = factor(dataset_labels, levels = dataset_labels),
  dataset_id = names(dataset_labels),
  y = 1
)

p4_top <- ggplot(dataset_annotation, aes(dataset_short, y, fill = dataset_id)) +
  geom_tile(colour = "white", linewidth = 0.4, height = 0.8) +
  scale_fill_manual(values = dataset_cols, guide = "none") +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_panel_void(7) +
  theme(
    axis.text.x.top = element_text(size = 6.8, face = "bold", lineheight = 0.9, colour = "#222222"),
    plot.margin = margin(0, 6, 1, 0)
  )

p4_left <- ggplot() +
  geom_rect(
    data = module_positions,
    aes(xmin = 0.92, xmax = 1.12, ymin = ymin, ymax = ymax, fill = module),
    colour = "white",
    linewidth = 0.45
  ) +
  geom_text(
    data = module_positions,
    aes(x = 0.86, y = y, label = short_label, colour = module),
    inherit.aes = FALSE,
    hjust = 1,
    size = 1.75,
    fontface = "bold",
    lineheight = 0.9
  ) +
  scale_fill_manual(values = module_cols, guide = "none") +
  scale_colour_manual(values = module_cols, guide = "none") +
  scale_x_continuous(limits = c(0.02, 1.18), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, length(levels(top_ko$ko_label)) + 0.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_panel_void(7) +
  theme(plot.margin = margin(0, 1, 5, 0))

p4_heat <- ggplot(top_ko, aes(dataset_short, ko_label, fill = prevalence)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradientn(
    colours = c("#17121A", "#33234C", "#6D3C75", "#B75568", "#E99070", "#F8E7BA"),
    limits = c(0, 1),
    name = "MAG\nprevalence"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_nature(7) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6.6, colour = "#222222"),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(0, 6, 5, 0)
  )

p4_title <- ggplot() +
  annotate("text", x = 0, y = 1, label = "Target KO prevalence across source catalogues",
           hjust = 0, vjust = 1, size = 3.3, fontface = "bold", colour = "#111111") +
  annotate("text", x = 0, y = 0.52,
           label = "Unadjusted MAG prevalence; taxonomy-controlled effects are shown in Fig. S1",
           hjust = 0, vjust = 1, size = 2.35, colour = "#555555") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_panel_void(7) +
  theme(plot.margin = margin(2, 6, 0, 0))

p4 <- p4_title / (plot_spacer() + p4_top + plot_layout(widths = c(0.18, 1))) /
  (p4_left + p4_heat + plot_layout(widths = c(0.25, 1))) +
  plot_layout(heights = c(0.12, 0.12, 1))

save_plot(p4, "Fig4_variable_KO_prevalence_heatmap", 7.6, 5.8)

# -------------------------
# Figure S1: GTDB r226 taxonomy composition and taxonomy-adjusted dataset effects
# -------------------------
top_phyla <- gtdb_taxonomy %>%
  count(phylum, sort = TRUE) %>%
  slice_head(n = 9) %>%
  pull(phylum)

taxonomy_composition <- gtdb_taxonomy %>%
  mutate(
    phylum_plot = if_else(phylum %in% top_phyla, phylum, "Other phyla"),
    dataset_short = factor(recode(dataset_id, !!!dataset_labels), levels = dataset_labels)
  ) %>%
  count(dataset_id, dataset_short, phylum_plot, name = "n_genomes") %>%
  group_by(dataset_id) %>%
  mutate(fraction = n_genomes / sum(n_genomes)) %>%
  ungroup()
write_csv(taxonomy_composition, file.path(datadir, "gtdb_r226_phylum_composition_plot_data.csv"))

phylum_palette <- c(
  "Actinomycetota" = "#436B9A", "Pseudomonadota" = "#D99A2B",
  "Acidobacteriota" = "#5A7F3C", "Halobacteriota" = "#B95B3F",
  "Bacteroidota" = "#4F9A94", "Bacteroidota_A" = "#73A8A2",
  "Chloroflexota" = "#8B6F9B", "Cyanobacteriota" = "#C6A83F",
  "Verrucomicrobiota" = "#777777", "Other phyla" = "#D9D9D9"
)

p_s1a <- ggplot(taxonomy_composition, aes(dataset_short, fraction, fill = phylum_plot)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = phylum_palette, name = "GTDB r226 phylum") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Genome fraction", title = "a  Catalogue taxonomy differs strongly") +
  theme_nature(7) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  theme(legend.position = "bottom")

taxonomy_effect_plot <- taxonomy_models %>%
  left_join(ko_prev %>% distinct(ko, module), by = "ko") %>%
  mutate(
    significant = fit_status == "ok" & p_adj_bh < 0.05,
    ko_label = paste0(gene, " (", ko, ")"),
    ko_label = factor(ko_label, levels = rev(ko_label[order(pseudo_r2_dataset_given_class, na.last = TRUE)])),
    module = factor(module, levels = names(module_cols))
  )
write_csv(taxonomy_effect_plot, file.path(datadir, "taxonomy_adjusted_ko_dataset_effect_plot_data.csv"))

p_s1b <- ggplot(taxonomy_effect_plot, aes(pseudo_r2_dataset_given_class, ko_label, colour = module)) +
  geom_segment(aes(x = 0, xend = pseudo_r2_dataset_given_class, yend = ko_label), colour = "#D2D2D2", linewidth = 0.3) +
  geom_point(aes(shape = significant), size = 1.8, stroke = 0.35) +
  scale_colour_manual(values = module_cols, labels = module_labels, name = "Adaptive module") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "BH-adjusted P < 0.05", `FALSE` = "Not significant / nonconverged"), name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Additional deviance explained by dataset after GTDB-class control",
    y = NULL,
    title = "b  Residual dataset association varies among focal KOs",
    subtitle = "Filled points: BH-adjusted P < 0.05; open points: not significant or nonconverged"
  ) +
  theme_nature(7) +
  guides(
    colour = guide_legend(nrow = 3, byrow = TRUE),
    shape = "none"
  ) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6.2))

fig_s1 <- p_s1a + p_s1b + plot_layout(widths = c(1, 1.35))
save_plot(fig_s1, "FigS1_GTDB_r226_taxonomy_control", 9.0, 5.8)

# -------------------------
# Figure 5: adaptive strategy marker panels + mDeepFRI support
# -------------------------
strategy_map <- tribble(
  ~strategy, ~gene,
  "Light, repair,\nDNA protection", "bop",
  "Light, repair,\nDNA protection", "hop",
  "Light, repair,\nDNA protection", "phr",
  "Light, repair,\nDNA protection", "dps",
  "Trace-gas and\nchemolithotrophy", "hyaA",
  "Trace-gas and\nchemolithotrophy", "hyaB",
  "Trace-gas and\nchemolithotrophy", "hydA",
  "Trace-gas and\nchemolithotrophy", "hydB",
  "Trace-gas and\nchemolithotrophy", "coxL",
  "Trace-gas and\nchemolithotrophy", "coxM",
  "Trace-gas and\nchemolithotrophy", "coxS",
  "Trace-gas and\nchemolithotrophy", "soxZ",
  "Osmotic and\nion homeostasis", "kdpA",
  "Osmotic and\nion homeostasis", "kdpC",
  "Osmotic and\nion homeostasis", "betA",
  "Osmotic and\nion homeostasis", "betB",
  "Osmotic and\nion homeostasis", "otsA",
  "Osmotic and\nion homeostasis", "otsB",
  "Surface colonization\nand EPS", "algD",
  "Surface colonization\nand EPS", "bcsA",
  "Surface colonization\nand EPS", "pgaC",
  "Oxidative defense,\ndormancy-resuscitation", "katG",
  "Oxidative defense,\ndormancy-resuscitation", "katE",
  "Oxidative defense,\ndormancy-resuscitation", "sodB",
  "Oxidative defense,\ndormancy-resuscitation", "relA_spoT"
)

strategy_prev <- ko_prev %>%
  inner_join(strategy_map, by = "gene") %>%
  select(strategy, gene, ko, all_of(ko_cols)) %>%
  pivot_longer(all_of(ko_cols), names_to = "dataset_label_full", values_to = "prevalence") %>%
  mutate(
    dataset_id = names(dataset_labels)[match(dataset_label_full, c("Alaska permafrost", "Atacama halite rainfall", "Atacama salt crust", "Qaidam cold hyperarid desert", "Stordalen thaw gradient"))],
    dataset_short = factor(recode(dataset_id, !!!dataset_labels), levels = dataset_labels),
    marker = paste0(gene, " (", ko, ")")
  )
write_csv(strategy_prev, file.path(datadir, "adaptive_strategy_marker_prevalence_long.csv"))

strategy_summary <- strategy_prev %>%
  group_by(strategy, dataset_id, dataset_short) %>%
  summarise(mean_prevalence = mean(prevalence, na.rm = TRUE), .groups = "drop")
write_csv(strategy_summary, file.path(datadir, "adaptive_strategy_mean_prevalence.csv"))

strategy_summary <- strategy_summary %>%
  mutate(
    strategy = factor(strategy, levels = names(strategy_cols)),
    dataset_short = factor(dataset_short, levels = rev(dataset_labels))
  )

p5a <- ggplot(strategy_summary, aes(strategy, dataset_short, fill = mean_prevalence)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = percent(mean_prevalence, accuracy = 1)), size = 2.15, colour = "#111111") +
  scale_fill_gradientn(
    colours = c("#F7F7F4", "#D7CEC0", "#A89472", "#6F5738", "#2B2118"),
    limits = c(0, 1),
    name = "Mean marker\nprevalence",
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Habitat-specific weighting of shared adaptive modules",
    subtitle = "Mean MAG prevalence across marker KOs in each module"
  ) +
  theme_nature(7) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", lineheight = 0.9),
    panel.grid = element_blank(),
    legend.position = "right"
  )

mdeep_plot <- mdeep %>%
  mutate(
    label = paste0(gene, " / ", toupper(mode)),
    module = factor(module, levels = names(module_cols))
  ) %>%
  arrange(module, desc(proteins_with_predictions))
write_csv(mdeep_plot, file.path(datadir, "mdeepfri_support_plot_data.csv"))

p5b <- ggplot(mdeep_plot, aes(proteins_with_predictions, reorder(label, proteins_with_predictions), colour = module)) +
  geom_segment(aes(x = 0, xend = proteins_with_predictions, yend = reorder(label, proteins_with_predictions)),
               linewidth = 0.45, alpha = 0.7) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = module_cols, labels = module_labels[names(module_cols)], name = "Module") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Proteins with mDeepFRI-CNN support", y = NULL, title = "Deep-learning annotation support") +
  theme_nature(7) +
  theme(legend.position = "right")

fig5 <- p5a / p5b + plot_layout(heights = c(1.05, 1))
save_plot(fig5, "Fig5_adaptive_strategies_mDeepFRI_support", 7.4, 5.7)

# -------------------------
# Figure 6: Mars biosignature prioritization framework
# -------------------------
target_map <- tribble(
  ~martian_setting, ~analogue_support, ~gene_module_targets, ~molecular_or_preservation_targets, ~evidence_level, ~setting_colour,
  "Salts / evaporites / brines", "Atacama halite;\nAtacama salt crust", "bop/hop; phr; dps;\nect/bet/ots", "retinal pigments;\ncarotenoids; solutes; EPS", "KOfam + MAG;\nmDeepFRI subset", "#D99A2B",
  "Rock interiors / weathering crust", "Qaidam;\nendolithic logic", "DNA repair;\noxidative defense;\nEPS/biofilm", "pigments;\nmineral-bound organics;\nEPS residues", "MAG pattern;\nliterature translation", "#8E6B3D",
  "Cold oligotrophic regolith", "Qaidam; Alaska;\nStordalen", "hya/hyd; cox;\nsox; cbb", "Fe-S/redox cofactors;\nsulfur intermediates;\nlow-biomass organics", "KOfam + MAG;\nmDeepFRI subset", "#436B9A",
  "Ancient aqueous sediments / clays", "Preservation logic", "stress-response fragments;\nbiofilm genes", "lipids; hopanoids;\npigments; EPS;\nmineral-bound organics", "Hypothesis-generating\ntranslation", "#6D5FA7",
  "Freeze-thaw / transient hydration", "Permafrost;\nthaw-gradient analogues", "cold-shock;\nmembrane adaptation;\ntrehalose; betaine; EPS", "unsaturated fatty acids;\ncompatible solutes; EPS", "MAG pattern;\nanalogue inference", "#5A7F3C"
)
write_csv(target_map, file.path(datadir, "mars_biosignature_prioritization_framework.csv"))

target_cards <- target_map %>%
  mutate(row = rev(row_number()))

target_headers <- tibble(
  x = c(1.74, 3.0, 4.32, 5.64),
  label = c("Earth analogue", "Gene/module targets", "Molecular or preserved traces", "Evidence")
)

p6 <- ggplot() +
  geom_rect(
    data = target_cards,
    aes(xmin = 0.04, xmax = 6.28, ymin = row - 0.43, ymax = row + 0.43),
    fill = "#F7F7F4",
    colour = "white",
    linewidth = 0.8
  ) +
  geom_rect(
    data = target_cards,
    aes(xmin = 0.04, xmax = 0.14, ymin = row - 0.43, ymax = row + 0.43, fill = martian_setting),
    colour = NA
  ) +
  geom_text(
    data = target_cards,
    aes(x = 0.22, y = row, label = str_wrap(martian_setting, 20)),
    hjust = 0,
    size = 2.35,
    fontface = "bold",
    lineheight = 0.88,
    colour = "#111111"
  ) +
  geom_text(
    data = target_cards,
    aes(x = 1.74, y = row, label = analogue_support),
    size = 2.25,
    lineheight = 0.9,
    colour = "#222222"
  ) +
  geom_text(
    data = target_cards,
    aes(x = 3.0, y = row, label = gene_module_targets),
    size = 2.25,
    lineheight = 0.9,
    colour = "#222222"
  ) +
  geom_text(
    data = target_cards,
    aes(x = 4.32, y = row, label = molecular_or_preservation_targets),
    size = 2.15,
    lineheight = 0.9,
    colour = "#222222"
  ) +
  geom_text(
    data = target_cards,
    aes(x = 5.64, y = row, label = evidence_level),
    size = 2.05,
    lineheight = 0.9,
    colour = "#444444"
  ) +
  geom_text(
    data = target_headers,
    aes(x = x, y = max(target_cards$row) + 0.62, label = label),
    size = 2.45,
    fontface = "bold",
    colour = "#222222"
  ) +
  scale_fill_manual(values = setNames(target_cards$setting_colour, target_cards$martian_setting), guide = "none") +
  coord_cartesian(xlim = c(0, 6.34), ylim = c(0.42, max(target_cards$row) + 0.78), clip = "off") +
  labs(
    title = "Earth-analogue-informed, habitat-stratified biosignature prioritization",
    subtitle = "Observed MAG/KOfam signals are translated conservatively into location-specific Mars search targets"
  ) +
  theme_panel_void(7) +
  theme(plot.margin = margin(8, 8, 8, 8))
save_plot(p6, "Fig6_Mars_biosignature_prioritization_framework", 8.2, 4.1)

write_csv(figure_index, file.path(datadir, "figure_file_index.csv"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "analysis_session_info.txt"))

message("Done. Figures: ", figdir)
message("Plot data: ", datadir)
