#!/usr/bin/env Rscript
# Display-only redraw of the frozen pooled CRC species panel.
# No filtering, statistical testing, normalization, or role reassignment.
# Workflow: inspect and cross-check both frozen sources -> arrange the original
# descending effect-size order -> render -> save displayed values and provenance.
set.seed(42) # Fixed for reproducibility; this deterministic rendering uses no RNG.
suppressPackageStartupMessages({ library(grid); library(here) })
here::i_am("공공_Metabolomics&Metagenomics_분석모음/HGMT_CRC_WGS-Healthy_vs_Cancer/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/125_CRC_WGS_species_readable_26.09.06.R")
base <- file.path("공공_Metabolomics&Metagenomics_분석모음", "HGMT_CRC_WGS-Healthy_vs_Cancer", "Healthy_vs_Cancer_4_CRC_cohorts_integrated")
input_rel <- file.path(base, "results_integrated", "Fig2d_SupFig8_taxonomic_complete_case_26.08.12", "sarcosine_bacteria_diff_abundance_pooled_corrected.csv")
tagged_rel <- file.path(base, "results_integrated", "CRC_WGS_POOLED_SARCOSINE_SPECIES_PATHWAY_TAGS_26.09.01", "CRC_WGS_pooled_sarcosine_associated_species_with_pathway_tags_values.csv")
out_rel <- file.path(base, "results_integrated", "CRC_WGS_SPECIES_READABLE_26.09.06")
out_dir <- here::here(out_rel)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source <- read.csv(here::here(input_rel), check.names = FALSE)
tagged <- read.csv(here::here(tagged_rel), check.names = FALSE)
required <- c("Species", "Role", "log2FC", "p_adj", "Direction", "N_Healthy", "N_Cancer", "cohort")
stopifnot(all(required %in% names(source)), nrow(source) == 20L,
          !anyDuplicated(source$Species), !anyNA(source[required]),
          all(source$N_Healthy == 745L), all(source$N_Cancer == 904L),
          all(source$cohort == "pooled"), all(is.finite(source$log2FC)),
          all(source$p_adj >= 0 & source$p_adj <= 1),
          setequal(source$Species, tagged$Species))
matched <- tagged[match(source$Species, tagged$Species), ]
stopifnot(identical(source[required], matched[required] |> `rownames<-`(NULL)))
roles <- c("Degradation", "Production", "Degradation & Production")
stopifnot(all(source$Role %in% roles),
          all(source$Direction == ifelse(source$log2FC < 0, "Healthy-enriched", "Cancer-enriched")))

# The source Role field records pathway association, not experimental activity.
d <- source[order(source$log2FC, decreasing = TRUE), ]
rownames(d) <- NULL
d$Species_display <- gsub("_", " ", d$Species)
d$Degradation_associated <- d$Role %in% c("Degradation", "Degradation & Production")
d$Production_associated <- d$Role %in% c("Production", "Degradation & Production")
# Retain all q values, including nonsignificant results. Stored zero is printed
# as <0.001 (the original display convention), never as an exact probability 0.
d$q_label <- ifelse(d$p_adj < 0.001, "<0.001", sprintf("%.3f", d$p_adj))
d$wrapped_label <- vapply(d$Species_display, function(s) paste(strwrap(s, width = 34), collapse = "\n"), "")
line_count <- lengths(strsplit(d$wrapped_label, "\n"))
row_height <- ifelse(line_count == 1L, 7.0, 10.4) # mm; space for 11 pt labels.
width_mm <- 200
header_mm <- 34
height_mm <- header_mm + sum(row_height) + 26
y_top <- height_mm - header_mm
row_y <- y_top - cumsum(row_height) + row_height / 2
body_bottom <- y_top - sum(row_height)
col_healthy <- "#1B9E8F"
col_crc <- "#C43C3C"
col_deg <- "#237CB3"
col_prod <- "#C88616"
axis_min <- -1.2
axis_max <- 4.7
chart_left <- 123
chart_right <- 178
map_x <- function(v) chart_left + (v - axis_min) / (axis_max - axis_min) * (chart_right - chart_left)
stopifnot(all(d$log2FC > axis_min & d$log2FC < axis_max))

draw <- function() {
  grid.newpage()
  pushViewport(viewport(xscale = c(0, width_mm), yscale = c(0, height_mm)))
  txt <- function(label, x, y, size = 11, colour = "#20252A", face = "plain", just = "centre") {
    grid.text(label, x, y, default.units = "native", just = just,
              gp = gpar(fontfamily = "Arial", fontsize = size, col = colour,
                        fontface = face, lineheight = 1.05))
  }
  line <- function(x0, y0, x1, y1, colour = "#DCE1E4", lwd = 0.5) {
    grid.segments(x0, y0, x1, y1, default.units = "native", gp = gpar(col = colour, lwd = lwd))
  }
  txt("Sarcosine-associated species", 6, height_mm - 7, 15, face = "bold", just = "left")
  txt("Pooled CRC cohorts  |  Healthy n = 745; CRC n = 904", 6, height_mm - 14, 9.5,
      colour = "#606970", just = "left")
  txt("Species", 6, height_mm - 27, 11, face = "bold", just = "left")
  txt("Pathway association", 100, height_mm - 23, 10, face = "bold")
  txt("Deg.", 92, height_mm - 29, 10.5, col_deg, "bold")
  txt("Prod.", 108, height_mm - 29, 10.5, col_prod, "bold")
  txt("Healthy higher", 137, height_mm - 23, 8.7, col_healthy, "bold")
  txt("CRC higher", 171, height_mm - 23, 8.7, col_crc, "bold")
  txt("q value", 190, height_mm - 27, 10.5, face = "bold")
  line(6, y_top + 0.6, 195, y_top + 0.6, "#AEB7BD", 0.8)
  for (i in seq_len(nrow(d))) {
    if (i %% 2L == 0L) {
      grid.rect(x = 100.5, y = row_y[i], width = 189, height = row_height[i],
                default.units = "native", gp = gpar(fill = "#F4F6F7", col = NA))
    }
  }
  for (tick in c(-1, 1, 2, 3, 4)) {
    line(map_x(tick), body_bottom, map_x(tick), y_top, "#E2E6E9", 0.45)
  }
  line(map_x(0), body_bottom, map_x(0), y_top, "#7B858C", 0.8)
  for (i in seq_len(nrow(d))) {
    txt(d$wrapped_label[i], 6, row_y[i], 11, face = "italic", just = "left")
    associations <- c(d$Degradation_associated[i], d$Production_associated[i])
    for (j in 1:2) {
      grid.circle(x = c(92, 108)[j], y = row_y[i], r = unit(1.35, "mm"),
                  default.units = "native", gp = gpar(
                    fill = if (associations[j]) c(col_deg, col_prod)[j] else NA,
                    col = if (associations[j]) c(col_deg, col_prod)[j] else "#D0D6DA", lwd = 0.8))
    }
    effect_x <- map_x(d$log2FC[i])
    zero_x <- map_x(0)
    grid.rect(x = (effect_x + zero_x) / 2, y = row_y[i],
              width = abs(effect_x - zero_x), height = 3.1,
              default.units = "native", gp = gpar(col = NA,
                fill = if (d$log2FC[i] < 0) col_healthy else col_crc))
    txt(d$q_label[i], 190, row_y[i], 10)
  }
  line(6, body_bottom, 195, body_bottom, "#AEB7BD", 0.8)
  for (tick in -1:4) {
    line(map_x(tick), body_bottom, map_x(tick), body_bottom - 1.2, "#667078", 0.65)
    txt(as.character(tick), map_x(tick), body_bottom - 4, 9)
  }
  txt(expression(log[2]*"(CRC / Healthy)"), (chart_left + chart_right) / 2, body_bottom - 10, 10)
  txt("Filled dots: degradation (Deg.) and/or production (Prod.) pathway association.",
      6, 9, 8.5, colour = "#4B555C", just = "left")
  txt("KO-based associations; q values: Benjamini-Hochberg adjusted P values.",
      6, 4.5, 8.5, colour = "#4B555C", just = "left")
  popViewport()
}

stem <- "Fig2_Sarcosine_species_26.09.06"
ragg::agg_png(file.path(out_dir, paste0(stem, ".png")), width = width_mm,
              height = height_mm, units = "mm", res = 600, background = "white")
draw(); invisible(dev.off())
if (capabilities("aqua")) {
  # Native macOS PDF device embeds Arial without requiring an XQuartz install.
  quartz(type = "pdf", file = file.path(out_dir, paste0(stem, ".pdf")),
         width = width_mm / 25.4, height = height_mm / 25.4,
         family = "Arial", bg = "white")
} else {
  cairo_pdf(file.path(out_dir, paste0(stem, ".pdf")), width = width_mm / 25.4,
            height = height_mm / 25.4, family = "Arial", onefile = TRUE, bg = "white")
}
draw(); invisible(dev.off())
write.csv(d, file.path(out_dir, "displayed_values.csv"), row.names = FALSE)
source_manifest <- data.frame(
  file = c(input_rel, tagged_rel),
  sha256 = vapply(c(input_rel, tagged_rel), function(p) digest::digest(file = here::here(p), algo = "sha256"), "")
)
write.csv(source_manifest, file.path(out_dir, "source_sha256.csv"), row.names = FALSE)
writeLines(c(
  "# CRC sarcosine-associated species: display-only redraw (2026-09-06)",
  "", "All 20 taxa, log2 fold changes, BH-adjusted P values and Role labels are retained from the frozen sources.",
  "No statistics were recomputed; no samples or taxa were removed. Display order is descending log2FC, as in the original panel.",
  "Healthy n = 745 and CRC n = 904 are the species-table counts, not the separate KO-table counts.",
  "Filled dots display the existing Degradation / Production / Degradation & Production Role field.",
  "These are KO-based pathway associations, not direct measurements of production or degradation activity.",
  "A hollow dot means the corresponding association is not assigned in the source; it does not establish biological absence.",
  "Bar colors show the direction of the mean abundance difference, regardless of q-value significance.",
  "Upstream p_adj = 0 is retained in the values file and displayed as <0.001, following the source convention.",
  "", paste0("Export: white-background 600 dpi PNG and vector PDF, ", width_mm, " x ", height_mm, " mm; Arial; species labels 11 pt."),
  "The figure has no panel letter so the author can assign it in PowerPoint.",
  "", "## Reproduction", "Run the accompanying scripts/125_CRC_WGS_species_readable_26.09.06.R from within the project tree.",
  "Uses relative project paths via here. This is deterministic; seed 42 is recorded but unused.",
  "See source_sha256.csv, displayed_values.csv and sessionInfo.txt for provenance and software details."
), file.path(out_dir, "README.md"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
transport <- here::here("Manuscript작업", "FigDesign&Manuscript")
stopifnot(dir.exists(transport))
for (ext in c("png", "pdf")) {
  name <- paste0(stem, ".", ext)
  stopifnot(file.copy(file.path(out_dir, name), file.path(transport, name), overwrite = TRUE))
}
cat("Rendered", nrow(d), "species;", width_mm, "x", height_mm, "mm.\n")
print(table(d$Role))
cat("Outputs:", file.path(out_rel, stem), "\n")
