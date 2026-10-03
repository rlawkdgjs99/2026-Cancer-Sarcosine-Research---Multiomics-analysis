#!/usr/bin/env Rscript
# Display-only redraw of the frozen pooled CRC species panel.
# No filtering, statistical testing, normalization, or role reassignment.
# Workflow: inspect and cross-check both frozen sources -> arrange the original
# descending effect-size order -> render -> save displayed values and provenance.
set.seed(42) # Fixed for reproducibility; this deterministic rendering uses no RNG.
suppressPackageStartupMessages({ library(grid); library(here) })
here::i_am("공공_Metabolomics&Metagenomics_분석모음/HGMT_CRC_WGS-Healthy_vs_Cancer/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/126_CRC_WGS_species_bar_style_26.09.06.R")
base <- file.path("공공_Metabolomics&Metagenomics_분석모음", "HGMT_CRC_WGS-Healthy_vs_Cancer", "Healthy_vs_Cancer_4_CRC_cohorts_integrated")
input_rel <- file.path(base, "results_integrated", "Fig2d_SupFig8_taxonomic_complete_case_26.08.12", "sarcosine_bacteria_diff_abundance_pooled_corrected.csv")
tagged_rel <- file.path(base, "results_integrated", "CRC_WGS_POOLED_SARCOSINE_SPECIES_PATHWAY_TAGS_26.09.01", "CRC_WGS_pooled_sarcosine_associated_species_with_pathway_tags_values.csv")
out_rel <- file.path(base, "results_integrated", "CRC_WGS_SPECIES_BAR_STYLE_26.09.06")
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
role_labels <- c("Degradation" = "Deg", "Production" = "Prod",
                 "Degradation & Production" = "Deg & Prod")
role_colours <- c("Degradation" = "#2878C8", "Production" = "#CF9B00",
                  "Degradation & Production" = "#AD3FC4")
d$Role_label <- unname(role_labels[d$Role])
width_mm <- 210
row_height <- 6.2 # mm; single-line labels, consistent spacing.
header_mm <- 24
height_mm <- header_mm + nrow(d) * row_height + 35
y_top <- height_mm - header_mm
row_y <- y_top - (seq_len(nrow(d)) - 0.5) * row_height
body_bottom <- y_top - nrow(d) * row_height
col_healthy <- "#1B9E8F"
col_crc <- "#C43C3C"
axis_min <- -2.05
axis_max <- 5.5
chart_left <- 91
chart_right <- 176
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
  line <- function(x0, y0, x1, y1, colour = "#E7EBEE", lwd = 0.5) {
    grid.segments(x0, y0, x1, y1, default.units = "native", gp = gpar(col = colour, lwd = lwd))
  }
  square <- function(x, y, colour, size = 1.7) {
    grid.rect(x, y, width = size, height = size, default.units = "native",
              gp = gpar(fill = colour, col = NA))
  }
  # Fail explicitly if a species label would extend beyond the left margin.
  label_widths <- vapply(d$Species_display, function(label) {
    convertWidth(grobWidth(textGrob(label, gp = gpar(fontfamily = "Arial",
                 fontface = "italic", fontsize = 11.5))), "mm", valueOnly = TRUE)
  }, 0.0)
  stopifnot(max(label_widths) < 83)
  txt("Sarcosine-associated species", width_mm / 2, height_mm - 6.5, 15, face = "bold")
  txt("Pooled CRC cohorts  |  Healthy n = 745; CRC n = 904", width_mm / 2,
      height_mm - 13, 9.5, colour = "#667078")
  txt("Pathway\nassociation", 193, y_top + 5, 9, colour = "#5D656C")
  for (tick in c(-2, -1, 1, 2, 3, 4, 5)) {
    line(map_x(tick), body_bottom, map_x(tick), y_top)
  }
  line(map_x(0), body_bottom, map_x(0), y_top, "#747F86", 0.85)
  for (i in seq_len(nrow(d))) {
    txt(d$Species_display[i], 87, row_y[i], 11.5, face = "italic", just = "right")
    effect_x <- map_x(d$log2FC[i])
    zero_x <- map_x(0)
    grid.rect(x = (effect_x + zero_x) / 2, y = row_y[i],
              width = abs(effect_x - zero_x), height = 3.65,
              default.units = "native", gp = gpar(col = NA,
                fill = if (d$log2FC[i] < 0) col_healthy else col_crc))
    txt(d$q_label[i], effect_x + if (d$log2FC[i] < 0) -1.1 else 1.1,
        row_y[i], 9.5, colour = "#434B50",
        just = if (d$log2FC[i] < 0) "right" else "left")
    role_colour <- unname(role_colours[d$Role[i]])
    square(182, row_y[i], role_colour)
    txt(d$Role_label[i], 185, row_y[i], 10, role_colour, just = "left")
  }
  line(chart_left, body_bottom, chart_right, body_bottom, "#5B656C", 0.8)
  for (tick in -2:5) {
    line(map_x(tick), body_bottom, map_x(tick), body_bottom - 1.3, "#5B656C", 0.7)
    txt(as.character(tick), map_x(tick), body_bottom - 4.2, 9.5)
  }
  txt(expression(log[2]*" fold change (CRC / Healthy)"),
      (chart_left + chart_right) / 2, body_bottom - 10, 10.5)
  square(83, body_bottom - 16, col_crc, 2.2)
  txt("CRC higher", 86, body_bottom - 16, 10, just = "left")
  square(116, body_bottom - 16, col_healthy, 2.2)
  txt("Healthy higher", 119, body_bottom - 16, 10, just = "left")
  square(43, body_bottom - 23, role_colours[["Degradation"]])
  txt("Deg: degradation", 46, body_bottom - 23, 9.5, just = "left")
  square(92, body_bottom - 23, role_colours[["Degradation & Production"]])
  txt("Deg & Prod: both", 95, body_bottom - 23, 9.5, just = "left")
  square(144, body_bottom - 23, role_colours[["Production"]])
  txt("Prod: production", 147, body_bottom - 23, 9.5, just = "left")
  txt("Labels indicate KO-based pathway associations; numbers are BH-adjusted P values (q).",
      width_mm / 2, 4.5, 8.5, colour = "#596269")
  popViewport()
}

stem <- "Fig2_Sarcosine_species_bar_style_26.09.06"
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
  "# CRC sarcosine-associated species: display-only redraw, reference bar-plot style (2026-09-06)",
  "", "All 20 taxa, log2 fold changes, BH-adjusted P values and Role labels are retained from the frozen sources.",
  "No statistics were recomputed; no samples or taxa were removed. Display order is descending log2FC, as in the original panel.",
  "Healthy n = 745 and CRC n = 904 are the species-table counts, not the separate KO-table counts.",
  "Colored square markers and right-side text labels display the existing Degradation / Production / Degradation & Production Role field.",
  "These are KO-based pathway associations, not direct measurements of production or degradation activity.",
  "Bar colors show the direction of the mean abundance difference, regardless of q-value significance.",
  "Upstream p_adj = 0 is retained in the values file and displayed as <0.001, following the source convention.",
  "", paste0("Export: white-background 600 dpi PNG and vector PDF, ", width_mm, " x ", height_mm, " mm; Arial; species labels 11.5 pt."),
  "The figure has no panel letter so the author can assign it in PowerPoint.",
  "", "## Reproduction", "Run the accompanying scripts/126_CRC_WGS_species_bar_style_26.09.06.R from within the project tree.",
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
