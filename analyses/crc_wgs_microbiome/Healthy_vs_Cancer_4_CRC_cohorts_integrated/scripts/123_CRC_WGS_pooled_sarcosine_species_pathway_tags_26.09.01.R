#!/usr/bin/env Rscript

# Add per-species sarcosine-pathway tags to the frozen pooled CRC-WGS
# sarcosine-associated-species panel.
#
# This is a display-only change. Species, log2FC, q values, directions and
# ordering are taken from the frozen complete-case pooled table. The Role
# field is converted using the same convention as the individual-cohort plots:
#   Degradation -> Deg
#   Production -> Prod
#   Degradation & Production -> Deg & Prod

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ragg)
  library(magick)
  library(png)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path.")
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_arg)))
BASE_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
RESULTS_DIR <- file.path(BASE_DIR, "results_integrated")
PROJECT_ROOT <- normalizePath(file.path(BASE_DIR, "..", "..", ".."))

INPUT_FILE <- normalizePath(file.path(
  RESULTS_DIR, "Fig2d_SupFig8_taxonomic_complete_case_26.08.12",
  "sarcosine_bacteria_diff_abundance_pooled_corrected.csv"
))
OUT_DIR <- file.path(
  RESULTS_DIR, "CRC_WGS_POOLED_SARCOSINE_SPECIES_PATHWAY_TAGS_26.09.01"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PNG <- file.path(
  OUT_DIR, "CRC_WGS_pooled_sarcosine_associated_species_with_pathway_tags.png"
)
OUT_PDF <- file.path(
  OUT_DIR, "CRC_WGS_pooled_sarcosine_associated_species_with_pathway_tags.pdf"
)
OUT_VALUES <- file.path(
  OUT_DIR, "CRC_WGS_pooled_sarcosine_associated_species_with_pathway_tags_values.csv"
)
OUT_README <- file.path(OUT_DIR, "README.md")
OUT_INPUT <- file.path(OUT_DIR, "input_md5.csv")
OUT_RENDER <- file.path(OUT_DIR, "render_spec.csv")
OUT_SESSION <- file.path(OUT_DIR, "sessionInfo.txt")
OUT_MANIFEST <- file.path(OUT_DIR, "output_md5.csv")
TRANSPORT_PNG <- file.path(
  PROJECT_ROOT, "Manuscript작업", "FigDesign&Manuscript",
  "PPT_INSERT_CRC_WGS_pooled_sarcosine_species_PATHWAY_TAGS.png"
)

COL_HEALTHY <- "#1B9E8F"
COL_CANCER <- "#C43C3C"
GROUP_COLORS <- c(Healthy = COL_HEALTHY, Cancer = COL_CANCER)
PATHWAY_COLORS <- c(
  "Deg" = "#2196F3",
  "Deg & Prod" = "#9C27B0",
  "Prod" = "#FF9800"
)

TITLE_PT <- 20
AXIS_TITLE_PT <- 17
AXIS_TEXT_PT <- 15
SPECIES_PT <- 16
SPECIES_STRIP_PT <- 14
PVAL_PT <- 15
PVAL_MM <- PVAL_PT / ggplot2::.pt
WIDTH_IN <- 11.5
HEIGHT_IN <- 8.8
DPI <- 300

required <- c(
  "Species", "Role", "Mean_Healthy", "Mean_Cancer", "log2FC", "Direction",
  "max_abs_rho", "N_Healthy", "N_Cancer", "p_adj", "cohort"
)
dat <- read.csv(INPUT_FILE, check.names = FALSE, stringsAsFactors = FALSE)
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Frozen pooled table is missing: ", paste(missing, collapse = ", "))
if (nrow(dat) != 20L || anyDuplicated(dat$Species)) {
  stop("Expected 20 unique pooled sarcosine-associated species.")
}
if (!all(dat$cohort == "pooled") || !all(dat$N_Healthy == 745L) || !all(dat$N_Cancer == 904L)) {
  stop("Unexpected pooled cohort/sample-count metadata.")
}
if (!all(dat$Direction == ifelse(dat$log2FC < 0, "Healthy-enriched", "Cancer-enriched"))) {
  stop("Direction column disagrees with log2FC sign.")
}
if (!all(dat$Role %in% c("Degradation", "Production", "Degradation & Production"))) {
  stop("Unexpected sarcosine pathway Role value.")
}

role_map <- c(
  "Degradation" = "Deg",
  "Production" = "Prod",
  "Degradation & Production" = "Deg & Prod"
)

format_prob <- function(x) {
  ifelse(
    is.na(x), "",
    ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  )
}

plot_dat <- dat |>
  mutate(
    Species_display = gsub("_", " ", Species),
    Role_label = unname(role_map[Role]),
    direction = ifelse(log2FC < 0, "Healthy", "Cancer"),
    q_label = format_prob(p_adj)
  ) |>
  arrange(log2FC)

if (anyNA(plot_dat$Role_label) || !all(plot_dat$Role_label %in% names(PATHWAY_COLORS))) {
  stop("Role-to-label conversion failed.")
}
plot_dat$Species_display <- factor(
  plot_dat$Species_display, levels = plot_dat$Species_display
)

x_pad <- max(abs(plot_dat$log2FC), na.rm = TRUE) * 0.28

theme_pub <- function() {
  theme_classic(base_family = "Arial", base_size = AXIS_TEXT_PT) +
    theme(
      plot.title = element_text(
        family = "Arial", size = TITLE_PT, face = "bold",
        hjust = 0, margin = margin(b = 8)
      ),
      axis.title = element_text(family = "Arial", size = AXIS_TITLE_PT),
      axis.text = element_text(family = "Arial", size = AXIS_TEXT_PT, colour = "black"),
      legend.text = element_text(family = "Arial", size = AXIS_TEXT_PT),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank()
    )
}

p <- ggplot(plot_dat, aes(log2FC, Species_display, fill = direction)) +
  geom_col(width = 0.68, alpha = 0.88, colour = "black", linewidth = 0.25) +
  geom_text(
    aes(
      x = log2FC + ifelse(log2FC >= 0, x_pad * 0.08, -x_pad * 0.08),
      label = q_label, hjust = ifelse(log2FC >= 0, 0, 1)
    ),
    family = "Arial", size = PVAL_MM
  ) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.45) +
  geom_point(
    data = plot_dat,
    aes(x = Inf, y = Species_display, colour = Role_label),
    inherit.aes = FALSE, shape = 15, size = 4.2,
    show.legend = c(colour = TRUE, fill = FALSE)
  ) +
  geom_text(
    data = plot_dat,
    aes(x = Inf, y = Species_display, label = Role_label, colour = Role_label),
    inherit.aes = FALSE, hjust = -0.36, family = "Arial",
    size = SPECIES_STRIP_PT / ggplot2::.pt,
    fontface = "bold", show.legend = FALSE
  ) +
  scale_fill_manual(values = GROUP_COLORS, name = NULL) +
  scale_colour_manual(values = PATHWAY_COLORS, name = "Sarcosine pathway") +
  scale_x_continuous(
    breaks = c(-2, 0, 2, 4),
    expand = expansion(mult = c(0.25, 0.28))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Sarcosine-associated species",
    x = expression(log[2]~fold~change), y = NULL
  ) +
  theme_pub() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(
      family = "Arial", size = AXIS_TEXT_PT, face = "bold"
    ),
    axis.text.y = element_text(
      family = "Arial", size = SPECIES_PT, colour = "black"
    ),
    plot.margin = margin(10, 145, 10, 10)
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(colour = NA)),
    colour = guide_legend(
      order = 2, override.aes = list(shape = 15, size = 4.2)
    )
  )

write.csv(plot_dat, OUT_VALUES, row.names = FALSE, quote = TRUE)
write.csv(
  data.frame(
    Input = "Frozen pooled complete-case sarcosine-associated species table",
    File = INPUT_FILE,
    MD5 = unname(tools::md5sum(INPUT_FILE)),
    stringsAsFactors = FALSE
  ),
  OUT_INPUT, row.names = FALSE, quote = TRUE
)
write.csv(
  data.frame(
    width_in = WIDTH_IN, height_in = HEIGHT_IN, dpi = DPI,
    width_px = as.integer(WIDTH_IN * DPI), height_px = as.integer(HEIGHT_IN * DPI),
    background = "white", colorspace = "sRGB", bit_depth = 8,
    stringsAsFactors = FALSE
  ),
  OUT_RENDER, row.names = FALSE, quote = TRUE
)

ragg::agg_png(
  filename = OUT_PNG, width = WIDTH_IN, height = HEIGHT_IN,
  units = "in", res = DPI, background = "white"
)
print(p)
dev.off()
ggsave(OUT_PDF, p, width = WIDTH_IN, height = HEIGHT_IN, device = cairo_pdf, bg = "white")

SRGB_PROFILE <- "/System/Library/ColorSync/Profiles/sRGB Profile.icc"
if (!file.exists(SRGB_PROFILE)) stop("System sRGB profile not found.")

standardize_png <- function(path, dpi) {
  img <- magick::image_read(path)
  img <- magick::image_background(img, "white", flatten = TRUE)
  img <- magick::image_convert(
    img, format = "png", colorspace = "sRGB", depth = 8, type = "TrueColor"
  )
  img <- magick::image_strip(img)
  magick::image_write(img, path = path, format = "png", depth = 8)
  before <- png::readPNG(path, native = FALSE, info = TRUE)
  tmp <- tempfile(fileext = ".png")
  sips_output <- system2(
    "/usr/bin/sips",
    args = c(
      "-s", "dpiWidth", as.character(dpi),
      "-s", "dpiHeight", as.character(dpi),
      "-e", shQuote(SRGB_PROFILE), shQuote(path), "--out", shQuote(tmp)
    ),
    stdout = TRUE, stderr = TRUE
  )
  if (!file.exists(tmp) || file.info(tmp)$size <= 5000) {
    stop("sips failed: ", paste(sips_output, collapse = "\n"))
  }
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not install standardized PNG.")
  unlink(tmp)
  after <- png::readPNG(path, native = FALSE, info = TRUE)
  if (!identical(dim(before), dim(after)) || !identical(as.numeric(before), as.numeric(after))) {
    stop("Embedding sRGB/dpi metadata changed decoded pixels.")
  }
}

standardize_png(OUT_PNG, DPI)
if (!file.copy(OUT_PNG, TRANSPORT_PNG, overwrite = TRUE)) {
  stop("Could not create PowerPoint transport PNG.")
}
for (path in c(OUT_PNG, OUT_PDF, OUT_VALUES, OUT_README, OUT_INPUT,
               OUT_RENDER, OUT_SESSION, TRANSPORT_PNG)) {
  if (file.exists(path)) {
    invisible(system2("/usr/bin/xattr", args = c("-c", shQuote(path)), stdout = TRUE, stderr = TRUE))
    invisible(system2("/bin/chmod", args = c("644", shQuote(path)), stdout = TRUE, stderr = TRUE))
    invisible(system2("/usr/bin/chflags", args = c("nohidden", shQuote(path)), stdout = TRUE, stderr = TRUE))
  }
}
if (!identical(unname(tools::md5sum(OUT_PNG)), unname(tools::md5sum(TRANSPORT_PNG)))) {
  stop("PowerPoint transport PNG is not byte-identical to the analytical PNG.")
}

role_counts <- table(plot_dat$Role_label)
writeLines(
  c(
    "# CRC WGS pooled sarcosine-associated species with pathway tags",
    "",
    "- Display-only revision of the frozen pooled complete-case 20-species panel.",
    "- No species, log2FC, q value, direction or ordering was changed.",
    "- Per-species tags use the same mapping and colours as the individual-cohort panels:",
    "  - `Deg`: Degradation, blue `#2196F3`",
    "  - `Prod`: Production, orange `#FF9800`",
    "  - `Deg & Prod`: Degradation & Production, purple `#9C27B0`",
    sprintf("- Role counts: Deg=%d, Deg & Prod=%d, Prod=%d.", role_counts[["Deg"]], role_counts[["Deg & Prod"]], role_counts[["Prod"]]),
    "- Group colours remain Healthy teal `#1B9E8F` and Cancer red `#C43C3C`.",
    sprintf("- Analytical PNG: `%s`", OUT_PNG),
    sprintf("- PowerPoint transport PNG: `%s`", TRANSPORT_PNG),
    "- The pre-existing selected figure was not overwritten."
  ),
  OUT_README
)
capture.output(sessionInfo(), file = OUT_SESSION)

manifest_files <- c(OUT_PNG, OUT_PDF, OUT_VALUES, OUT_README, OUT_INPUT, OUT_RENDER, OUT_SESSION)
write.csv(
  data.frame(
    File = basename(manifest_files),
    Bytes = file.info(manifest_files)$size,
    MD5 = unname(tools::md5sum(manifest_files)),
    stringsAsFactors = FALSE
  ),
  OUT_MANIFEST, row.names = FALSE, quote = TRUE
)

cat("Wrote analytical PNG: ", OUT_PNG, "\n", sep = "")
cat("Wrote PowerPoint PNG: ", TRANSPORT_PNG, "\n", sep = "")
cat("Wrote plotted values: ", OUT_VALUES, "\n", sep = "")
