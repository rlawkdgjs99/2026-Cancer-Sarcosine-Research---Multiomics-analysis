#!/usr/bin/env Rscript

# Candidate CRC-WGS targeted seven-KO lollipop plot.
#
# Scientific scope
#   - Input is the frozen pooled four-cohort targeted sarcosine-KO table.
#   - All seven detected predefined KOs are shown, including four KOs below the
#     original >=10% overall-prevalence threshold.
#   - Wilcoxon P values and BH-adjusted q values are read from the frozen table.
#     These q values were adjusted across the seven targeted KOs, and therefore
#     are not mixed with the genome-wide 7,106-KO BH values used in the older
#     three-KO panel.
#   - log2FC is reconstructed only from the frozen group means using the same
#     pseudocount (1e-8) as the pooled KO differential-abundance pipeline.
#   - Open circles explicitly flag overall prevalence <10%.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtext)
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
SHARED_THEME <- normalizePath(file.path(
  BASE_DIR, "..", "..", "_shared", "theme_nc_26.08.18.R"
))
source(SHARED_THEME)

INPUT_TARGET <- normalizePath(file.path(
  RESULTS_DIR, "sarcosine", "sarcosine_KO_comparison_filtered_pooled.csv"
))
INPUT_GENOME <- normalizePath(file.path(
  RESULTS_DIR, "kegg", "diff_KO_abundance_pooled.csv"
))
OUT_DIR <- file.path(RESULTS_DIR, "CRC_WGS_SARCOSINE_7KO_LOLLIPOP_26.09.01")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PNG <- file.path(OUT_DIR, "CRC_WGS_sarcosine_7KO_lollipop.png")
OUT_PDF <- file.path(OUT_DIR, "CRC_WGS_sarcosine_7KO_lollipop.pdf")
OUT_VALUES <- file.path(OUT_DIR, "CRC_WGS_sarcosine_7KO_lollipop_values.csv")
OUT_INPUTS <- file.path(OUT_DIR, "input_md5.csv")
OUT_RENDER <- file.path(OUT_DIR, "render_spec.csv")
OUT_SESSION <- file.path(OUT_DIR, "sessionInfo.txt")
OUT_README <- file.path(OUT_DIR, "README.md")
OUT_MANIFEST <- file.path(OUT_DIR, "output_md5.csv")
TRANSPORT_PNG <- file.path(
  PROJECT_ROOT, "Manuscript작업", "FigDesign&Manuscript",
  "PPT_INSERT_CRC_WGS_sarcosine_7KO_lollipop.png"
)

EXPECTED_KO <- c("K00301", "K00302", "K00303", "K00305", "K00306", "K00315", "K08688")
EXPECTED_ROLE <- c(rep("Degradation", 5), rep("Production", 2))
PSEUDOCOUNT <- 1e-8
PREVALENCE_THRESHOLD <- 10
WIDTH_IN <- 4.50
HEIGHT_IN <- 4.45
DPI <- 600
EXPECTED_PX <- c(width = as.integer(WIDTH_IN * DPI), height = as.integer(HEIGHT_IN * DPI))

target <- read.csv(INPUT_TARGET, check.names = FALSE, stringsAsFactors = FALSE)
genome <- read.csv(INPUT_GENOME, check.names = FALSE, stringsAsFactors = FALSE)

required_target <- c(
  "KO", "Enzyme", "EC", "Role", "Mean_Healthy", "Mean_Cancer",
  "p_value", "p_adj", "Prev_Overall", "Prev_Healthy", "Prev_Cancer",
  "Pass_Prevalence"
)
required_genome <- c("KO", "Mean_Healthy", "Mean_Cancer", "log2FC", "p_value", "p_adj")
if (length(setdiff(required_target, names(target)))) {
  stop("Targeted KO table is missing: ", paste(setdiff(required_target, names(target)), collapse = ", "))
}
if (length(setdiff(required_genome, names(genome)))) {
  stop("Genome-wide KO table is missing: ", paste(setdiff(required_genome, names(genome)), collapse = ", "))
}
if (nrow(target) != 7L || anyDuplicated(target$KO)) {
  stop("Expected exactly seven unique detected KOs in the frozen targeted table.")
}
target <- target[match(EXPECTED_KO, target$KO), , drop = FALSE]
if (anyNA(target$KO) || !identical(target$KO, EXPECTED_KO)) {
  stop("The frozen targeted table does not contain the expected seven KOs.")
}
if (!identical(target$Role, EXPECTED_ROLE)) stop("Unexpected Production/Degradation role mapping.")

# Verify, but do not replace, the stored targeted-family BH q values.
recomputed_q_for_check <- p.adjust(target$p_value, method = "BH")
if (!isTRUE(all.equal(target$p_adj, recomputed_q_for_check, tolerance = 1e-14))) {
  stop("Stored targeted-family BH q values do not reproduce from the seven frozen P values.")
}

target$log2FC_Cancer_vs_Healthy <- log2(
  (target$Mean_Cancer + PSEUDOCOUNT) / (target$Mean_Healthy + PSEUDOCOUNT)
)

# The three prevalence-qualified rows must agree exactly with the genome-wide table.
qualified <- target$KO[target$Pass_Prevalence]
if (!identical(qualified, c("K00302", "K00303", "K08688"))) {
  stop("Unexpected set of KOs passing the saved >=10% prevalence filter.")
}
g <- genome[match(qualified, genome$KO), required_genome, drop = FALSE]
t <- target[match(qualified, target$KO), , drop = FALSE]
if (anyNA(g$KO) || max(abs(g$Mean_Healthy - t$Mean_Healthy)) > 1e-15 ||
    max(abs(g$Mean_Cancer - t$Mean_Cancer)) > 1e-15 ||
    max(abs(g$p_value - t$p_value)) > 1e-15 ||
    max(abs(g$log2FC - t$log2FC_Cancer_vs_Healthy)) > 1e-14) {
  stop("Shared values disagree between the targeted and genome-wide frozen tables.")
}
if (any(target$KO[!target$Pass_Prevalence] %in% genome$KO)) {
  stop("A <10% prevalence KO unexpectedly appears in the prevalence-filtered genome-wide table.")
}

short_name <- c(
  K00301 = "Sarcosine oxidase (monomeric)",
  K00302 = "<i>soxA</i>",
  K00303 = "<i>soxB</i>",
  K00305 = "<i>soxG</i>",
  K00306 = "PIPOX",
  K00315 = "DMGDH",
  K08688 = "Creatinase"
)

plot_dat <- data.frame(
  KO = target$KO,
  Enzyme = target$Enzyme,
  EC = target$EC,
  Role = target$Role,
  Mean_Healthy = target$Mean_Healthy,
  Mean_Cancer = target$Mean_Cancer,
  log2FC_Cancer_vs_Healthy = target$log2FC_Cancer_vs_Healthy,
  wilcox_p = target$p_value,
  targeted_7KO_BH_q = target$p_adj,
  Prevalence_percent = target$Prev_Overall,
  Prevalence_Healthy_percent = target$Prev_Healthy,
  Prevalence_Cancer_percent = target$Prev_Cancer,
  Pass_10pct_prevalence = target$Pass_Prevalence,
  stringsAsFactors = FALSE
)
plot_dat$Status <- ifelse(
  plot_dat$targeted_7KO_BH_q < 0.05 & plot_dat$log2FC_Cancer_vs_Healthy < 0,
  "Healthy enriched",
  ifelse(
    plot_dat$targeted_7KO_BH_q < 0.05 & plot_dat$log2FC_Cancer_vs_Healthy > 0,
    "Cancer enriched", "Not significant"
  )
)
plot_dat$Prevalence_display <- ifelse(
  plot_dat$Pass_10pct_prevalence,
  sprintf("prevalence %.1f%%", plot_dat$Prevalence_percent),
  sprintf("prevalence %.1f%% (&lt;10%%)", plot_dat$Prevalence_percent)
)
plot_dat$y <- rev(seq_len(nrow(plot_dat)))

format_q <- function(x) {
  ifelse(x < 0.001, sprintf("q = %.2e", x), sprintf("q = %.3g", x))
}
plot_dat$q_label <- format_q(plot_dat$targeted_7KO_BH_q)

EXPECTED_LOG2FC <- c(
  -2.343561508170084, 0.432516736783865, -0.581999817446284,
  2.358075889973686, 0.662554070790937, 1.843046933267022,
  0.892971415794552
)
EXPECTED_STATUS <- c(
  "Not significant", "Not significant", "Healthy enriched",
  "Cancer enriched", "Not significant", "Not significant", "Cancer enriched"
)
if (!isTRUE(all.equal(plot_dat$log2FC_Cancer_vs_Healthy, EXPECTED_LOG2FC, tolerance = 1e-13)) ||
    !identical(plot_dat$Status, EXPECTED_STATUS)) {
  stop("Derived seven-KO effect values/statuses differ from the frozen expected values.")
}

axis_labels <- setNames(
  sprintf(
    "%s (%s)<br><span style='font-size:5.2pt'>%s · %s</span>",
    unname(short_name[plot_dat$KO]), plot_dat$KO, plot_dat$Role,
    plot_dat$Prevalence_display
  ),
  as.character(plot_dat$y)
)

status_colours <- c(
  "Healthy enriched" = COL_HEALTHY,
  "Cancer enriched" = COL_CANCER,
  "Not significant" = "#8A8A8A"
)
high_prev <- plot_dat[plot_dat$Pass_10pct_prevalence, , drop = FALSE]
low_prev <- plot_dat[!plot_dat$Pass_10pct_prevalence, , drop = FALSE]

p <- ggplot(plot_dat, aes(x = log2FC_Cancer_vs_Healthy, y = y)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 2.5,
           fill = "#F7F7F7", colour = NA) +
  geom_hline(yintercept = 2.5, linewidth = pt_lw(0.7), colour = "#C8C8C8") +
  geom_vline(xintercept = 0, colour = "#2D2D2D", linewidth = NC_AXIS_LW) +
  geom_segment(
    aes(x = 0, xend = log2FC_Cancer_vs_Healthy, yend = y, colour = Status),
    linewidth = pt_lw(1.35), lineend = "round"
  ) +
  geom_point(
    data = high_prev, aes(fill = Status, size = Prevalence_percent), shape = 21,
    colour = "#2D2D2D", stroke = pt_lw(1.0)
  ) +
  geom_point(
    data = low_prev, aes(colour = Status, size = Prevalence_percent), shape = 21,
    fill = "white", stroke = pt_lw(1.3), show.legend = FALSE
  ) +
  geom_text(
    aes(label = q_label, colour = Status),
    position = position_nudge(y = 0.27),
    family = "Arial", size = mm_text(5.6), vjust = 0.5,
    show.legend = FALSE
  ) +
  annotate(
    "text", x = -2.08, y = 7.67, label = "Healthy-enriched",
    family = "Arial", fontface = "bold", size = mm_text(5.8),
    colour = COL_HEALTHY
  ) +
  annotate(
    "text", x = 1.98, y = 7.67, label = "Cancer-enriched",
    family = "Arial", fontface = "bold", size = mm_text(5.8),
    colour = COL_CANCER
  ) +
  scale_colour_manual(values = status_colours, guide = "none") +
  scale_fill_manual(
    values = status_colours,
    name = "Direction",
    breaks = c("Healthy enriched", "Cancer enriched", "Not significant"),
    labels = c("Higher in Healthy", "Higher in Cancer", "Not significant")
  ) +
  scale_size_continuous(
    name = "Overall prevalence (%)",
    range = c(2.1, 5.4), breaks = c(10, 50, 90), limits = c(0, 100)
  ) +
  scale_x_continuous(
    breaks = c(-2, -1, 0, 1, 2),
    labels = c("−2", "−1", "0", "1", "2"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = plot_dat$y, labels = axis_labels,
    limits = c(0.48, 7.82), expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(-2.73, 2.73), clip = "off") +
  labs(
    x = expression(log[2]~fold~change~"(Cancer/Healthy)"), y = NULL,
    caption = "BH correction across 7 targeted KOs; open circles: prevalence <10%"
  ) +
  guides(
    fill = guide_legend(
      order = 1, nrow = 1, title.position = "top",
      override.aes = list(shape = 21, size = 3.2, colour = "#2D2D2D")
    ),
    size = guide_legend(
      order = 2, nrow = 1, title.position = "top",
      override.aes = list(shape = 21, fill = "#A8A8A8", colour = "#2D2D2D")
    )
  ) +
  theme_nc() +
  theme(
    axis.text.y = ggtext::element_markdown(
      family = "Arial", size = NC_TICK_PT, colour = "black",
      lineheight = 0.92, hjust = 1, margin = margin(r = 4)
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title = element_text(
      family = "Arial", face = "bold", size = NC_LEGEND_PT,
      colour = "#222222", hjust = 0.5
    ),
    legend.text = element_text(
      family = "Arial", size = NC_LEGEND_PT, colour = "#222222"
    ),
    legend.key.width = unit(8, "pt"),
    legend.spacing.y = unit(0.5, "pt"),
    legend.box.spacing = unit(1.5, "pt"),
    plot.caption = element_text(
      family = "Arial", size = 5.2, colour = "#4D4D4D", hjust = 0.5,
      margin = margin(t = 2, b = 0)
    ),
    plot.margin = margin(t = 2, r = 3, b = 1, l = 1)
  )

write.csv(plot_dat, OUT_VALUES, row.names = FALSE, quote = TRUE)
write.csv(
  data.frame(
    Input = c("Targeted pooled seven-KO table", "Prevalence-filtered genome-wide pooled KO table"),
    File = c(INPUT_TARGET, INPUT_GENOME),
    MD5 = unname(tools::md5sum(c(INPUT_TARGET, INPUT_GENOME))),
    stringsAsFactors = FALSE
  ),
  OUT_INPUTS, row.names = FALSE, quote = TRUE
)
write.csv(
  data.frame(
    width_in = WIDTH_IN, height_in = HEIGHT_IN, dpi = DPI,
    width_px = EXPECTED_PX[["width"]], height_px = EXPECTED_PX[["height"]],
    background = "white", colorspace = "sRGB", bit_depth = 8,
    stringsAsFactors = FALSE
  ),
  OUT_RENDER, row.names = FALSE, quote = TRUE
)

save_nc(p, OUT_PNG, width_in = WIDTH_IN, height_in = HEIGHT_IN, dpi = DPI)
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
    stop("Embedding sRGB/dpi metadata changed decoded PNG pixels.")
  }
}

standardize_png(OUT_PNG, DPI)
if (!file.copy(OUT_PNG, TRANSPORT_PNG, overwrite = TRUE)) {
  stop("Could not create the PowerPoint transport PNG.")
}

for (path in c(OUT_PNG, OUT_PDF, OUT_VALUES, OUT_INPUTS, OUT_RENDER,
               OUT_SESSION, OUT_README, TRANSPORT_PNG)) {
  if (file.exists(path)) {
    invisible(system2("/usr/bin/xattr", args = c("-c", shQuote(path)), stdout = TRUE, stderr = TRUE))
    invisible(system2("/bin/chmod", args = c("644", shQuote(path)), stdout = TRUE, stderr = TRUE))
    invisible(system2("/usr/bin/chflags", args = c("nohidden", shQuote(path)), stdout = TRUE, stderr = TRUE))
  }
}
if (!identical(unname(tools::md5sum(OUT_PNG)), unname(tools::md5sum(TRANSPORT_PNG)))) {
  stop("PowerPoint transport PNG is not byte-identical to the analytical PNG.")
}

writeLines(
  c(
    "# CRC WGS pooled targeted sarcosine 7-KO lollipop",
    "",
    "- Four CRC WGS cohorts were pooled in the frozen upstream analysis.",
    "- All seven detected predefined sarcosine-metabolism KOs are displayed.",
    "- log2FC is Cancer/Healthy and was reconstructed from the frozen group means with pseudocount 1e-8, matching the pooled KO pipeline.",
    "- Wilcoxon P and BH q values are taken from `sarcosine_KO_comparison_filtered_pooled.csv`.",
    "- BH correction is across the seven targeted KOs, not across the genome-wide 7,106-KO family.",
    "- Circle area encodes overall prevalence; the legend reports 10%, 50% and 90% reference sizes.",
    "- Filled circles passed the saved >=10% overall-prevalence filter; open circles did not.",
    "- The direction legend reports Higher in Healthy, Higher in Cancer and Not significant.",
    "- Low-prevalence KOs are displayed because the author requested all seven; they require cautious interpretation.",
    "- The older three-KO figure remains unchanged and used genome-wide q values only for prevalence-qualified KOs.",
    sprintf("- Final analytical PNG: `%s`", OUT_PNG),
    sprintf("- PowerPoint transport PNG: `%s`", TRANSPORT_PNG),
    sprintf("- Render size: %.2f x %.2f inches at %d dpi (%d x %d px).", WIDTH_IN, HEIGHT_IN, DPI, EXPECTED_PX[["width"]], EXPECTED_PX[["height"]])
  ),
  OUT_README
)
capture.output(sessionInfo(), file = OUT_SESSION)

manifest_files <- c(OUT_PNG, OUT_PDF, OUT_VALUES, OUT_INPUTS, OUT_RENDER, OUT_README, OUT_SESSION)
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
