# =============================================================================
# Figure 1j display variants — DISPLAY ONLY, no statistic is recomputed
# Date: 2026-08-25
#
# Reads the frozen BioSample-level source data and result table written by
# 105_Fig1j_boxplot_BioSample_26.08.25.R and redraws the same numbers three ways so
# the author can choose. Every variant shows the identical Wilcoxon P values.
#
#   A  boxplot, full axis            (as first produced)
#   B  boxplot, axis cropped to the 1st-99th percentile of positive values
#   C  empirical cumulative distribution
#
# Why B and C exist: the full axis has to span 3.3 decades for soxB and 4.5 for
# creatinase plus a zero floor, in a 1.32 in tall panel, in order to display a shift
# of roughly 0.2 log units. Outliers therefore set the scale and visually flatten the
# very difference the panel is about. B crops the DISPLAY only — the test still uses
# every observation, and the cropped fraction is printed on the panel. C avoids the
# problem entirely: the gap between two ECDF curves is the shift at every abundance.
#
# Run: LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript 106_Fig1j_display_variants_26.08.25.R
# =============================================================================

set.seed(42)
options(warn = 1)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(ragg) })
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"), fixed = TRUE)) stop("Need a UTF-8 locale.")

script_file <- sub("^--file=", "",
                   grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
proj <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
root <- normalizePath(file.path(proj, "..", "..", ".."), mustWork = TRUE)
source(file.path(root, "공공_Metabolomics&Metagenomics_분석모음", "_shared",
                 "theme_nc_26.08.18.R"))
dir_in <- file.path(proj, "results_integrated", "Fig1j_boxplot_BioSample_26.08.25")
bios <- read.csv(file.path(dir_in, "Fig1j_source_data_BioSample.csv"),
                 stringsAsFactors = FALSE)
prim <- read.csv(file.path(dir_in, "Fig1j_primary_BioSample.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(bios) == 560L, nrow(prim) == 2L)

GENES <- c("K00303", "K08688")
LABEL <- c(K00303 = "soxB (K00303)", K08688 = "Creatinase (K08688)")
ROLE  <- c(K00303 = "Degradation",   K08688 = "Production")
FLEV  <- sprintf("%s\n%s", LABEL[GENES], ROLE[GENES])
W_IN <- 2.034; H_IN <- 1.320

long <- bios |>
  tidyr::pivot_longer(all_of(GENES), names_to = "KO", values_to = "abund") |>
  mutate(Group = factor(Group, levels = c("Healthy", "Cancer")),
         facet = factor(sprintf("%s\n%s", LABEL[KO], ROLE[KO]), levels = FLEV))

fmt_p <- function(p) if (p < 0.001) {
  e <- floor(log10(p)); sprintf("italic(P) == %.1f %%*%% 10^{%d}", p / 10^e, e)
} else sprintf("italic(P) == %.3f", p)
ann <- prim |>
  mutate(facet = factor(sprintf("%s\n%s", gene, role), levels = FLEV),
         lab = vapply(p_value, fmt_p, character(1)))

save_variant <- function(pl, stem) {
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, pl, width = W_IN, height = H_IN, dpi = 600, units = "in",
         bg = "white", device = ragg::agg_png)
  out <- file.path(dir_in, paste0(stem, ".png"))
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp)
  dm <- dim(png::readPNG(out))
  stopifnot(dm[2] == round(W_IN * 600), dm[1] == round(H_IN * 600))
  tmp2 <- tempfile(fileext = ".pdf")
  ggsave(tmp2, pl, width = W_IN, height = H_IN, units = "in", bg = "white",
         device = cairo_pdf)
  file.copy(tmp2, file.path(dir_in, paste0(stem, ".pdf")), overwrite = TRUE); unlink(tmp2)
  cat(sprintf("   %-34s %d x %d px\n", stem, dm[2], dm[1]))
}

# ---- B: cropped boxplot ------------------------------------------------------
# Crop bounds come from positive values only, per gene, and the proportion of
# observations falling outside the drawn window is printed so nothing is hidden.
lims <- long |> filter(abund > 0) |> group_by(facet) |>
  summarise(lo = quantile(abund, 0.01), hi = quantile(abund, 0.99), .groups = "drop")
LO <- min(lims$lo); HI <- max(lims$hi)
outside <- mean(long$abund < LO | long$abund > HI)
cat(sprintf("variant B: display window %.2e - %.2e, %.1f%% of observations outside\n",
            LO, HI, outside * 100))

pB <- ggplot(long, aes(Group, abund, fill = Group)) +
  box_points_layers(width = 0.56, pt_size = 0.30, pt_alpha = 0.22,
                    jitter_w = 0.13, seed = 250825) +
  geom_text(data = ann, aes(x = 1.5, y = HI * 1.9, label = lab), parse = TRUE,
            inherit.aes = FALSE, family = "Arial", size = mm_text(NC_ANNOT_PT - 0.6)) +
  facet_wrap(~facet, nrow = 1) +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_y_log10(labels = function(x) format(x, scientific = TRUE, digits = 1)) +
  coord_cartesian(ylim = c(LO, HI * 3.2), clip = "off") +
  labs(x = NULL, y = "Relative abundance") +
  theme_nc() +
  theme(plot.margin = margin(1, 2, 1, 1),
        strip.text = element_text(size = NC_STRIP_PT - 1, lineheight = 0.92))
save_variant(pB, "Fig1j_variantB_boxplot_cropped")

# ---- C: empirical cumulative distribution ------------------------------------
ecdf_dat <- long |> filter(abund > 0) |> arrange(facet, Group, abund) |>
  group_by(facet, Group) |> mutate(F = seq_len(dplyr::n()) / dplyr::n()) |> ungroup()

pC <- ggplot(ecdf_dat, aes(abund, F, colour = Group)) +
  geom_step(linewidth = pt_lw(1.1)) +
  geom_text(data = ann, aes(x = HI, y = 0.06, label = lab), parse = TRUE, hjust = 1,
            inherit.aes = FALSE, family = "Arial", size = mm_text(NC_ANNOT_PT - 0.6)) +
  facet_wrap(~facet, nrow = 1) +
  scale_colour_manual(values = GROUP_COLORS, name = NULL) +
  scale_x_log10(breaks = c(1e-5, 1e-3, 1e-1),
                labels = c(expression(10^-5), expression(10^-3), expression(10^-1))) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = expansion(mult = 0.03)) +
  coord_cartesian(xlim = c(LO, HI), clip = "off") +
  labs(x = "Relative abundance", y = "Cumulative fraction") +
  theme_nc() +
  theme(legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, 1, 0), legend.box.spacing = unit(1, "pt"),
        legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = NC_LEGEND_PT - 0.5),
        plot.margin = margin(1, 3, 1, 1),
        strip.text = element_text(size = NC_STRIP_PT - 1, lineheight = 0.92))
save_variant(pC, "Fig1j_variantC_ecdf")

write.csv(data.frame(
  variant = c("A", "B", "C"),
  description = c("boxplot, full axis with disclosed zero floor",
                  "boxplot, display cropped to the 1st-99th percentile of positive values",
                  "empirical cumulative distribution of positive values"),
  display_window_low = c(NA, LO, LO), display_window_high = c(NA, HI, HI),
  pct_observations_outside_window = c(0, outside * 100, outside * 100),
  note = c("zeros shown at a disclosed floor; prevalence annotated",
           "test uses every observation; cropping is display-only",
           "zeros omitted from the curve; prevalence differs by 0.9-8.4% between genes"),
  stringsAsFactors = FALSE),
  file.path(dir_in, "Fig1j_display_variants.csv"), row.names = FALSE)
cat("variants written to:", dir_in, "\n")
