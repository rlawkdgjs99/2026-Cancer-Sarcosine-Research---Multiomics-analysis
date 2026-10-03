# =============================================================================
# Figure 1j as a two-group boxplot at the RUN level — display change only
# Date: 2026-08-25
#
# AUTHOR DECISION: keep Figure 1j on the same observational unit as the rest of the
# Figure 1 h-o block, which the manuscript legend declares as "KO n = 1,647
# (745/902)". The BioSample-level alternative (n = 560) was built first and is kept
# at results_integrated/Fig1j_boxplot_BioSample_26.08.25/ as the sensitivity; the
# author chose run level for block consistency after being shown both.
#
# NO STATISTIC IS RECOMPUTED FOR DISPLAY. The panel prints the same pooled Wilcoxon
# P and genome-wide BH q that the current lollipop prints and that the manuscript
# already states, so this change requires no edit to the Results, Methods or legend
# statistics. The script recomputes the Wilcoxon P from the per-run table purely as a
# gate and stops unless it reproduces the stored values exactly.
#
# STILL OPEN, and deliberately not fixed here: 1,647 runs come from 560 BioSamples
# (PRJEB6070 1,066/157; PRJEB27928 260/82), so a run-level boxplot draws repeated
# sequencing of the same biological sample as separate points. Section 9 Priority 1
# still requires the whole CRC-WGS block to move to BioSample units before
# submission, at which point this panel must be rebuilt from the sensitivity version.
#
# Run: LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript 107_Fig1j_boxplot_runlevel_26.08.25.R
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

per_run <- file.path(proj, "results_integrated", "sarcosine",
                     "sarcosine_KO_per_sample_pooled.csv")
frozen <- file.path(proj, "results_integrated",
                    "publication_style_plots_PROPOSED_FIG1H_KO_26.08.19",
                    "Proposed_Fig1h_targeted_sarcosine_KOs_values.csv")
out_dir <- file.path(proj, "results_integrated", "Fig1j_boxplot_runlevel_26.08.25")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(per_run), file.exists(frozen))

GENES <- c("K00303", "K08688")
LABEL <- c(K00303 = "soxB (K00303)", K08688 = "Creatinase (K08688)")
ROLE  <- c(K00303 = "Degradation",   K08688 = "Production")

dat <- read.csv(per_run, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(dat) == 1647L, !anyDuplicated(dat$Run.ID))
stopifnot(identical(as.integer(table(dat$Group)[c("Cancer", "Healthy")]), c(902L, 745L)))

fz <- read.csv(frozen, stringsAsFactors = FALSE)
fz <- fz[match(GENES, fz$KO), ]
stopifnot(nrow(fz) == 2L)

# ---- reproduction gate: the displayed P/q must be the already published ones ----
recomputed <- vapply(GENES, function(ko) {
  h <- dat[[ko]][dat$Group == "Healthy"]; k <- dat[[ko]][dat$Group == "Cancer"]
  suppressWarnings(wilcox.test(k, h, exact = FALSE))$p.value
}, numeric(1))
if (max(abs(recomputed - fz$wilcox_p) / fz$wilcox_p) > 1e-9) {
  stop("Run-level Wilcoxon P does not reproduce the stored panel values.")
}
message("   reproduction gate: PASS (Wilcoxon P matches the published panel values)")

# per-group descriptive numbers for the panel and its source data
desc <- do.call(rbind, lapply(GENES, function(ko) {
  h <- dat[[ko]][dat$Group == "Healthy"]; k <- dat[[ko]][dat$Group == "Cancer"]
  data.frame(KO = ko, gene = LABEL[[ko]], role = ROLE[[ko]],
             n_healthy = length(h), n_cancer = length(k),
             median_healthy = median(h), median_cancer = median(k),
             mean_healthy = mean(h), mean_cancer = mean(k),
             prevalence_healthy = mean(h > 0), prevalence_cancer = mean(k > 0),
             prevalence_overall = mean(c(h, k) > 0),
             log2FC_cancer_vs_healthy = fz$log2FC_Cancer_vs_Healthy[fz$KO == ko],
             wilcox_p = fz$wilcox_p[fz$KO == ko],
             genomewide_BH_q = fz$genomewide_BH_q[fz$KO == ko],
             stringsAsFactors = FALSE)
}))
write.csv(desc, file.path(out_dir, "Fig1j_runlevel_summary.csv"), row.names = FALSE)
write.csv(dat[, c("Run.ID", "Cohort", "Group", GENES)],
          file.path(out_dir, "Fig1j_runlevel_source_data.csv"), row.names = FALSE)

# ---- panel -------------------------------------------------------------------
# 7.2% of soxB and 44.5% of creatinase runs are exactly zero. Those zeros are drawn
# at a disclosed floor rather than cropped away: cropping to a percentile window
# would push nearly half of the creatinase observations off the panel.
#
# The two genes differ by roughly two orders of magnitude in absolute abundance and
# their positive values span 3.3 and 4.6 decades, so a shared axis would force both
# facets onto ~5.5 decades and flatten a 1.5-fold shift into invisibility. Each facet
# therefore carries its own axis, with its own floor placed just under that gene's
# smallest positive value. Axis ticks are labelled per facet so the change of scale
# between the two panels is explicit.
FLOOR <- vapply(GENES, function(ko) {
  v <- dat[[ko]]; 10^(log10(min(v[v > 0])) - 0.4)
}, numeric(1))

FLEV <- sprintf("%s\n%s", LABEL[GENES], ROLE[GENES])
fac <- function(ko) factor(sprintf("%s\n%s", LABEL[ko], ROLE[ko]), levels = FLEV)
long <- dat |>
  tidyr::pivot_longer(all_of(GENES), names_to = "KO", values_to = "abund") |>
  mutate(y = pmax(abund, FLOOR[KO]),
         Group = factor(Group, levels = c("Healthy", "Cancer")),
         facet = fac(KO))

sci10 <- function(x) {                       # 1e-03 -> 10^-3, Nature tick style
  ifelse(is.na(x), NA, sprintf("10^%d", round(log10(x))))
}
q_lab <- function(q) {
  e <- floor(log10(q)); m <- q / 10^e
  if (q < 0.001) sprintf("italic(q) == %.2f %%*%% 10^{%d}", m, e)
  else sprintf("italic(q) == %.4f", q)
}

rng <- long |> group_by(KO, facet) |>
  summarise(lo = min(y), hi = max(y), .groups = "drop")
ann <- desc |> mutate(facet = fac(KO), lab = vapply(genomewide_BH_q, q_lab, character(1))) |>
  left_join(rng, by = c("KO", "facet"))
prev <- desc |>
  select(KO, prevalence_healthy, prevalence_cancer) |>
  tidyr::pivot_longer(starts_with("prevalence"), names_to = "g", values_to = "prev") |>
  mutate(Group = factor(ifelse(grepl("healthy", g), "Healthy", "Cancer"),
                        levels = c("Healthy", "Cancer")),
         facet = fac(KO), lab = sprintf("%.0f%%", prev * 100),
         y = FLOOR[KO] * 0.45)
# headroom for the q annotation and clearance for the prevalence labels
hdr <- rng |> mutate(lo = lo * 0.24, hi = hi * 9) |>
  tidyr::pivot_longer(c(lo, hi), values_to = "y") |> mutate(Group = "Healthy")

p <- ggplot(long, aes(Group, y, fill = Group)) +
  geom_hline(data = rng, aes(yintercept = lo), linetype = "22",
             linewidth = NC_AXIS_LW, colour = "grey70", inherit.aes = FALSE) +
  box_points_layers(width = 0.56, pt_size = 0.20, pt_alpha = 0.14,
                    jitter_w = 0.14, seed = 250825) +
  geom_blank(data = hdr, aes(x = Group, y = y), inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = 1.5, y = hi * 3.4, label = lab), parse = TRUE,
            inherit.aes = FALSE, family = "Arial", size = mm_text(NC_ANNOT_PT - 0.6)) +
  geom_text(data = prev, aes(x = Group, y = y, label = lab),
            inherit.aes = FALSE, family = "Arial",
            size = mm_text(NC_ANNOT_PT - 1.4), colour = "grey35") +
  facet_wrap(~facet, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_y_log10(labels = function(x) parse(text = sci10(x)),
                expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "Relative abundance") +
  theme_nc() +
  theme(plot.margin = margin(1, 2, 1, 1),
        panel.spacing.x = unit(4, "pt"),
        strip.text = element_text(size = NC_STRIP_PT - 1, lineheight = 0.92))

W_IN <- 2.034; H_IN <- 1.320       # measured Figure 1j slot; no PowerPoint scaling
png_path <- file.path(out_dir, "Fig1j_boxplot_runlevel.png")
tmp <- tempfile(fileext = ".png")
ggsave(tmp, p, width = W_IN, height = H_IN, dpi = 600, units = "in", bg = "white",
       device = ragg::agg_png)
stopifnot(file.copy(tmp, png_path, overwrite = TRUE)); unlink(tmp)
dm <- dim(png::readPNG(png_path))
stopifnot(dm[2] == round(W_IN * 600), dm[1] == round(H_IN * 600))
tmp2 <- tempfile(fileext = ".pdf")
ggsave(tmp2, p, width = W_IN, height = H_IN, units = "in", bg = "white", device = cairo_pdf)
file.copy(tmp2, file.path(out_dir, "Fig1j_boxplot_runlevel.pdf"), overwrite = TRUE)
unlink(tmp2)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

outs <- setdiff(list.files(out_dir, full.names = TRUE),
                file.path(out_dir, "output_sha256.tsv"))
write.table(data.frame(file = basename(outs),
  sha256 = vapply(outs, function(f) strsplit(system2("shasum", c("-a", "256", shQuote(f)),
    stdout = TRUE), "[[:space:]]+")[[1]][1], character(1)), row.names = NULL),
  file.path(out_dir, "output_sha256.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=============== Figure 1j, run level (n = 1,647) ===============\n")
print(desc[, c("gene", "role", "n_healthy", "n_cancer", "median_healthy",
               "median_cancer", "prevalence_healthy", "prevalence_cancer",
               "log2FC_cancer_vs_healthy", "wilcox_p", "genomewide_BH_q")],
      row.names = FALSE, digits = 3)
cat(sprintf("\ndisplay floor = %.2e ; panel %.3f x %.3f in (%d x %d px)\n",
            FLOOR, W_IN, H_IN, dm[2], dm[1]))
cat("Outputs:", out_dir, "\n")
