# =============================================================================
# Figure 1j redesigned as a two-group boxplot at BioSample resolution
# Date: 2026-08-25
# Frozen plan: results_integrated/Fig1j_boxplot_BioSample_26.08.25/ANALYSIS_PLAN_FROZEN.md
#              SHA-256 8ed6c49f69b71d20483d45cc26f2c8137e8f89a08d4421b51e0ee5ca4f70fea2
#
# WHY THE UNIT CHANGES
#   sarcosine_KO_per_sample_pooled.csv is keyed by Run.ID: 1,647 rows from only 560
#   BioSamples (PRJEB6070 1,066 runs / 157 samples; PRJEB27928 260 / 82). A lollipop
#   hid that inside a P value; a boxplot draws every point, so the run level would
#   display ~7 repeats of the same patient as independent observations. Runs are
#   therefore averaged within ENA sample_accession before any test, and the run level
#   is retained only as an explicit sensitivity.
#
# Genes: the two BH-significant sarcosine KOs in the Figure 1i genome-wide volcano,
#   soxB (K00303, degradation) and creatinase (K08688, production). soxA (K00302,
#   q = 0.735) is excluded by author instruction.
#
# Run: LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript 105_Fig1j_boxplot_BioSample_26.08.25.R
# =============================================================================

set.seed(42)
options(warn = 1)
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(ragg)
})
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"), fixed = TRUE)) {
  stop("Refusing to run outside a UTF-8 locale.")
}

script_file <- sub("^--file=", "",
                   grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
stopifnot(length(script_file) == 1L)
proj <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
root <- normalizePath(file.path(proj, "..", "..", ".."), mustWork = TRUE)

per_sample <- file.path(proj, "results_integrated", "sarcosine",
                        "sarcosine_KO_per_sample_pooled.csv")
map_dir <- file.path(root, "공공_Metabolomics&Metagenomics_분석모음",
                     "HGMT_NSCLC_ICI_RvsNR_WGS", "pooled_analysis", "resources",
                     "ENA_CRC_run_sample_mapping_26.08.23")
theme_file <- file.path(root, "공공_Metabolomics&Metagenomics_분석모음", "_shared",
                        "theme_nc_26.08.18.R")
hist_cmp <- file.path(proj, "results_integrated", "sarcosine",
                      "sarcosine_KO_comparison_pooled.csv")
out_dir <- file.path(proj, "results_integrated", "Fig1j_boxplot_BioSample_26.08.25")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(per_sample), dir.exists(map_dir), file.exists(theme_file),
          file.exists(hist_cmp))
source(theme_file)

GENES <- c("K00303", "K08688")
LABEL <- c(K00303 = "soxB (K00303)", K08688 = "Creatinase (K08688)")
ROLE  <- c(K00303 = "Degradation",   K08688 = "Production")

# ---- inputs ------------------------------------------------------------------
dat <- read.csv(per_sample, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(dat) == 1647L, !anyDuplicated(dat$Run.ID),
          all(GENES %in% names(dat)), setequal(dat$Group, c("Healthy", "Cancer")))

maps <- do.call(rbind, lapply(list.files(map_dir, pattern = "_ena_runs[.]tsv$",
                                         full.names = TRUE), function(f) {
  m <- read.delim(f, stringsAsFactors = FALSE)
  stopifnot(all(c("run_accession", "sample_accession") %in% names(m)))
  m[, c("run_accession", "sample_accession")]
}))
stopifnot(!anyDuplicated(maps$run_accession))
dat$BioSample <- maps$sample_accession[match(dat$Run.ID, maps$run_accession)]
if (anyNA(dat$BioSample)) stop("Unmapped runs: ", sum(is.na(dat$BioSample)))

# A BioSample must not straddle Healthy/Cancer or two cohorts.
chk <- dat |> group_by(BioSample) |>
  summarise(g = n_distinct(Group), c = n_distinct(Cohort), .groups = "drop")
if (any(chk$g > 1) || any(chk$c > 1)) stop("BioSample spans multiple groups/cohorts.")

unit_tab <- dat |> group_by(Cohort) |>
  summarise(runs = n(), biosamples = n_distinct(BioSample), .groups = "drop")
message("== run -> BioSample ==")
print(as.data.frame(unit_tab), row.names = FALSE)
stopifnot(sum(unit_tab$biosamples) == 560L)

# ---- BioSample aggregation (primary) ----------------------------------------
bios <- dat |>
  group_by(BioSample, Cohort, Group) |>
  summarise(across(all_of(GENES), mean), n_runs = n(), .groups = "drop")
stopifnot(nrow(bios) == 560L)

# ---- testing -----------------------------------------------------------------
rb_indep <- function(a, b) {                 # oriented Cancer - Healthy
  2 * (mean(outer(a, b, ">")) + 0.5 * mean(outer(a, b, "=="))) - 1
}
boot_rb <- function(a, b, n = 10000L) {
  v <- vapply(seq_len(n), function(i)
    rb_indep(sample(a, length(a), TRUE), sample(b, length(b), TRUE)), numeric(1))
  unname(quantile(v[is.finite(v)], c(0.025, 0.975)))
}

test_one <- function(df, ko, unit) {
  h <- df[[ko]][df$Group == "Healthy"]; k <- df[[ko]][df$Group == "Cancer"]
  w <- suppressWarnings(wilcox.test(k, h, alternative = "two.sided", exact = FALSE))
  rb <- rb_indep(k, h); ci <- boot_rb(k, h)
  data.frame(
    unit = unit, KO = ko, gene = LABEL[[ko]], role = ROLE[[ko]],
    n_healthy = length(h), n_cancer = length(k),
    median_healthy = median(h), median_cancer = median(k),
    mean_healthy = mean(h), mean_cancer = mean(k),
    prevalence_healthy = mean(h > 0), prevalence_cancer = mean(k > 0),
    wilcoxon_W = unname(w$statistic), p_value = w$p.value,
    rank_biserial_cancer_higher = rb, rb_ci_low = ci[1], rb_ci_high = ci[2],
    stringsAsFactors = FALSE)
}

primary <- do.call(rbind, lapply(GENES, function(g) test_one(bios, g, "BioSample")))
sens_run <- do.call(rbind, lapply(GENES, function(g) test_one(dat, g, "Run (historical)")))

# cohort-stratified sensitivity at BioSample level
strat <- do.call(rbind, lapply(GENES, function(ko) {
  per <- lapply(split(bios, bios$Cohort), function(s) {
    h <- s[[ko]][s$Group == "Healthy"]; k <- s[[ko]][s$Group == "Cancer"]
    if (length(h) < 5 || length(k) < 5) return(NULL)
    w <- suppressWarnings(wilcox.test(k, h, exact = FALSE))
    data.frame(KO = ko, Cohort = s$Cohort[1], n_h = length(h), n_c = length(k),
               rb = rb_indep(k, h), p = w$p.value, stringsAsFactors = FALSE)
  })
  do.call(rbind, per)
}))

hist <- read.csv(hist_cmp, stringsAsFactors = FALSE)
hist <- hist[hist$KO %in% GENES, c("KO", "p_value", "p_adj")]
names(hist) <- c("KO", "historical_run_p", "historical_run_BH_q")
primary <- merge(primary, hist, by = "KO", sort = FALSE)
primary <- primary[match(GENES, primary$KO), ]

write.csv(primary, file.path(out_dir, "Fig1j_primary_BioSample.csv"), row.names = FALSE)
write.csv(sens_run, file.path(out_dir, "Fig1j_sensitivity_run_level.csv"), row.names = FALSE)
write.csv(strat, file.path(out_dir, "Fig1j_sensitivity_cohort_stratified.csv"), row.names = FALSE)
write.csv(as.data.frame(unit_tab), file.path(out_dir, "Fig1j_run_to_biosample_audit.csv"), row.names = FALSE)
write.csv(bios, file.path(out_dir, "Fig1j_source_data_BioSample.csv"), row.names = FALSE)

# ---- panel -------------------------------------------------------------------
# Zeros are 7% of soxB and 45% of creatinase, so they must stay visible. A disclosed
# display floor half a decade below the smallest positive value keeps them on the
# log scale instead of silently dropping them.
# Place the floor just under the smallest positive value rather than under the
# decade boundary; the latter left an empty half-decade of dead space in a 1.3 in panel.
pos_min <- min(unlist(bios[GENES])[unlist(bios[GENES]) > 0])
FLOOR <- 10^(log10(pos_min) - 0.35)

plot_dat <- bios |>
  tidyr::pivot_longer(all_of(GENES), names_to = "KO", values_to = "abund") |>
  mutate(
    y = pmax(abund, FLOOR),
    is_zero = abund == 0,
    Group = factor(Group, levels = c("Healthy", "Cancer")),
    facet = factor(sprintf("%s\n%s", LABEL[KO], ROLE[KO]),
                   levels = sprintf("%s\n%s", LABEL[GENES], ROLE[GENES])))

fmt_p <- function(p) if (p < 0.001) {
  e <- floor(log10(p)); sprintf("italic(P) == %.1f %%*%% 10^{%d}", p / 10^e, e)
} else sprintf("italic(P) == %.3f", p)

ann <- primary |>
  mutate(facet = factor(sprintf("%s\n%s", gene, role),
                        levels = levels(plot_dat$facet)),
         lab = vapply(p_value, fmt_p, character(1)))

prev <- primary |>
  select(gene, role, prevalence_healthy, prevalence_cancer) |>
  tidyr::pivot_longer(starts_with("prevalence"), names_to = "g", values_to = "prev") |>
  mutate(Group = factor(ifelse(grepl("healthy", g), "Healthy", "Cancer"),
                        levels = c("Healthy", "Cancer")),
         facet = factor(sprintf("%s\n%s", gene, role), levels = levels(plot_dat$facet)),
         lab = sprintf("%.0f%%", prev * 100))

ytop <- max(plot_dat$y)
p <- ggplot(plot_dat, aes(Group, y, fill = Group)) +
  geom_hline(yintercept = FLOOR, linetype = "22", linewidth = NC_AXIS_LW,
             colour = "grey70") +
  box_points_layers(width = 0.56, pt_size = 0.28, pt_alpha = 0.22,
                    jitter_w = 0.13, seed = 250825) +
  geom_text(data = ann, aes(x = 1.5, y = ytop * 4.5, label = lab), parse = TRUE,
            inherit.aes = FALSE, family = "Arial", size = mm_text(NC_ANNOT_PT - 0.6)) +
  geom_text(data = prev, aes(x = Group, y = FLOOR * 0.55, label = lab),
            inherit.aes = FALSE, family = "Arial",
            size = mm_text(NC_ANNOT_PT - 1.4), colour = "grey35") +
  facet_wrap(~facet, nrow = 1) +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_y_log10(labels = function(x) format(x, scientific = TRUE, digits = 1),
                expand = expansion(mult = c(0.14, 0.16))) +
  labs(x = NULL, y = "Relative abundance") +
  theme_nc() +
  theme(plot.margin = margin(1, 2, 1, 1),
        strip.text = element_text(size = NC_STRIP_PT - 1, lineheight = 0.92))

W_IN <- 2.034; H_IN <- 1.320
save_png_safe <- function(pl, path, w, h) {
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, pl, width = w, height = h, dpi = 600, units = "in", bg = "white",
         device = ragg::agg_png)
  stopifnot(file.copy(tmp, path, overwrite = TRUE)); unlink(tmp)
  dm <- dim(png::readPNG(path))
  if (dm[2] != round(w * 600) || dm[1] != round(h * 600))
    stop(basename(path), " is ", dm[2], "x", dm[1])
  invisible(path)
}
png_path <- file.path(out_dir, "Fig1j_boxplot_BioSample.png")
save_png_safe(p, png_path, W_IN, H_IN)
tmp <- tempfile(fileext = ".pdf")
ggsave(tmp, p, width = W_IN, height = H_IN, units = "in", bg = "white", device = cairo_pdf)
file.copy(tmp, file.path(out_dir, "Fig1j_boxplot_BioSample.pdf"), overwrite = TRUE)
unlink(tmp)

outs <- setdiff(list.files(out_dir, full.names = TRUE),
                file.path(out_dir, "output_sha256.tsv"))
write.table(data.frame(file = basename(outs),
  sha256 = vapply(outs, function(f) strsplit(system2("shasum", c("-a", "256", shQuote(f)),
    stdout = TRUE), "[[:space:]]+")[[1]][1], character(1)), row.names = NULL),
  file.path(out_dir, "output_sha256.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("\n=============== PRIMARY (BioSample, n = 560) ===============\n")
print(primary[, c("gene", "role", "n_healthy", "n_cancer", "median_healthy",
                  "median_cancer", "prevalence_healthy", "prevalence_cancer",
                  "p_value", "rank_biserial_cancer_higher", "historical_run_BH_q")],
      row.names = FALSE, digits = 3)
cat("\n=============== SENSITIVITY (run level, n = 1,647) ===============\n")
print(sens_run[, c("gene", "n_healthy", "n_cancer", "p_value",
                   "rank_biserial_cancer_higher")], row.names = FALSE, digits = 3)
cat("\n=============== SENSITIVITY (cohort-stratified, BioSample) ===============\n")
print(strat, row.names = FALSE, digits = 3)
cat(sprintf("\ndisplay floor = %.2e ; panel %.3f x %.3f in\n", FLOOR, W_IN, H_IN))
cat("Outputs:", out_dir, "\n")
