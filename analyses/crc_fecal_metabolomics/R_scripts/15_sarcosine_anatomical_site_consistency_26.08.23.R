#!/usr/bin/env Rscript

# Proposed Figure 1c: fecal sarcosine across CRC anatomical sites.
#
# Biological question
#   Is the CRC-associated increase in the MTBLS10232 NEG-mode feature
#   putatively annotated as sarcosine confined to one anatomical site, or is
#   the direction reproduced in right-colon, left-colon, and rectal cancers?
#
# Statistical design
#   - Unit: one deposited patient fecal sample (unique sample identifiers).
#   - Design: independent Control, right-colon, left-colon, and rectal groups.
#   - Outcome: log2 quantile-normalized NEG-mode feature intensity. Quantile
#     normalization is applied across samples with features in rows, matching
#     the corrected Figure 1b pipeline.
#   - Covariate: deposited younger/older age stratum. Individual ages and other
#     clinical covariates are unavailable in the deposited metadata.
#   - Primary contrasts: each specified CRC site versus Control from one linear
#     model. HC3 heteroscedasticity-robust 95% CIs and P values are reported;
#     the three planned P values are Holm-adjusted.
#   - Consistency test: a two-df HC3 Wald test asks whether the adjusted means
#     differ among the three specified CRC sites. Failure to reject this test
#     is not proof of equality and is interpreted only with effect estimates.
#   - Sensitivity checks: 5,000 cell-stratified bootstrap fits, age-stratified
#     site-versus-Control models, and a site-by-age interaction test.
#   - Exclusions: 62 QC injections and six CRC samples with unspecified site.
#
# Structural-identification caveat
#   The public MAF annotation is putative (HMDB0000271; m/z 88.03885, RT
#   7.4811 min), lacks deposited MS/MS confirmation, and is isobaric with
#   structural isomers. The analysis therefore concerns this deposited feature
#   annotation and does not independently confirm molecular identity.

suppressPackageStartupMessages({
  library(digest)
  library(ggplot2)
  library(here)
  library(limma)
  library(ragg)
  library(svglite)
})

set.seed(
  42L,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path.")
SCRIPT_PATH <- normalizePath(sub("^--file=", "", script_arg))
BASE_DIR <- normalizePath(file.path(dirname(SCRIPT_PATH), ".."))
setwd(BASE_DIR)
here::i_am("R_scripts/15_sarcosine_anatomical_site_consistency_26.08.23.R")

SHARED_THEME <- normalizePath(file.path(
  here::here(), "..", "_shared", "theme_nc_26.08.18.R"
))
source(SHARED_THEME)

INPUT_RDS <- here::here("analysis", "intensity_matrices.rds")
INPUT_METADATA <- here::here("metadata_merged_619samples.tsv")
OUT_DIR <- here::here("Manuscript_Final_Panels_26.08.23")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_BASE <- file.path(
  OUT_DIR,
  "Proposed_Fig1c_sarcosine_anatomical_site_consistency"
)
OUT_PNG <- paste0(OUT_BASE, ".png")
OUT_PDF <- paste0(OUT_BASE, ".pdf")
OUT_SVG <- paste0(OUT_BASE, ".svg")
OUT_PRIMARY <- paste0(OUT_BASE, "_site_vs_control_statistics.csv")
OUT_HETEROGENEITY <- paste0(OUT_BASE, "_site_heterogeneity_statistics.csv")
OUT_PAIRWISE <- paste0(OUT_BASE, "_crc_site_pairwise_statistics.csv")
OUT_AGE_SENSITIVITY <- paste0(OUT_BASE, "_age_stratified_sensitivity.csv")
OUT_SOURCE <- paste0(OUT_BASE, "_source_data.csv")
OUT_CHECKSUM <- paste0(OUT_BASE, "_input_sha256.tsv")
OUT_VALIDATION <- paste0(OUT_BASE, "_validation_report.txt")
OUT_SESSION <- paste0(OUT_BASE, "_sessionInfo.txt")

if (!file.exists(INPUT_RDS)) stop("Missing input: ", INPUT_RDS)
if (!file.exists(INPUT_METADATA)) stop("Missing input: ", INPUT_METADATA)

# HC3 covariance without relying on an additional package.
vcov_hc3 <- function(fit) {
  x <- model.matrix(fit)
  e <- residuals(fit)
  h <- hatvalues(fit)
  if (any(h >= 1)) stop("A leverage value is >= 1; HC3 is undefined.")
  bread <- solve(crossprod(x))
  meat <- crossprod(x, x * as.numeric(e^2 / (1 - h)^2))
  bread %*% meat %*% bread
}

contrast_table <- function(fit, contrast_matrix, labels) {
  if (nrow(contrast_matrix) != length(labels)) stop("Contrast label mismatch.")
  b <- coef(fit)
  v <- vcov_hc3(fit)
  est <- as.numeric(contrast_matrix %*% b)
  se <- sqrt(diag(contrast_matrix %*% v %*% t(contrast_matrix)))
  df <- df.residual(fit)
  crit <- qt(0.975, df = df)
  data.frame(
    Contrast = labels,
    estimate_log2 = est,
    hc3_se = se,
    hc3_ci_95_low = est - crit * se,
    hc3_ci_95_high = est + crit * se,
    hc3_t = est / se,
    hc3_df = df,
    hc3_p = 2 * pt(abs(est / se), df = df, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

joint_wald_hc3 <- function(fit, contrast_matrix) {
  b <- coef(fit)
  v <- vcov_hc3(fit)
  rb <- as.numeric(contrast_matrix %*% b)
  rv <- contrast_matrix %*% v %*% t(contrast_matrix)
  q <- nrow(contrast_matrix)
  f_stat <- as.numeric(t(rb) %*% solve(rv, rb)) / q
  df2 <- df.residual(fit)
  c(F = f_stat, df1 = q, df2 = df2, p = pf(f_stat, q, df2, lower.tail = FALSE))
}

format_p <- function(p) {
  if (p < 0.001) {
    format(p, scientific = TRUE, digits = 2)
  } else {
    sprintf("%.3f", p)
  }
}

# -----------------------------------------------------------------------------
# Data inspection and construction of the corrected QN outcome
# -----------------------------------------------------------------------------
dat <- readRDS(INPUT_RDS)
metadata <- read.delim(
  INPUT_METADATA,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_rds <- c("intensity_neg", "feature_meta_neg")
if (!all(required_rds %in% names(dat))) stop("Required objects missing from RDS.")
required_meta <- c("MTBLS_Sample_Name", "age_cohort", "rcc_lcc")
missing_meta <- setdiff(required_meta, names(metadata))
if (length(missing_meta)) stop("Metadata columns missing: ", paste(missing_meta, collapse = ", "))

intensity <- dat$intensity_neg
feature_meta <- dat$feature_meta_neg
if (!is.matrix(intensity) || storage.mode(intensity) != "double") {
  stop("NEG intensity object is not a double matrix.")
}
if (nrow(intensity) != nrow(feature_meta)) stop("Feature metadata is misaligned.")
if (anyNA(intensity) || any(!is.finite(intensity)) || any(intensity < 0)) {
  stop("NEG matrix contains missing, non-finite, or negative values.")
}
if (anyDuplicated(colnames(intensity))) stop("Duplicated matrix sample IDs.")
if (anyDuplicated(metadata$MTBLS_Sample_Name)) stop("Duplicated metadata sample IDs.")
if (!setequal(colnames(intensity), metadata$MTBLS_Sample_Name)) {
  stop("Matrix and metadata sample-ID sets differ.")
}

positive_values <- intensity[intensity > 0]
if (!length(positive_values)) stop("No positive intensities.")
intensity_imputed <- intensity
intensity_imputed[intensity_imputed == 0] <- min(positive_values) / 2
qn_neg <- limma::normalizeBetweenArrays(log2(intensity_imputed), method = "quantile")
if (!identical(dim(qn_neg), dim(intensity))) stop("Normalization changed dimensions.")

sarc_idx <- which(
  feature_meta$metabolite_identification == "Sarcosine" &
    feature_meta$database_identifier == "HMDB0000271"
)
if (length(sarc_idx) != 1L) stop("Expected exactly one sarcosine/HMDB0000271 feature.")

metadata_aligned <- metadata[
  match(colnames(qn_neg), metadata$MTBLS_Sample_Name),
  ,
  drop = FALSE
]
if (!identical(colnames(qn_neg), metadata_aligned$MTBLS_Sample_Name)) {
  stop("Metadata alignment failed.")
}

analysis_all <- data.frame(
  sample_id = metadata_aligned$MTBLS_Sample_Name,
  sarcosine_log2_qn = as.numeric(qn_neg[sarc_idx, ]),
  site = metadata_aligned$rcc_lcc,
  age_stratum = ifelse(
    grepl("^o", metadata_aligned$age_cohort),
    "Older",
    ifelse(grepl("^y", metadata_aligned$age_cohort), "Younger", NA_character_)
  ),
  stringsAsFactors = FALSE
)

main_sites <- c("Control", "RCC", "LCC", "Rectum")
analysis_data <- analysis_all[
  analysis_all$site %in% main_sites & !is.na(analysis_all$age_stratum),
  ,
  drop = FALSE
]
analysis_data$site <- factor(analysis_data$site, levels = main_sites)
analysis_data$age_stratum <- factor(
  analysis_data$age_stratum,
  levels = c("Younger", "Older")
)

expected_counts <- c(Control = 246L, RCC = 59L, LCC = 87L, Rectum = 159L)
observed_counts <- table(analysis_data$site)
if (!identical(as.integer(observed_counts), as.integer(expected_counts))) {
  stop("Unexpected main-analysis sample counts.")
}
if (any(table(analysis_data$site, analysis_data$age_stratum) == 0L)) {
  stop("A site-by-age cell is empty.")
}
if (anyNA(analysis_data) || any(!is.finite(analysis_data$sarcosine_log2_qn))) {
  stop("Missing or non-finite main-analysis values.")
}

# -----------------------------------------------------------------------------
# Primary model: each anatomical site versus Control, adjusted for age stratum
# -----------------------------------------------------------------------------
fit_primary <- lm(sarcosine_log2_qn ~ site + age_stratum, data = analysis_data)
if (fit_primary$rank != ncol(model.matrix(fit_primary))) stop("Primary design is rank deficient.")

coef_names <- names(coef(fit_primary))
primary_coef <- c("siteRCC", "siteLCC", "siteRectum")
if (!all(primary_coef %in% coef_names)) stop("Unexpected primary-model coefficients.")
primary_l <- matrix(0, nrow = 3L, ncol = length(coef_names), dimnames = list(NULL, coef_names))
primary_l[cbind(seq_len(3L), match(primary_coef, coef_names))] <- 1
primary <- contrast_table(
  fit_primary,
  primary_l,
  c("Right colon vs Control", "Left colon vs Control", "Rectum vs Control")
)
primary$site <- c("RCC", "LCC", "Rectum")
primary$n_site <- as.integer(observed_counts[primary$site])
primary$n_control <- as.integer(observed_counts["Control"])
primary$hc3_p_holm_across_3 <- p.adjust(primary$hc3_p, method = "holm")
primary$fold_change <- 2^primary$estimate_log2
primary$fold_change_ci_95_low <- 2^primary$hc3_ci_95_low
primary$fold_change_ci_95_high <- 2^primary$hc3_ci_95_high

# Cell-stratified bootstrap preserves every site-by-age sample count.
BOOT_REPS <- 5000L
BOOT_SEED <- 42L
set.seed(BOOT_SEED)
strata <- interaction(analysis_data$site, analysis_data$age_stratum, drop = TRUE)
stratum_rows <- split(seq_len(nrow(analysis_data)), strata)
boot_estimates <- matrix(NA_real_, nrow = BOOT_REPS, ncol = 3L)
colnames(boot_estimates) <- primary$site
for (b in seq_len(BOOT_REPS)) {
  sampled_rows <- unlist(
    lapply(stratum_rows, function(i) sample(i, length(i), replace = TRUE)),
    use.names = FALSE
  )
  boot_fit <- lm(
    sarcosine_log2_qn ~ site + age_stratum,
    data = analysis_data[sampled_rows, , drop = FALSE]
  )
  boot_estimates[b, ] <- coef(boot_fit)[primary_coef]
}
if (anyNA(boot_estimates)) stop("Bootstrap produced missing estimates.")
primary$bootstrap_ci_95_low <- apply(boot_estimates, 2, quantile, probs = 0.025, names = FALSE)
primary$bootstrap_ci_95_high <- apply(boot_estimates, 2, quantile, probs = 0.975, names = FALSE)
primary$bootstrap_reps <- BOOT_REPS
primary$bootstrap_seed <- BOOT_SEED

# -----------------------------------------------------------------------------
# Are adjusted means detectably different among the three specified CRC sites?
# -----------------------------------------------------------------------------
crc_data <- droplevels(analysis_data[analysis_data$site != "Control", , drop = FALSE])
crc_data$site <- relevel(crc_data$site, ref = "RCC")
fit_crc_null <- lm(sarcosine_log2_qn ~ age_stratum, data = crc_data)
fit_crc_site <- lm(sarcosine_log2_qn ~ site + age_stratum, data = crc_data)
crc_coef_names <- names(coef(fit_crc_site))
crc_site_coef <- c("siteLCC", "siteRectum")
crc_joint_l <- matrix(
  0,
  nrow = 2L,
  ncol = length(crc_coef_names),
  dimnames = list(NULL, crc_coef_names)
)
crc_joint_l[cbind(seq_len(2L), match(crc_site_coef, crc_coef_names))] <- 1
crc_wald <- joint_wald_hc3(fit_crc_site, crc_joint_l)
crc_classical <- anova(fit_crc_null, fit_crc_site)

heterogeneity <- data.frame(
  Test = c(
    "HC3 joint Wald: RCC = LCC = Rectum",
    "Classical partial F: add site to age-stratum model"
  ),
  F = c(unname(crc_wald["F"]), crc_classical$F[2]),
  df1 = c(unname(crc_wald["df1"]), crc_classical$Df[2]),
  df2 = c(unname(crc_wald["df2"]), crc_classical$Res.Df[2]),
  p = c(unname(crc_wald["p"]), crc_classical$`Pr(>F)`[2]),
  stringsAsFactors = FALSE
)

pairwise_l <- matrix(
  0,
  nrow = 3L,
  ncol = length(crc_coef_names),
  dimnames = list(NULL, crc_coef_names)
)
pairwise_l[1, "siteLCC"] <- 1
pairwise_l[2, "siteRectum"] <- 1
pairwise_l[3, c("siteLCC", "siteRectum")] <- c(-1, 1)
pairwise <- contrast_table(
  fit_crc_site,
  pairwise_l,
  c("Left colon - Right colon", "Rectum - Right colon", "Rectum - Left colon")
)
pairwise$hc3_p_holm_across_3 <- p.adjust(pairwise$hc3_p, method = "holm")

# -----------------------------------------------------------------------------
# Sensitivity: site effects within each deposited age stratum
# -----------------------------------------------------------------------------
age_results <- list()
for (age_level in levels(analysis_data$age_stratum)) {
  z <- droplevels(analysis_data[analysis_data$age_stratum == age_level, , drop = FALSE])
  z$site <- relevel(z$site, ref = "Control")
  fit_age <- lm(sarcosine_log2_qn ~ site, data = z)
  age_coef_names <- names(coef(fit_age))
  age_l <- matrix(0, nrow = 3L, ncol = length(age_coef_names), dimnames = list(NULL, age_coef_names))
  age_l[cbind(seq_len(3L), match(primary_coef, age_coef_names))] <- 1
  age_tab <- contrast_table(
    fit_age,
    age_l,
    c("Right colon vs Control", "Left colon vs Control", "Rectum vs Control")
  )
  age_tab$age_stratum <- age_level
  age_tab$site <- c("RCC", "LCC", "Rectum")
  age_tab$n_site <- as.integer(table(z$site)[age_tab$site])
  age_tab$n_control <- as.integer(table(z$site)["Control"])
  age_tab$hc3_p_holm_within_stratum <- p.adjust(age_tab$hc3_p, method = "holm")
  age_results[[age_level]] <- age_tab
}
age_sensitivity <- do.call(rbind, age_results)
rownames(age_sensitivity) <- NULL

# Site-by-age interaction as a secondary sensitivity test.
fit_interaction <- lm(sarcosine_log2_qn ~ site * age_stratum, data = analysis_data)
interaction_terms <- grep(":", names(coef(fit_interaction)), value = TRUE)
interaction_l <- matrix(
  0,
  nrow = length(interaction_terms),
  ncol = length(coef(fit_interaction)),
  dimnames = list(NULL, names(coef(fit_interaction)))
)
interaction_l[cbind(seq_along(interaction_terms), match(interaction_terms, names(coef(fit_interaction))))] <- 1
interaction_wald <- joint_wald_hc3(fit_interaction, interaction_l)
interaction_classical <- anova(fit_primary, fit_interaction)

# -----------------------------------------------------------------------------
# Descriptive/source tables and exports
# -----------------------------------------------------------------------------
site_summary <- do.call(
  rbind,
  lapply(levels(analysis_data$site), function(site_level) {
    z <- analysis_data[analysis_data$site == site_level, , drop = FALSE]
    data.frame(
      site = site_level,
      n = nrow(z),
      n_younger = sum(z$age_stratum == "Younger"),
      n_older = sum(z$age_stratum == "Older"),
      mean_log2_qn = mean(z$sarcosine_log2_qn),
      sd_log2_qn = sd(z$sarcosine_log2_qn),
      median_log2_qn = median(z$sarcosine_log2_qn),
      q1_log2_qn = unname(quantile(z$sarcosine_log2_qn, 0.25)),
      q3_log2_qn = unname(quantile(z$sarcosine_log2_qn, 0.75)),
      stringsAsFactors = FALSE
    )
  })
)

primary <- merge(primary, site_summary, by = "site", all.x = TRUE, sort = FALSE)
primary <- primary[match(c("RCC", "LCC", "Rectum"), primary$site), ]
write.csv(primary, OUT_PRIMARY, row.names = FALSE, quote = TRUE)
write.csv(heterogeneity, OUT_HETEROGENEITY, row.names = FALSE, quote = TRUE)
write.csv(pairwise, OUT_PAIRWISE, row.names = FALSE, quote = TRUE)
write.csv(age_sensitivity, OUT_AGE_SENSITIVITY, row.names = FALSE, quote = TRUE)
write.csv(analysis_data, OUT_SOURCE, row.names = FALSE, quote = TRUE)

write.table(
  data.frame(
    File = c("analysis/intensity_matrices.rds", "metadata_merged_619samples.tsv"),
    SHA256 = c(
      digest::digest(INPUT_RDS, algo = "sha256", file = TRUE),
      digest::digest(INPUT_METADATA, algo = "sha256", file = TRUE)
    ),
    stringsAsFactors = FALSE
  ),
  OUT_CHECKSUM,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------------------------------------------------------
# Compact main-figure panel: adjusted site-versus-Control effects
# -----------------------------------------------------------------------------
plot_data <- primary
plot_data$site_label <- c(
  sprintf("Right colon\n(n = %d)", plot_data$n_site[plot_data$site == "RCC"]),
  sprintf("Left colon\n(n = %d)", plot_data$n_site[plot_data$site == "LCC"]),
  sprintf("Rectum\n(n = %d)", plot_data$n_site[plot_data$site == "Rectum"])
)
plot_data$site_label <- factor(plot_data$site_label, levels = rev(plot_data$site_label))
plot_data$holm_label <- paste0("Holm P = ", vapply(plot_data$hc3_p_holm_across_3, format_p, character(1)))

heterogeneity_p <- heterogeneity$p[heterogeneity$Test == "HC3 joint Wald: RCC = LCC = Rectum"]

p <- ggplot(plot_data, aes(x = estimate_log2, y = site_label)) +
  geom_vline(
    xintercept = 0,
    linewidth = NC_AXIS_LW,
    linetype = "22",
    colour = "#606060"
  ) +
  geom_errorbar(
    aes(xmin = hc3_ci_95_low, xmax = hc3_ci_95_high),
    orientation = "y",
    width = 0,
    linewidth = NC_AXIS_LW,
    colour = "#303030"
  ) +
  geom_point(
    shape = 21,
    size = 3.2,
    stroke = NC_AXIS_LW,
    fill = COL_CANCER,
    colour = "#202020"
  ) +
  geom_label(
    aes(x = 0.80, label = holm_label),
    family = "Helvetica",
    size = mm_text(NC_ANNOT_PT),
    hjust = 1,
    colour = "black",
    fill = "white",
    linewidth = 0,
    label.padding = unit(0.7, "pt")
  ) +
  annotate(
    "text",
    x = 0.80,
    y = 3.52,
    label = paste0("CRC-site heterogeneity P = ", format_p(heterogeneity_p)),
    family = "Helvetica",
    size = mm_text(NC_ANNOT_PT),
    hjust = 1,
    colour = "black"
  ) +
  scale_x_continuous(
    limits = c(-0.08, 0.82),
    breaks = seq(0, 0.8, by = 0.2),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Adjusted difference in log2 intensity vs Control (n = 246)",
    y = NULL
  ) +
  theme_nc(base_pt = NC_TICK_PT, grid_y = TRUE) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica", size = NC_TITLE_PT),
    axis.text = element_text(family = "Helvetica", size = NC_TICK_PT),
    strip.text = element_text(family = "Helvetica", size = NC_STRIP_PT),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(lineheight = 0.92, margin = margin(r = 3)),
    plot.margin = margin(10, 4, 2, 2)
  )

WIDTH_IN <- 3.25
HEIGHT_IN <- 2.35
save_nc(p, OUT_PNG, WIDTH_IN, HEIGHT_IN, dpi = 600)
ggsave(
  OUT_PDF,
  p,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  device = grDevices::pdf,
  family = "Helvetica",
  useDingbats = FALSE,
  bg = "white"
)
svglite::svglite(OUT_SVG, width = WIDTH_IN, height = HEIGHT_IN, bg = "white")
print(p)
dev.off()

# -----------------------------------------------------------------------------
# Validation report and session snapshot
# -----------------------------------------------------------------------------
fligner <- fligner.test(sarcosine_log2_qn ~ site, data = analysis_data)
cook_threshold <- 4 / nrow(analysis_data)
cook_n <- sum(cooks.distance(fit_primary) > cook_threshold)

validation_lines <- c(
  "Proposed Figure 1c: anatomical-site consistency of fecal sarcosine",
  "",
  sprintf("Input NEG matrix: %d features x %d samples", nrow(intensity), ncol(intensity)),
  sprintf("Input metadata: %d rows; unique sample IDs: %s", nrow(metadata), !anyDuplicated(metadata$MTBLS_Sample_Name)),
  sprintf("Matrix/metadata ID sets identical: %s", setequal(colnames(intensity), metadata$MTBLS_Sample_Name)),
  sprintf("Sarcosine feature rows matching name + HMDB0000271: %d", length(sarc_idx)),
  sprintf("Feature m/z: %s", feature_meta$mass_to_charge[sarc_idx]),
  sprintf("Feature RT (min): %s", feature_meta$retention_time[sarc_idx]),
  "Normalization: zeros -> half global minimum positive; log2; quantile normalization with features in rows and samples in columns.",
  "Main model: sarcosine_log2_qn ~ site + age_stratum.",
  sprintf("Main-analysis n: %d (Control 246; RCC 59; LCC 87; Rectum 159)", nrow(analysis_data)),
  sprintf("Excluded: QC = %d; CRC unspecified site = %d", sum(analysis_all$site == "QC"), sum(analysis_all$site == "CRC_unspecified")),
  sprintf("Primary design full rank: %s (%d/%d)", fit_primary$rank == ncol(model.matrix(fit_primary)), fit_primary$rank, ncol(model.matrix(fit_primary))),
  sprintf("HC3 CRC-site heterogeneity P: %.8g", heterogeneity_p),
  sprintf("Classical partial-F CRC-site heterogeneity P: %.8g", heterogeneity$p[2]),
  sprintf("HC3 site-by-age interaction P: %.8g", unname(interaction_wald["p"])),
  sprintf("Classical site-by-age interaction P: %.8g", interaction_classical$`Pr(>F)`[2]),
  sprintf("Fligner-Killeen variance test across four groups P: %.8g", fligner$p.value),
  sprintf("Observations with Cook distance > 4/n (threshold %.5f): %d", cook_threshold, cook_n),
  sprintf("Bootstrap: %d site-by-age-cell-stratified repetitions; seed %d", BOOT_REPS, BOOT_SEED),
  sprintf("All three HC3 CIs exclude zero: %s", all(primary$hc3_ci_95_low > 0)),
  sprintf("All three bootstrap CIs exclude zero: %s", all(primary$bootstrap_ci_95_low > 0)),
  sprintf("All six age-stratified point estimates are positive: %s", all(age_sensitivity$estimate_log2 > 0)),
  sprintf("PNG bytes: %d", file.info(OUT_PNG)$size),
  sprintf("PDF bytes: %d", file.info(OUT_PDF)$size),
  sprintf("SVG bytes: %d", file.info(OUT_SVG)$size),
  "",
  "Interpretation safeguard: non-significant CRC-site heterogeneity is compatible with similar effects but does not prove equality.",
  "Identification safeguard: the deposited sarcosine annotation is putative and isobaric; no structural confirmation is claimed."
)
writeLines(validation_lines, OUT_VALIDATION)
capture.output(sessionInfo(), file = OUT_SESSION)

required_outputs <- c(
  OUT_PNG, OUT_PDF, OUT_SVG, OUT_PRIMARY, OUT_HETEROGENEITY,
  OUT_PAIRWISE, OUT_AGE_SENSITIVITY, OUT_SOURCE, OUT_CHECKSUM,
  OUT_VALIDATION, OUT_SESSION
)
if (!all(file.exists(required_outputs))) stop("One or more outputs were not written.")
if (any(file.info(required_outputs)$size == 0)) stop("An output file is empty.")

cat("Completed anatomical-site consistency analysis.\n")
cat("Primary statistics:\n")
print(primary[, c(
  "site", "n_site", "estimate_log2", "hc3_ci_95_low", "hc3_ci_95_high",
  "hc3_p_holm_across_3", "bootstrap_ci_95_low", "bootstrap_ci_95_high"
)])
cat("\nCRC-site heterogeneity:\n")
print(heterogeneity)
cat("\nOutputs:\n", paste(required_outputs, collapse = "\n"), "\n")
