#!/usr/bin/env Rscript

# Independent verification of script 113 outputs. This script reloads the frozen
# within-cohort inputs, independently recomputes every eligible random-effects
# meta-analysis, BH correction, display selection and disease-enrichment overlap,
# then checks the final PNG and PowerPoint transport copy.

suppressPackageStartupMessages(library(metafor))

argv <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", argv[grepl("^--file=", argv)])
if (length(script_file) != 1L) stop("Run with Rscript")
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
project_root <- normalizePath(file.path(root, "..", "..", ".."), mustWork = TRUE)

cross_dir <- file.path(root, "results_integrated", "cross_cohort")
out_dir <- file.path(root, "results_integrated", "CRC_WGS_PROD_DEG_META_HEATMAP_26.08.31")
prod_file <- file.path(cross_dir, "cross_cohort_prod_assoc_correlations.csv")
deg_file <- file.path(cross_dir, "cross_cohort_deg_assoc_correlations.csv")
cancer_file <- file.path(cross_dir, "cross_cohort_CRC_enriched_membership_matrix.csv")
healthy_file <- file.path(root, "results_integrated", "SupFig3de_HEALTHY_ENRICHED_26.08.27", "healthy_enriched_membership_matrix.csv")
transport_file <- file.path(project_root, "Manuscript작업", "FigDesign&Manuscript", "PPT_INSERT_CRC_WGS_ProdDeg_Heatmap.png")
png_file <- file.path(out_dir, "Fig_CRC_WGS_species_Production_Degradation_meta_heatmap.png")

checks <- character()
record <- function(name, condition) {
  if (!isTRUE(condition)) stop("FAILED: ", name)
  checks <<- c(checks, paste0("PASS\t", name))
}

prod <- read.csv(prod_file, check.names = FALSE)
deg <- read.csv(deg_file, check.names = FALSE)
meta_saved <- read.csv(file.path(out_dir, "meta_correlations_all_species.csv"), check.names = FALSE)
selected_saved <- read.csv(file.path(out_dir, "selected_top15_per_axis.csv"), check.names = FALSE)
plot_saved <- read.csv(file.path(out_dir, "heatmap_plot_values_long.csv"), check.names = FALSE)
overlap_saved <- read.csv(file.path(out_dir, "selected_species_disease_enrichment_overlap.csv"), check.names = FALSE)
input_md5_saved <- read.csv(file.path(out_dir, "input_file_md5.csv"), check.names = FALSE)

required_inputs <- c(prod_file, deg_file, cancer_file, healthy_file)
md5_now <- unname(tools::md5sum(required_inputs))
record("all four input files match the recorded MD5 values", identical(md5_now, input_md5_saved$md5))
record("frozen production input has 1,815 rows", nrow(prod) == 1815L)
record("frozen degradation input has 1,815 rows", nrow(deg) == 1815L)

to_long <- function(x, axis) {
  data.frame(
    Axis = axis,
    Cohort = x$Cohort,
    Species_full = x$Species_full,
    Species = x$Species,
    rho = x$rho,
    n_samples = x$n_samples,
    stringsAsFactors = FALSE
  )
}
cor_long <- rbind(to_long(prod, "Production"), to_long(deg, "Degradation"))
key <- paste(cor_long$Axis, cor_long$Species_full, sep = "|||IDX|||")
split_input <- split(cor_long, key)
eligible <- split_input[vapply(split_input, function(x) length(unique(x$Cohort)) >= 3L, logical(1))]

recomputed <- lapply(eligible, function(dat) {
  dat <- dat[order(dat$Cohort), , drop = FALSE]
  rho <- pmin(pmax(dat$rho, -0.999999), 0.999999)
  fit <- metafor::rma.uni(yi = atanh(rho), vi = 1 / (dat$n_samples - 3), method = "REML")
  data.frame(
    Axis = dat$Axis[1],
    Species_full = dat$Species_full[1],
    Species = dat$Species[1],
    k_cohorts = length(unique(dat$Cohort)),
    meta_rho_recomputed = tanh(as.numeric(fit$b)),
    meta_p_recomputed = as.numeric(fit$pval),
    I2_recomputed = as.numeric(fit$I2),
    stringsAsFactors = FALSE
  )
})
recomputed <- do.call(rbind, recomputed)
recomputed$meta_q_recomputed <- ave(
  recomputed$meta_p_recomputed,
  recomputed$Axis,
  FUN = function(x) p.adjust(x, method = "BH")
)

cmp <- merge(
  meta_saved[, c("Axis", "Species_full", "Species", "k_cohorts", "meta_rho", "meta_p", "meta_q", "I2")],
  recomputed,
  by = c("Axis", "Species_full", "Species", "k_cohorts"),
  all = TRUE,
  sort = TRUE
)
record("independent meta-analysis has the same 656 eligible axis/species rows", nrow(cmp) == 656L && !anyNA(cmp$Axis))
record("all meta-rho values independently match within 1e-12", max(abs(cmp$meta_rho - cmp$meta_rho_recomputed)) < 1e-12)
record("all meta P values independently match within 1e-12", max(abs(cmp$meta_p - cmp$meta_p_recomputed)) < 1e-12)
record("all BH q values independently match within 1e-12", max(abs(cmp$meta_q - cmp$meta_q_recomputed)) < 1e-12)
record("all I2 values independently match within 1e-10", max(abs(cmp$I2 - cmp$I2_recomputed)) < 1e-10)

select_independently <- function(axis) {
  x <- recomputed[recomputed$Axis == axis & recomputed$meta_rho_recomputed > 0 & recomputed$meta_q_recomputed < 0.05, ]
  x <- x[order(-x$meta_rho_recomputed, x$meta_q_recomputed, x$Species), ]
  head(x$Species_full, 15L)
}
prod_selected <- select_independently("Production")
deg_selected <- select_independently("Degradation")
record("Production top-15 selection independently matches", identical(prod_selected, selected_saved$Species_full[selected_saved$Axis == "Production"]))
record("Degradation top-15 selection independently matches", identical(deg_selected, selected_saved$Species_full[selected_saved$Axis == "Degradation"]))
record("the display contains exactly 30 unique species and 60 cells", length(unique(plot_saved$Species_full)) == 30L && nrow(plot_saved) == 60L)
record("all displayed heatmap values are finite", all(is.finite(plot_saved$meta_rho)))

cancer <- read.csv(cancer_file, check.names = FALSE)
healthy <- read.csv(healthy_file, check.names = FALSE)
selected_ids <- overlap_saved$Species_full
expected_cancer <- selected_ids %in% cancer$Species_full[cancer$n_cohorts_significant >= 3L]
expected_healthy <- selected_ids %in% healthy$Species_full[healthy$n_cohorts_strict >= 3L]
record("recurrent Cancer-enrichment overlap is independently mapped", identical(overlap_saved$Cancer_recurrent_n >= 3L, expected_cancer))
record("recurrent Healthy-enrichment overlap is independently mapped", identical(overlap_saved$Healthy_recurrent_n >= 3L, expected_healthy))
record("exactly three selected species are recurrently Cancer-enriched", sum(expected_cancer) == 3L)
record("no selected species is recurrently Healthy-enriched", sum(expected_healthy) == 0L)

record("analytical and PowerPoint PNG files exist", file.exists(png_file) && file.exists(transport_file))
record("analytical and PowerPoint PNG hashes are identical", identical(unname(tools::md5sum(png_file)), unname(tools::md5sum(transport_file))))
file_info <- system2("/usr/bin/file", shQuote(transport_file), stdout = TRUE, stderr = TRUE)
record("PowerPoint PNG is non-interlaced 8-bit RGB at 3800x4800", any(grepl("3800 x 4800, 8-bit/color RGB, non-interlaced", file_info, fixed = TRUE)))
sips_info <- system2(
  "/usr/bin/sips",
  c("-g", "profile", "-g", "dpiWidth", "-g", "dpiHeight", shQuote(transport_file)),
  stdout = TRUE,
  stderr = TRUE
)
record("PowerPoint PNG has embedded sRGB IEC61966-2.1", any(grepl("profile: sRGB IEC61966-2.1", sips_info, fixed = TRUE)))
record("PowerPoint PNG has exact 400 dpi metadata", any(grepl("dpiWidth: 400.000", sips_info, fixed = TRUE)) && any(grepl("dpiHeight: 400.000", sips_info, fixed = TRUE)))
xattrs <- system2("/usr/bin/xattr", c("-l", shQuote(transport_file)), stdout = TRUE, stderr = TRUE)
record("PowerPoint PNG has no extended attributes", length(xattrs) == 0L)

report <- c(
  "Independent verification: CRC WGS Production/Degradation meta-heatmap",
  paste0("Verified: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  checks,
  paste0("TOTAL_PASS\t", length(checks)),
  "TOTAL_FAIL\t0"
)
writeLines(report, file.path(out_dir, "independent_verification_report.txt"))
cat(paste(report, collapse = "\n"), "\n")
