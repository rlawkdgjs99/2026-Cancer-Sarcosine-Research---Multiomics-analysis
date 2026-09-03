#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
data_root <- normalizePath(file.path(analysis_root, ".."), mustWork = TRUE)
prior_dir <- file.path(analysis_root, "02_lineage_reannotation")
.libPaths(c(file.path(analysis_dir, "R_libs"), file.path(prior_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})

options(stringsAsFactors = FALSE, warn = 1, future.globals.maxSize = 32 * 1024^3)
set.seed(260825)

table_dir <- file.path(analysis_dir, "results", "tables")
pub_dir <- file.path(analysis_dir, "results", "figures_publication")
log_dir <- file.path(analysis_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) sub(" .*", "", system2("shasum", c("-a", "256", path), stdout = TRUE))
record <- list()
check <- function(name, condition, observed, expected) {
  ok <- isTRUE(condition)
  record[[length(record) + 1L]] <<- data.table(
    check = name, status = if (ok) "PASS" else "FAIL",
    observed = as.character(observed), expected = as.character(expected)
  )
  if (!ok) stop("Verification failed: ", name, "; observed=", observed, "; expected=", expected)
}

counts_file <- file.path(prior_dir, "intermediate", "01_full_counts_mt20_qc.rds")
lineage_file <- file.path(prior_dir, "results", "tables", "09_final_cell_lineages_FROZEN.csv")
score_file <- file.path(analysis_dir, "intermediate", "01_cell_paper_style_scores.rds")
deg_file <- file.path(table_dir, "02_high_vs_low_cell_FindMarkers.csv")
go_file <- file.path(table_dir, "03_high_low_GO_BP_overrepresentation.csv")

check("count SHA-256", sha256(counts_file) == "7552ba882c306fa68380390a8172e0fc96af2c7be7cb3f5778345204e9fd61b2",
      sha256(counts_file), "7552ba882c306fa68380390a8172e0fc96af2c7be7cb3f5778345204e9fd61b2")
check("lineage SHA-256", sha256(lineage_file) == "97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1",
      sha256(lineage_file), "97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1")
check("twice-reproduced DEG SHA-256",
      sha256(deg_file) == "af5b154e2d16434be435204ab511449a1713206773fc4d62087a8325e669dc09",
      sha256(deg_file), "af5b154e2d16434be435204ab511449a1713206773fc4d62087a8325e669dc09")

message("Loading frozen inputs and saved scores")
counts <- readRDS(counts_file)
lin <- fread(lineage_file)
score <- readRDS(score_file)
check("count dimensions", identical(dim(counts), c(24292L, 92053L)), paste(dim(counts), collapse = "x"), "24292x92053")
check("unique count cell IDs", anyDuplicated(colnames(counts)) == 0L, anyDuplicated(colnames(counts)), 0)
check("unique lineage cell IDs", anyDuplicated(lin$cell_id) == 0L, anyDuplicated(lin$cell_id), 0)
check("count/lineage exact ID set", setequal(colnames(counts), lin$cell_id),
      length(intersect(colnames(counts), lin$cell_id)), 92053)
check("score/count exact cell order", identical(score$cell_id, colnames(counts)),
      sum(score$cell_id == colnames(counts)), 92053)

targets <- c("GNMT", "DMGDH", "SARDH", "PIPOX")
check("all four module genes unique", all(targets %in% rownames(counts)) && !anyDuplicated(rownames(counts)),
      paste(targets[targets %in% rownames(counts)], collapse = ";"), paste(targets, collapse = ";"))

# Independent LogNormalize verification on deterministic cells and the four module genes.
verify_idx <- unique(c(1:8, seq.int(1000L, ncol(counts), length.out = 12L), (ncol(counts) - 7L):ncol(counts)))
verify_idx <- as.integer(round(verify_idx))
lib <- Matrix::colSums(counts[, verify_idx, drop = FALSE])
expected_norm <- log1p(sweep(as.matrix(counts[targets, verify_idx, drop = FALSE]), 2, lib / 10000, "/"))
saved_norm <- rbind(
  GNMT = score$GNMT_log1pCP10k[verify_idx],
  DMGDH = score$DMGDH_log1pCP10k[verify_idx],
  SARDH = score$SARDH_log1pCP10k[verify_idx],
  PIPOX = score$PIPOX_log1pCP10k[verify_idx]
)
lognorm_diff <- max(abs(expected_norm - saved_norm))
check("independent LogNormalize subset", lognorm_diff <= 1e-12,
      format(lognorm_diff, scientific = TRUE), "<=1e-12")

message("Rebuilding a clean Seurat object for independent AddModuleScore verification")
obj <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize",
                     scale.factor = 10000, margin = 1, verbose = FALSE)
obj <- AddModuleScore(
  obj,
  features = list(Production = c("GNMT", "DMGDH"), Degradation = c("SARDH", "PIPOX")),
  pool = rownames(obj), nbin = 24, ctrl = 100, k = FALSE, assay = "RNA",
  name = "SarcosineModule", seed = 260825, search = FALSE, slot = "data"
)
recalc <- obj[[]]
prod_diff <- max(abs(recalc$SarcosineModule1 - score$production_module_raw))
deg_diff <- max(abs(recalc$SarcosineModule2 - score$degradation_module_raw))
check("clean-object Production AddModuleScore", prod_diff <= 1e-12,
      format(prod_diff, scientific = TRUE), "<=1e-12")
check("clean-object Degradation AddModuleScore", deg_diff <= 1e-12,
      format(deg_diff, scientific = TRUE), "<=1e-12")
rm(obj, recalc)
invisible(gc())

expected_prod_shift <- score$production_module_raw - min(score$production_module_raw)
expected_deg_shift <- score$degradation_module_raw - min(score$degradation_module_raw)
check("Production minimum shift", max(abs(expected_prod_shift - score$production_module_shifted)) <= 1e-15,
      max(abs(expected_prod_shift - score$production_module_shifted)), "<=1e-15")
check("Degradation minimum shift", max(abs(expected_deg_shift - score$degradation_module_shifted)) <= 1e-15,
      max(abs(expected_deg_shift - score$degradation_module_shifted)), "<=1e-15")

defined <- expected_deg_shift > 0
expected_ratio <- expected_prod_shift[defined] / expected_deg_shift[defined]
ratio_diff <- max(abs(expected_ratio - score$production_degradation_ratio[defined]))
check("ratio formula", ratio_diff <= 1e-15, format(ratio_diff, scientific = TRUE), "<=1e-15")
check("undefined denominator count", sum(!defined) == 1L, sum(!defined), 1)
ratio_median <- median(expected_ratio)
expected_group <- ifelse(!defined, "Undefined denominator",
                         ifelse(score$production_degradation_ratio > ratio_median, "High",
                                ifelse(score$production_degradation_ratio < ratio_median, "Low", "At median")))
check("median ratio", abs(ratio_median - 0.343967676527376) <= 1e-15,
      format(ratio_median, digits = 17), "0.343967676527376")
check("High/Low membership", identical(expected_group, score$ratio_group),
      sum(expected_group == score$ratio_group), 92053)
check("High cell count", sum(score$ratio_group == "High") == 46026L, sum(score$ratio_group == "High"), 46026)
check("Low cell count", sum(score$ratio_group == "Low") == 46026L, sum(score$ratio_group == "Low"), 46026)

deg <- fread(deg_file)
expected_direction <- fifelse(
  deg$p_val_adj < 0.05 & deg$avg_log2FC >= 0.25, "High",
  fifelse(deg$p_val_adj < 0.05 & deg$avg_log2FC <= -0.25, "Low", "Not selected")
)
check("DEG direction rule", identical(expected_direction, deg$direction),
      sum(expected_direction == deg$direction), nrow(deg))
check("selected High genes", sum(deg$direction == "High") == 5625L, sum(deg$direction == "High"), 5625)
check("selected Low genes", sum(deg$direction == "Low") == 1185L, sum(deg$direction == "Low"), 1185)

go <- fread(go_file)
bh_high_diff <- max(abs(p.adjust(go$P.High, method = "BH") - go$BH_High))
bh_low_diff <- max(abs(p.adjust(go$P.Low, method = "BH") - go$BH_Low))
check("GO BH High family", bh_high_diff <= 1e-12, format(bh_high_diff, scientific = TRUE), "<=1e-12")
check("GO BH Low family", bh_low_diff <= 1e-12, format(bh_low_diff, scientific = TRUE), "<=1e-12")

spec <- fread(file.path(table_dir, "04_figure_specifications.csv"))
for (i in seq_len(nrow(spec))) {
  png_file <- file.path(pub_dir, paste0(spec$stem[i], ".png"))
  pdf_file <- file.path(pub_dir, paste0(spec$stem[i], ".pdf"))
  check(paste0(spec$stem[i], " PNG exists"), file.exists(png_file), file.exists(png_file), TRUE)
  check(paste0(spec$stem[i], " PDF exists"), file.exists(pdf_file), file.exists(pdf_file), TRUE)
  sips_out <- system2("sips", c("-g", "pixelWidth", "-g", "pixelHeight", png_file), stdout = TRUE)
  width_px <- as.integer(sub(".*: ", "", sips_out[grepl("pixelWidth", sips_out)]))
  height_px <- as.integer(sub(".*: ", "", sips_out[grepl("pixelHeight", sips_out)]))
  check(paste0(spec$stem[i], " width"), width_px == spec$png_expected_width_px[i],
        width_px, spec$png_expected_width_px[i])
  check(paste0(spec$stem[i], " height"), height_px == spec$png_expected_height_px[i],
        height_px, spec$png_expected_height_px[i])
}

verification <- rbindlist(record, fill = TRUE)
fwrite(verification, file.path(table_dir, "05_verification_summary.csv"))
capture.output(sessionInfo(), file = file.path(log_dir, "02_verification_sessionInfo.txt"))
cat("FIG3-STYLE VERIFICATION PASS\n")
cat("Checks:", nrow(verification), "\n")
cat("Max LogNormalize difference:", format(lognorm_diff, scientific = TRUE), "\n")
cat("Max AddModuleScore differences Production/Degradation:",
    format(prod_diff, scientific = TRUE), "/", format(deg_diff, scientific = TRUE), "\n")
