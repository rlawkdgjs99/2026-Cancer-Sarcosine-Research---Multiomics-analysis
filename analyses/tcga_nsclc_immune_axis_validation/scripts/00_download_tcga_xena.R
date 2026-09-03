#!/usr/bin/env Rscript
## =====================================================================
## Download TCGA NSCLC expression, survival, and clinical data from UCSC Xena.
## ---------------------------------------------------------------------
## Data source : UCSC Xena tcgaHub  (https://tcga.xenahubs.net)
## Access date : 2026-06-15
## Dataset IDs : verified from UCSCXenaTools::XenaData (not fabricated)
##   NSCLC = TCGA Lung Cancer (LUNG)  [= LUAD + LUSC]
## Expression unit: log2(norm_count+1)  (RSEM, IlluminaHiSeq)
## Survival       : Liu et al. 2018 curated endpoints (OS/DSS/DFI/PFI)
## =====================================================================

suppressMessages(suppressWarnings(library(UCSCXenaTools)))
set.seed(1)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
proj <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
ddir <- file.path(proj, "data", "raw_xena")
dir.create(ddir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(proj, "results", "logs"), recursive = TRUE, showWarnings = FALSE)

data(XenaData, package = "UCSCXenaTools")

## verified dataset IDs (host = tcgaHub)
datasets <- list(
  NSCLC = c(expr = "TCGA.LUNG.sampleMap/HiSeqV2",
            clin = "TCGA.LUNG.sampleMap/LUNG_clinicalMatrix",
            surv = "survival/LUNG_survival.txt")
)

## robust single-dataset getter: pre-filter XenaData (avoids non-standard eval),
## require exactly one match, then query/download/load.
get_one <- function(id) {
  sub <- XenaData[XenaData$XenaHostNames == "tcgaHub" & XenaData$XenaDatasets == id, , drop = FALSE]
  if (nrow(sub) != 1L) stop(sprintf("expected 1 dataset for '%s', found %d", id, nrow(sub)))
  xq  <- XenaQuery(XenaGenerate(XenaData = sub))
  xd  <- XenaDownload(xq, destdir = ddir, download_probeMap = FALSE, force = FALSE)
  obj <- XenaPrepare(xd)
  if (is.list(obj) && !is.data.frame(obj)) obj <- obj[[1]]
  obj
}

raw <- list()
for (co in names(datasets)) {
  raw[[co]] <- list()
  for (ty in names(datasets[[co]])) {
    id <- datasets[[co]][[ty]]
    cat(sprintf("[%-5s/%-4s] %s\n", co, ty, id))
    raw[[co]][[ty]] <- get_one(id)
  }
}

saveRDS(raw, file.path(proj, "data", "01_xena_raw.rds"))

## ----------------------------- verification ----------------------------
cat("\n================ VERIFICATION ================\n")
for (co in names(raw)) {
  cat("\n##########", co, "##########\n")
  e <- raw[[co]]$expr; s <- raw[[co]]$surv; cl <- raw[[co]]$clin

  cat("EXPR class:", class(e)[1], "| dim:", paste(dim(e), collapse = " x "),
      "| id-col:", colnames(e)[1], "\n")
  cat("     first genes  :", paste(utils::head(e[[1]], 3), collapse = ", "), "\n")
  bc <- colnames(e)[-1]
  cat("     first samples:", paste(utils::head(bc, 3), collapse = ", "), "\n")
  cat("     sample-type code (barcode pos 14-15):\n")
  print(table(substr(bc, 14, 15)))

  cat("SURV class:", class(s)[1], "| dim:", paste(dim(s), collapse = " x "), "\n")
  cat("     colnames:", paste(colnames(s), collapse = ", "), "\n")
  surv_id <- as.character(s[[1]])
  cat("     first surv id:", surv_id[1], "\n")

  cat("CLIN class:", class(cl)[1], "| dim:", paste(dim(cl), collapse = " x "), "\n")

  cat("OVERLAP expr<->surv  full barcode :", length(intersect(bc, surv_id)),
      "/", length(bc), "expr samples\n")
  cat("OVERLAP expr<->surv  12-char (pt) :",
      length(intersect(substr(bc, 1, 12), substr(surv_id, 1, 12))), "\n")
}
cat("\nSaved -> data/01_xena_raw.rds\n")
writeLines(capture.output(sessionInfo()),
           file.path(proj, "results", "logs", "sessionInfo_download.txt"))
