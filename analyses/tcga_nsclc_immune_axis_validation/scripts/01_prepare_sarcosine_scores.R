#!/usr/bin/env Rscript
## =====================================================================
## Prepare TCGA NSCLC primary tumors and the RNA-defined sarcosine scores.
## ---------------------------------------------------------------------
## Inputs : data/01_xena_raw.rds   (expression = log2(norm_count+1), symbols)
## Output : data/02_scores.rds     (per-cohort score data.frame, primary tumors)
##
## Score definitions & sources:
##  - Sarcosine axis (KEGG-verified): production = GNMT, DMGDH ;
##                                     degradation = SARDH, PIPOX
##      production_z / degradation_z = mean of per-gene z-scores
##      balance_z   = production_z - degradation_z   (higher => sarcosine-accumulating)
##      balance_log2ratio = mean(log2 prod genes) - mean(log2 deg genes)
##  - CYT score (Rooney et al., Cell 2015): geometric mean of GZMA & PRF1.
##      NOTE: Rooney used TPM; Xena HiSeqV2 is RSEM log2(norm_count+1).
##      We de-log to norm_count, take geometric mean, then log2(+1). Documented analog.
##  - ssGSEA (GSVA ssgseaParam): Cytotoxicity & CD8 T cells = Danaher et al.,
##      J Immunother Cancer 2017; CD8_cytotoxic_panel = user-specified markers;
##      Exhaustion = PDCD1/HAVCR2/LAG3/TOX (study design).
##  - ssGSEA is rank-based within sample => invariant to the log transform.
## =====================================================================

suppressMessages(suppressWarnings({library(GSVA)}))
set.seed(1)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
proj <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
raw  <- readRDS(file.path(proj, "data", "01_xena_raw.rds"))

## ---- gene sets ----
sarco <- list(production  = c("GNMT","DMGDH"),
              degradation = c("SARDH","PIPOX"))
markers_cd8  <- c("CD8A","GZMA","GZMB","PRF1","IFNG","NKG7")
markers_exh  <- c("PDCD1","HAVCR2","LAG3","TOX")
ssgsea_sets <- list(
  Cytotoxicity_Danaher = c("CTSW","GNLY","GZMA","GZMB","GZMH","KLRB1","KLRD1","KLRK1","NKG7","PRF1"),
  CD8_Tcells_Danaher   = c("CD8A","CD8B"),
  CD8_cytotoxic_panel  = markers_cd8,
  Exhaustion           = markers_exh
)
all_query_genes <- unique(c(unlist(sarco), markers_cd8, markers_exh, unlist(ssgsea_sets)))

z <- function(v) { s <- sd(v, na.rm = TRUE); if (is.na(s) || s == 0) return(rep(NA_real_, length(v))); (v - mean(v, na.rm = TRUE)) / s }

build_scores <- function(e, cohort) {
  g <- e[[1]]
  m <- as.matrix(e[, -1]); rownames(m) <- g; storage.mode(m) <- "double"
  ## primary tumors only (barcode sample-type 01), de-duplicate
  m <- m[, substr(colnames(m), 14, 15) == "01", drop = FALSE]
  m <- m[, !duplicated(colnames(m)), drop = FALSE]
  m <- m[!duplicated(rownames(m)), , drop = FALSE]
  cat(sprintf("\n[%s] primary-tumor matrix: %d genes x %d samples | total NA: %d\n",
              cohort, nrow(m), ncol(m), sum(is.na(m))))
  ## gene presence check
  miss <- setdiff(all_query_genes, rownames(m))
  cat("   missing query genes:", if (length(miss)) paste(miss, collapse = ", ") else "none", "\n")

  gv <- function(gene) if (gene %in% rownames(m)) m[gene, ] else rep(NA_real_, ncol(m))

  ## sarcosine axis
  prod_z <- colMeans(rbind(z(gv("GNMT")), z(gv("DMGDH"))), na.rm = TRUE)
  deg_z  <- colMeans(rbind(z(gv("SARDH")), z(gv("PIPOX"))), na.rm = TRUE)
  balance_z <- prod_z - deg_z
  balance_log2ratio <- colMeans(rbind(gv("GNMT"), gv("DMGDH")), na.rm = TRUE) -
                       colMeans(rbind(gv("SARDH"), gv("PIPOX")), na.rm = TRUE)

  ## CYT (Rooney 2015) geometric mean of GZMA, PRF1 on de-logged norm_count
  lin <- function(x) pmax(2^x - 1, 0)
  CYT <- log2(sqrt(lin(gv("GZMA")) * lin(gv("PRF1"))) + 1)

  ## exhaustion composite (mean z of 4 markers)
  exh_z <- colMeans(do.call(rbind, lapply(markers_exh, function(gg) z(gv(gg)))), na.rm = TRUE)

  ## ssGSEA (GSVA) - restrict sets to present genes
  gs <- lapply(ssgsea_sets, function(s) intersect(s, rownames(m)))
  par <- ssgseaParam(m, gs, minSize = 1, normalize = TRUE)
  es  <- gsva(par, verbose = FALSE)                       # signatures x samples
  es  <- t(es); colnames(es) <- paste0("ssgsea_", colnames(es))

  df <- data.frame(
    sample = colnames(m),
    GNMT = gv("GNMT"), DMGDH = gv("DMGDH"), SARDH = gv("SARDH"), PIPOX = gv("PIPOX"),
    production_z = prod_z, degradation_z = deg_z,
    balance_z = balance_z, balance_log2ratio = balance_log2ratio,
    CD8A = gv("CD8A"), GZMA = gv("GZMA"), GZMB = gv("GZMB"),
    PRF1 = gv("PRF1"), IFNG = gv("IFNG"), NKG7 = gv("NKG7"),
    CYT = CYT,
    PDCD1 = gv("PDCD1"), HAVCR2 = gv("HAVCR2"), LAG3 = gv("LAG3"), TOX = gv("TOX"),
    exhaustion_z = exh_z,
    es[colnames(m), , drop = FALSE],
    row.names = NULL, check.names = FALSE
  )
  ## sanity summary
  cat("   score summary (median [min, max]):\n")
  for (cc in c("balance_z","CYT","CD8A","ssgsea_CD8_cytotoxic_panel","exhaustion_z"))
    cat(sprintf("     %-28s %.3f [%.3f, %.3f]\n", cc,
                median(df[[cc]], na.rm = TRUE), min(df[[cc]], na.rm = TRUE), max(df[[cc]], na.rm = TRUE)))
  df
}

scores <- list()
for (co in names(raw)) scores[[co]] <- build_scores(raw[[co]]$expr, co)
saveRDS(scores, file.path(proj, "data", "02_scores.rds"))
cat("\nSaved -> data/02_scores.rds\n")
