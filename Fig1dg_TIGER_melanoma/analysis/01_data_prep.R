# =============================================================================
# 01_data_prep.R
# Project: Sarcosine metabolism genes vs ICI response/survival (TIGER, Melanoma PRJEB23709)
# Purpose: Load + verify expression and clinical data, build the analysis data frame.
#
# Target genes (user-specified):
#   Degradation enzymes : SARDH, PIPOX
#   Production enzymes   : GNMT, DMGDH
#
# CRITICAL data quirk (verified): both TSVs have an off-by-one header
#   (R write.table row.names style) -> header has one FEWER field than data rows.
#   We therefore read with header=FALSE and re-attach names manually, then PROVE
#   correct alignment against ground-truth values extracted independently via awk.
# =============================================================================

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
})

# ---- Paths (relative to project root; run with cwd = project root) ----------
expr_file <- "Melanoma-PRJEB23709_ExpressionData.tsv"
clin_file <- "Melanoma-PRJEB23709_ClinicalData.tsv"
stopifnot(file.exists(expr_file), file.exists(clin_file))

# =============================================================================
# 1. EXPRESSION MATRIX
#    Header line: GENE_SYMBOL + 91 sample IDs            (92 fields)
#    Data rows  : <int index> + gene symbol + 91 values  (93 fields)
# =============================================================================
expr_header <- strsplit(readLines(expr_file, n = 1L), "\t")[[1]]
cat("Expression header fields:", length(expr_header), "\n")
stopifnot(expr_header[1] == "GENE_SYMBOL")
sample_ids_expr <- expr_header[-1]                       # 91 sample IDs (ERR...)
cat("  -> sample columns (from header):", length(sample_ids_expr), "\n")

# Read DATA rows only (skip header), so columns are: V1=index, V2=gene, V3..=values
expr_raw <- data.table::fread(
  expr_file, header = FALSE, sep = "\t", skip = 1L,
  showProgress = FALSE, data.table = FALSE
)
cat("Expression data rows x cols:", nrow(expr_raw), "x", ncol(expr_raw), "\n")
# Expect: ncol = 1 (index) + 1 (gene) + 91 (samples) = 93
stopifnot(ncol(expr_raw) == length(sample_ids_expr) + 2L)

gene_symbols <- as.character(expr_raw[[2]])              # V2 = gene symbol
expr_mat <- as.matrix(expr_raw[, 3:ncol(expr_raw)])     # V3..V93 = 91 sample values
colnames(expr_mat) <- sample_ids_expr                   # header[-1] aligns to V3..V93
storage.mode(expr_mat) <- "double"
cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")
cat("Value range:", paste(round(range(expr_mat, na.rm = TRUE), 4), collapse = " .. "),
    "| any NA:", any(is.na(expr_mat)), "\n")

# ---- GROUND-TRUTH alignment check (values extracted independently via awk) ----
# awk gave, for the FIRST two sample columns (ERR2208887, ERR2208888):
#   GNMT  : 4.49209634703616 , 5.75870662414432
#   SARDH : 0.756248331863607, 1.06339753002665
#   PIPOX : 0.201461209221552, 0.094927278894717
#   DMGDH : 0.907233194559215, 0.330014124981244
get_gene <- function(g) {
  idx <- which(gene_symbols == g)
  if (length(idx) != 1L) stop(sprintf("Gene '%s' matched %d rows (expected 1).", g, length(idx)))
  expr_mat[idx, , drop = TRUE]
}
truth <- list(
  GNMT  = c(ERR2208887 = 4.49209634703616, ERR2208888 = 5.75870662414432),
  SARDH = c(ERR2208887 = 0.756248331863607, ERR2208888 = 1.06339753002665),
  PIPOX = c(ERR2208887 = 0.201461209221552, ERR2208888 = 0.094927278894717),
  DMGDH = c(ERR2208887 = 0.907233194559215, ERR2208888 = 0.330014124981244)
)
for (g in names(truth)) {
  v <- get_gene(g)
  chk <- v[c("ERR2208887", "ERR2208888")]
  if (!isTRUE(all.equal(unname(chk), unname(truth[[g]]), tolerance = 1e-6))) {
    stop(sprintf("ALIGNMENT FAIL for %s: got %s, expected %s",
                 g, paste(chk, collapse=","), paste(truth[[g]], collapse=",")))
  }
  cat(sprintf("  alignment OK: %-5s ERR2208887=%.5f ERR2208888=%.5f\n", g, chk[1], chk[2]))
}

target_genes <- c("SARDH", "PIPOX", "GNMT", "DMGDH")
stopifnot(all(target_genes %in% gene_symbols))
# Extract target gene FPKM matrix (genes x samples)
gene_fpkm <- t(sapply(target_genes, get_gene))          # 4 x 91
cat("Target gene FPKM matrix:", nrow(gene_fpkm), "x", ncol(gene_fpkm), "\n")

# =============================================================================
# 2. CLINICAL DATA  (same off-by-one: header 17 fields, data 18 fields)
#    data V1=index; V2..V18 map to clin_header[1..17]
# =============================================================================
clin_header <- strsplit(readLines(clin_file, n = 1L), "\t")[[1]]
cat("\nClinical header fields:", length(clin_header), "\n")
clin_raw <- data.table::fread(
  clin_file, header = FALSE, sep = "\t", skip = 1L,
  na.strings = c("NA", ""), showProgress = FALSE, data.table = FALSE
)
cat("Clinical data rows x cols:", nrow(clin_raw), "x", ncol(clin_raw), "\n")
stopifnot(ncol(clin_raw) == length(clin_header) + 1L)   # +1 for leading index col
colnames(clin_raw) <- c("row_index", clin_header)       # re-attach correct names

# Pull the columns we need by their proper names
clin <- data.frame(
  sample_id     = as.character(clin_raw[["sample_id"]]),
  patient_name  = as.character(clin_raw[["patient_name"]]),
  timepoint     = as.character(clin_raw[["Treatment"]]),          # PRE / EDT
  recist        = as.character(clin_raw[["response"]]),           # CR/PR/SD/PD
  response_NR   = as.character(clin_raw[["response_NR"]]),        # R / N
  os_days       = as.numeric(clin_raw[["overall survival (days)"]]),
  vital_status  = as.character(clin_raw[["vital status"]]),       # Alive / Dead
  gender        = as.character(clin_raw[["Gender"]]),
  therapy       = as.character(clin_raw[["Therapy"]]),
  age           = as.numeric(clin_raw[["age_start"]]),
  stringsAsFactors = FALSE
)
cat("Clinical parsed:", nrow(clin), "rows\n")

# ---- GROUND-TRUTH clinical check (row 1 = sample ERR2208887, from awk) -------
# sample_id=ERR2208887, patient=ipiPD1_1, timepoint=EDT, recist=PD, response_NR=N,
# os_days=689, vital=Alive, gender=Female, therapy=anti-CTLA-4+anti-PD-1, age=42
r1 <- clin[clin$sample_id == "ERR2208887", ]
stopifnot(nrow(r1) == 1L)
stopifnot(r1$patient_name == "ipiPD1_1", r1$timepoint == "EDT", r1$recist == "PD",
          r1$response_NR == "N", r1$os_days == 689, r1$vital_status == "Alive",
          r1$gender == "Female", r1$therapy == "anti-CTLA-4+anti-PD-1", r1$age == 42)
cat("  clinical alignment OK (ERR2208887 row matches ground truth)\n")

# ---- Sample ID concordance: expression <-> clinical -------------------------
stopifnot(setequal(colnames(gene_fpkm), clin$sample_id))
cat("Sample ID concordance OK: 91 samples match 1:1 between expression and clinical\n")

# =============================================================================
# 3. BUILD ANALYSIS DATA FRAME (PRE-only cohort)
# =============================================================================
# Order clinical to match expression columns, then attach gene FPKM
clin <- clin[match(colnames(gene_fpkm), clin$sample_id), ]
stopifnot(identical(clin$sample_id, colnames(gene_fpkm)))

# log2(FPKM + 1) transform for visualization, comparison, and z-scoring
log_fpkm <- log2(gene_fpkm + 1)                          # 4 genes x 91

dat <- clin
for (g in target_genes) dat[[g]] <- as.numeric(log_fpkm[g, ])   # log2(FPKM+1) per gene

# Response group factor: R = responder (CR/PR), NR = non-responder (SD/PD)
dat$response_group <- factor(ifelse(dat$response_NR == "R", "R", "NR"),
                             levels = c("R", "NR"))
# Survival event: Dead = 1, Alive = 0
dat$event <- ifelse(dat$vital_status == "Dead", 1L,
                    ifelse(dat$vital_status == "Alive", 0L, NA_integer_))
stopifnot(!any(is.na(dat$event)))
# Therapy short labels
dat$therapy_short <- ifelse(dat$therapy == "anti-PD-1", "antiPD1",
                     ifelse(dat$therapy == "anti-CTLA-4+anti-PD-1", "combo", NA))
stopifnot(!any(is.na(dat$therapy_short)))

# ---- Subset to PRE (baseline) cohort ----------------------------------------
cat("\nBefore timepoint filter:", nrow(dat), "samples\n")
dat_pre <- dat[dat$timepoint == "PRE", ]
cat("After PRE-only filter   :", nrow(dat_pre), "samples\n")
stopifnot(nrow(dat_pre) == 73L)
# PRE cohort: confirm one sample per patient (no pseudoreplication)
stopifnot(!any(duplicated(dat_pre$patient_name)))
cat("PRE cohort: all", length(unique(dat_pre$patient_name)), "patients unique (no duplicates)\n")

# =============================================================================
# 4. FUNCTIONAL MODULE SCORES (z-score averaging within PRE cohort)
#    Rationale: GNMT (FPKM ~4-5) >> DMGDH (~0.3-0.9); a simple mean would be
#    dominated by the higher-expressed gene. z-scoring each gene's log2(FPKM+1)
#    to mean 0 / sd 1 lets both genes contribute equally, then average per module.
#    Scores are computed WITHIN the PRE analysis cohort (the population analyzed).
# =============================================================================
z <- function(x) (x - mean(x)) / sd(x)
deg_genes  <- c("SARDH", "PIPOX")   # degradation
prod_genes <- c("GNMT", "DMGDH")    # production
dat_pre$Degradation_score <- rowMeans(sapply(deg_genes,  function(g) z(dat_pre[[g]])))
dat_pre$Production_score  <- rowMeans(sapply(prod_genes, function(g) z(dat_pre[[g]])))
cat("\nModule scores added (z-score mean of log2(FPKM+1)):\n")
cat("  Degradation (SARDH,PIPOX): mean =", round(mean(dat_pre$Degradation_score),4),
    " sd =", round(sd(dat_pre$Degradation_score),4), "\n")
cat("  Production  (GNMT,DMGDH) : mean =", round(mean(dat_pre$Production_score),4),
    " sd =", round(sd(dat_pre$Production_score),4), "\n")

# =============================================================================
# 5. SANITY SUMMARY of the analysis cohort
# =============================================================================
cat("\n================ PRE COHORT SANITY SUMMARY ================\n")
cat("n =", nrow(dat_pre), "\n")
cat("\nresponse_group:\n"); print(table(dat_pre$response_group))
cat("\ntherapy x response_group:\n"); print(table(dat_pre$therapy_short, dat_pre$response_group))
cat("\ntherapy x event (Dead=1):\n"); print(table(dat_pre$therapy_short, dat_pre$event))
cat("\nOS days summary:\n"); print(summary(dat_pre$os_days))
cat("\nlog2(FPKM+1) per-gene summary (PRE):\n")
print(round(t(sapply(target_genes, function(g) summary(dat_pre[[g]]))), 3))

# =============================================================================
# 6. SAVE
# =============================================================================
saveRDS(dat_pre, "results/analysis_data_PRE.rds")
saveRDS(dat,     "results/analysis_data_all91.rds")
write.csv(dat_pre, "results/tables/analysis_data_PRE.csv", row.names = FALSE)
cat("\nSaved: results/analysis_data_PRE.rds (n=73), results/analysis_data_all91.rds (n=91)\n")
cat("Saved: results/tables/analysis_data_PRE.csv\n")

cat("\n--- sessionInfo ---\n")
print(sessionInfo())
