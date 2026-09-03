## 01_load_prepare.R
## Load HGMT export -> subset to 16S -> build genus-level TSS relative-abundance matrix.
## Workflow rationale:
##  - The Bacteria file is LONG (Taxa | Run ID | Abundance) and mixes 235 16S + 25 WGS samples.
##  - 16S resolves only to genus (no s__/t__); WGS uses MetaPhlAn4 SGB levels.
##  - We keep ONLY 16S samples and ONLY genus-terminal rows, then TSS-normalize per sample.

source("analysis/00_setup.R")

log_con <- file(file.path(LOG_DIR, "01_load_prepare.log"), open = "wt")
sink(log_con, split = TRUE)
cat("==== 01_load_prepare ====", as.character(Sys.time()), "\n")

## ---------- Metadata ----------
## Layout: line 1 = "Selected project id: ..."; line 2 = header (23 fields);
## lines 3+ = data rows that carry ONE extra unnamed trailing column (24 fields,
## a verified export artifact mirroring QC State). We therefore skip BOTH the
## title and header lines and assign the 23 known names + a dummy 24th.
meta_cols <- c("Project ID","Run ID","Sample name","Assay type","Sequencing method",
               "Phenotype name","Sample description","Phenotype ID","country",
               "Geographic location","Longitude","lattitude","age","sex","BMI",
               "Antibiotic use","Antibiotic name","Bodysite","Batch","Batch2",
               "QC bacteria","QC fungi","QC State","extra_trailing")
meta_raw <- fread(META_FILE, sep = "\t", header = FALSE, skip = 2,
                  col.names = meta_cols, quote = "", na.strings = c("NA", ""))
cat("metadata raw dim:", nrow(meta_raw), "x", ncol(meta_raw), "\n")
stopifnot(nrow(meta_raw) == 260, "Run ID" %in% names(meta_raw))
## positional-integrity guard: assay must be only 16S/WGS (catches any column shift)
stopifnot(all(meta_raw$`Assay type` %in% c("16S", "WGS")))

meta <- meta_raw[, .(
  run_id    = `Run ID`,
  sample    = `Sample name`,
  assay     = `Assay type`,
  phenotype = `Phenotype name`,
  country   = country,
  age       = suppressWarnings(as.numeric(age)),
  sex       = sex,
  bmi       = suppressWarnings(as.numeric(BMI)),
  qc_bac    = `QC bacteria`
)]

meta16 <- meta[assay == "16S"]
cat("16S samples in metadata:", nrow(meta16), "\n")
stopifnot(nrow(meta16) == 235)

## group factor: Health = reference (control)
meta16[, group := fifelse(phenotype == "Health", "Health",
                  fifelse(phenotype == "Carcinoma, Non-Small-Cell Lung", "NSCLC", NA_character_))]
stopifnot(all(!is.na(meta16$group)))
meta16[, group := factor(group, levels = GROUP_LEVELS)]
cat("group counts:\n"); print(table(meta16$group))
stopifnot(all(meta16$qc_bac == 1))             # all 16S pass bacterial QC
stopifnot(!any(duplicated(meta16$run_id)))     # one row per run
cat("age NA:", sum(is.na(meta16$age)), "| bmi NA:", sum(is.na(meta16$bmi)),
    "| sex NA:", sum(is.na(meta16$sex)), "\n")

## ---------- Bacteria (taxonomic abundance) ----------
bac <- fread(BACTERIA_FILE, sep = "\t", header = TRUE, quote = "")
setnames(bac, c("taxa", "run_id", "abund"))
cat("bacteria rows (all assays):", nrow(bac), "\n")

bac <- bac[run_id %in% meta16$run_id]          # keep only the 235 16S samples
## genus-terminal rows only: contain |g__ and NOT |s__ / |t__ (drops higher partial lineages too)
bac_g <- bac[grepl("\\|g__", taxa) & !grepl("\\|s__|\\|t__", taxa)]
cat("genus-terminal rows (16S):", nrow(bac_g), "\n")
stopifnot(nrow(bac_g) > 0)

## wide: genera (rows) x samples (cols); sum any accidental duplicates explicitly
wide <- dcast(bac_g, taxa ~ run_id, value.var = "abund", fun.aggregate = sum, fill = 0)
gmat <- as.matrix(wide[, -1])
rownames(gmat) <- wide$taxa
cat("genus matrix:", nrow(gmat), "genera x", ncol(gmat), "samples\n")
stopifnot(ncol(gmat) == 235)

## align metadata rows to matrix columns
meta16 <- meta16[match(colnames(gmat), run_id)]
stopifnot(all(meta16$run_id == colnames(gmat)))

## ---------- TSS normalization (relative abundance of recognizable genera; HGMT paper) ----------
csum <- colSums(gmat)
cat("pre-TSS per-sample genus sums (%): min", round(min(csum), 2),
    "max", round(max(csum), 2), "median", round(median(csum), 2), "\n")
stopifnot(all(csum > 0))
gtss <- sweep(gmat, 2, csum, "/")              # each sample's genera sum to 1
stopifnot(all(abs(colSums(gtss) - 1) < 1e-8))

## ---------- feature lookup (safe id <-> genus <-> full lineage) ----------
genus_short <- sub(".*\\|g__", "", rownames(gtss))
miss <- (genus_short == "" | is.na(genus_short))
if (any(miss)) genus_short[miss] <- paste0("unclassified_genus_", seq_len(sum(miss)))
safe_id <- make.unique(gsub("[^A-Za-z0-9]", "_", genus_short))
feat_lookup <- data.table(feature = safe_id, genus = genus_short, lineage = rownames(gtss))
fwrite(feat_lookup, file.path(DERIVED_DIR, "feature_lookup.csv"))

## orient samples x genera for downstream (diversity / MaAsLin2)
genus_tss <- t(gtss); colnames(genus_tss) <- safe_id   # TSS proportions
genus_raw <- t(gmat); colnames(genus_raw) <- safe_id   # raw HGMT rel-abund (%)

saveRDS(list(meta = meta16, genus_tss = genus_tss, genus_raw = genus_raw,
             feature_lookup = feat_lookup),
        file.path(DERIVED_DIR, "prepared_16S.rds"))
fwrite(meta16, file.path(DERIVED_DIR, "metadata_16S.csv"))

## ---------- sanity summary ----------
cat("\n[SANITY]\n")
cat("samples:", nrow(genus_tss), "| genera:", ncol(genus_tss), "\n")
cat("group:", paste(names(table(meta16$group)), as.integer(table(meta16$group)),
                    collapse = " | "), "\n")
cat("any NA in TSS matrix:", any(is.na(genus_tss)), "\n")
cat("TSS row-sum range:", round(min(rowSums(genus_tss)), 6), "-",
    round(max(rowSums(genus_tss)), 6), "\n")
cat("genera present (>0) in >=10% samples:",
    sum(colMeans(genus_raw > 0) >= PREVALENCE_MIN), "\n")
cat("==== done ====\n")
sink(); close(log_con)
