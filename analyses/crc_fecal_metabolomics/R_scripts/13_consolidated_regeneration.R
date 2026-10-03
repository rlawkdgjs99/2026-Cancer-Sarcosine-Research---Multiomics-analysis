# Consolidated regeneration: MaAsLin2 merge, replication fix, clean outputs, updated report
# All steps reproducible from saved RDS files — no FTP, no .wiff
#
# 2026-04-29 fixes (review):
#   FIX 1: QN axis corrected to features × samples (was: samples × features)
#   FIX 2: age_cohort_bin added as MaAsLin2 + limma covariate (older/younger)
#   FIX 3: original metabolite names restored from MaAsLin2 sanitized output
#   FIX 4: extended drug filter (Amisulpride, Prilocaine, Trimethadione, etc.)
#   FIX 5: lipid-class replication (don't match different chain lengths via KEGG)
sessionInfo_path <- file.path("R_scripts", paste0("sessionInfo_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
log_file <- file.path("analysis", sprintf("09_consolidated_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
sink(log_file, split = TRUE)

library(limma)
library(ggplot2)
library(readxl)

# Output routing (2026 reorg): Main vs Supplementary deliverables.
# Main = headline figures + the differential tables they are built from.
# Supplementary = everything else (replication, MaAsLin2 raw, neighbour plots, QC).
dir_main <- "analysis/main"
dir_supp <- "analysis/supplementary"
dir.create(dir_main, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_supp, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 0. Load base data
# =============================================================================
cat("=== Loading data ===\n")
dat <- readRDS("analysis/intensity_matrices.rds")
metadata <- read.table("metadata_merged_619samples.tsv", header = TRUE, sep = "\t",
    check.names = FALSE, stringsAsFactors = FALSE)

load_maf <- function(path) {
    read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
               stringsAsFactors = FALSE, comment.char = "", quote = "")
}
pos_maf <- load_maf("MTBLS10232_inventory/m_MTBLS10232_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv")
neg_maf <- load_maf("MTBLS10232_inventory/m_MTBLS10232_LC-MS_negative_reverse-phase_metabolite_profiling_v2_maf.tsv")

# Build HMDB→KEGG cross-ref
pos_pairs <- pos_maf[!is.na(pos_maf$keggdatabase) & pos_maf$keggdatabase != "" & pos_maf$keggdatabase != "null" &
                     !is.na(pos_maf$database_identifier) & pos_maf$database_identifier != "" & pos_maf$database_identifier != "null",
                     c("database_identifier", "keggdatabase")]
pos_pairs <- pos_pairs[!duplicated(pos_pairs$database_identifier), ]
hmdb2kegg <- setNames(pos_pairs$keggdatabase, pos_pairs$database_identifier)
cat(sprintf("HMDB->KEGG cross-ref: %d pairs\n", length(hmdb2kegg)))

# Build maf lookup
build_lookup <- function(maf) {
    lookup <- data.frame(
        metabolite = maf$metabolite_identification,
        HMDB_ID = maf$database_identifier,
        KEGG_ID_maf = maf$keggdatabase,
        mass_to_charge = maf$mass_to_charge,
        retention_time = maf$retention_time,
        stringsAsFactors = FALSE
    )
    lookup <- lookup[!is.na(lookup$metabolite) & lookup$metabolite != "", ]
    lookup <- lookup[!duplicated(lookup$metabolite), ]
    # Fill KEGG from cross-ref where maf KEGG is missing
    missing_kegg <- is.na(lookup$KEGG_ID_maf) | lookup$KEGG_ID_maf == "" | lookup$KEGG_ID_maf == "null"
    lookup$KEGG_ID <- lookup$KEGG_ID_maf
    lookup$KEGG_ID[missing_kegg] <- hmdb2kegg[lookup$HMDB_ID[missing_kegg]]
    lookup$KEGG_ID <- as.character(lookup$KEGG_ID)
    lookup$KEGG_ID_maf <- NULL
    lookup
}
lookup_pos <- build_lookup(pos_maf)
lookup_neg <- build_lookup(neg_maf)

# =============================================================================
# 1. Annotation confidence flagging function
# =============================================================================
build_annotation_flag <- function(feature_names) {
    n <- length(feature_names)
    flag <- rep("endogenous_likely", n)
    
    drug_pats <- c("cilastatin", "ofloxacin", "amifloxacin", "pantoprazole", "ramiprilat",
        "antibiotic", "irbesartan$", "valsartan$", "eprosartan", "ibuprofen",
        "naproxen", "hydromorphone", "glimepiride", "metformin",
        "mycophenolic acid", "bexarotene", "tamoxifen", "fluticasone",
        "dexamethasone", "prednisolone$", "prednisone$", "tiagabine",
        "olopatadine", "indacaterol", "diltiazem", "verapamil",
        "lidocaine", "bupivacaine", "omeprazole", "esomeprazole",
        "cetirizine", "fexofenadine", "loratadine", "sertraline",
        "fluoxetine", "citalopram", "risperidone", "quetiapine",
        # FIX 4 (2026-04-29): drugs that appeared in top hits but were missed
        "amisulpride", "sulpiride", "tiapride",                    # antipsychotics
        "prilocaine", "procaine", "tetracaine",                    # anesthetics
        "trimethadione", "paramethadione", "ethadione",            # anticonvulsants
        "caffeine", "theobromine", "theophylline", "paraxanthine", # xanthines (food/drug)
        "acetaminophen", "paracetamol",                            # OTC analgesic
        "salicyl",                                                 # aspirin family
        "warfarin", "clopidogrel", "ticagrelor",                   # anticoagulants
        "atenolol", "metoprolol", "bisoprolol", "propranolol",     # beta-blockers
        "amlodipine", "nifedipine",                                # Ca-channel blockers
        "ranitidine", "famotidine", "cimetidine",                  # H2 antagonists
        "losartan",                                                # angiotensin antag
        "verapamil")
    
    contam_pats <- c("perfluoro", "perfluorooctane", "pfos$", "pfoa$", "phthalate",
        "bisphenol", "atrazine", "glyphosate", "chlorpyrifos",
        "imidacloprid", "ddt$", "dde$", "polychlorinated", "polybrominated",
        "dioxin", "furan", "cotinine", "nicotine", "n.nitroso",
        "roquefortine", "zearalenone", "aflatoxin", "ochratoxin",
        "deoxynivalenol", "fumonisin", "citrinin", "patulin",
        "c.i..food", "food.red", "food.yellow", "food.blue",
        "brilliant.blue", "tartrazine", "sunset.yellow",
        "sucralose", "aspartame", "saccharin", "acesulfame",
        "paraben", "triclosan", "benzophenone", "oxybenzone")
    
    plant_pats <- c("gibberellin", "lucidone", "ganoderic", "ganoderenic",
        "austalide", "cinncassiol", "cincassiol", "cinnzeylanine",
        "melleolide", "tanabalin", "kanzonol", "glyurallin", "glyzaglabrin",
        "glycyrrh", "artocarpetin", "artocarpesin", "moracin",
        "medicagenic acid", "soyasapogenol", "soyasaponin",
        "ustiloxin", "polyporusterone", "bolegrevilol",
        "monacolin", "lovastatin")
    
    for (pat in drug_pats) flag[grepl(pat, feature_names, ignore.case=TRUE, perl=TRUE)] <- "drug_related"
    for (pat in contam_pats) flag[grepl(pat, feature_names, ignore.case=TRUE, perl=TRUE)] <- "contaminant_xenobiotic"
    for (pat in plant_pats) flag[grepl(pat, feature_names, ignore.case=TRUE, perl=TRUE)] <- "plant_exogenous"
    
    flag
}

# =============================================================================
# 2. Re-run full differential analysis (all 3 methods) and build clean outputs
# =============================================================================
cat("\n=== Re-running differential analysis (all 3 methods) ===\n")

run_full_analysis <- function(intensity_mat, feature_meta, sample_names,
                                maf_lookup, polarity_label) {

    cat(sprintf("\n--- %s ---\n", polarity_label))

    # FIX 1: QN on features × samples (correct axis)
    intensity_imp <- intensity_mat
    min_nz <- min(intensity_mat[intensity_mat > 0])
    intensity_imp[intensity_imp == 0] <- min_nz / 2
    log_mat <- log2(intensity_imp)                        # features × samples
    qn_mat  <- normalizeBetweenArrays(log_mat, method = "quantile")

    # Transpose to samples × features for downstream
    t_mat_qn <- t(qn_mat)
    feat_names <- feature_meta$metabolite_identification
    colnames(t_mat_qn) <- feat_names
    rownames(t_mat_qn) <- sample_names

    # Diagnostic: confirm sample-level normalization (per-sample SD should ≈ 0)
    per_sample_mean <- rowMeans(t_mat_qn)
    cat(sprintf("  QN diagnostic — per-sample mean SD: %.4f (should be ~0 for correct sample QN)\n",
        sd(per_sample_mean)))

    # Match metadata
    m_idx <- match(rownames(t_mat_qn), metadata$MTBLS_Sample_Name)
    stopifnot(sum(!is.na(m_idx)) == nrow(t_mat_qn))
    meta <- metadata[m_idx, ]

    # Drop QCs for testing
    is_patient <- meta$sample_type == "Patient"
    t_patients <- t_mat_qn[is_patient, ]
    meta_patients <- meta[is_patient, ]

    group   <- ifelse(grepl("CRC", meta_patients$cohort_subset), "CRC", "CTRL")
    age_bin <- ifelse(grepl("^older_", meta_patients$cohort_subset), "older", "younger")
    n_crc <- sum(group == "CRC"); n_ctrl <- sum(group == "CTRL")
    cat(sprintf("  CRC=%d CTRL=%d  Older=%d Younger=%d\n", n_crc, n_ctrl,
        sum(age_bin=="older"), sum(age_bin=="younger")))

    # FIX 2: limma with age_cohort_bin covariate
    group_f <- factor(group, levels = c("CTRL", "CRC"))
    age_f   <- factor(age_bin, levels = c("younger", "older"))
    design <- model.matrix(~ group_f + age_f)
    fit <- lmFit(t(t_patients), design)
    fit2 <- eBayes(fit, trend = TRUE)
    res <- topTable(fit2, coef = "group_fCRC", number = Inf, sort.by = "none")
    res$metabolite <- colnames(t_patients)

    # Wilcoxon (unchanged — non-parametric, no covariate adjustment available cleanly)
    crc_idx <- which(group == "CRC"); ctrl_idx <- which(group == "CTRL")
    w_pvals <- apply(t_patients, 2, function(x) {
        suppressWarnings(wilcox.test(x[crc_idx], x[ctrl_idx], exact = FALSE)$p.value)
    })
    w_padj <- p.adjust(w_pvals, method = "BH")
    w_fc <- colMeans(t_patients[crc_idx, , drop = FALSE]) - colMeans(t_patients[ctrl_idx, , drop = FALSE])

    # FIX 2 + FIX 3: Re-run MaAsLin2 with covariate, restore original names
    maaslin_res <- NULL
    if (requireNamespace("Maaslin2", quietly = TRUE)) {
        suppressMessages(library(Maaslin2))
        feat_table <- data.frame(t(t_patients), check.names = FALSE)
        meta_maaslin <- data.frame(
            sample = rownames(t_patients),
            group = group,
            age_cohort_bin = age_bin,
            stringsAsFactors = FALSE
        )
        rownames(meta_maaslin) <- meta_maaslin$sample
        out_dir <- file.path(dir_supp, paste0("maaslin2_", polarity_label))
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
        # Clean previous outputs
        for (f in list.files(out_dir, full.names = TRUE)) try(file.remove(f), silent = TRUE)

        tryCatch({
            invisible(Maaslin2(
                input_data = feat_table, input_metadata = meta_maaslin, output = out_dir,
                fixed_effects = c("group", "age_cohort_bin"),
                reference = c("group,CTRL", "age_cohort_bin,younger"),
                normalization = "NONE", transform = "NONE", standardize = FALSE,
                plot_heatmap = FALSE, plot_scatter = FALSE, cores = 1
            ))
            mfile <- file.path(out_dir, "all_results.tsv")
            if (file.exists(mfile)) {
                m2 <- read.table(mfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
                # Keep only the group=CRC contrast
                m2 <- m2[m2$metadata == "group" & m2$value == "CRC", ]

                # FIX 3: restore original feature names (MaAsLin2 sanitizes via make.names)
                name_map <- setNames(colnames(t_patients), make.names(colnames(t_patients)))
                orig_names <- name_map[m2$feature]
                missing_n <- sum(is.na(orig_names))
                if (missing_n > 0) {
                    orig_names[is.na(orig_names)] <- m2$feature[is.na(orig_names)]
                }
                cat(sprintf("  MaAsLin2 name restoration: %d/%d mapped, %d unmapped\n",
                    sum(!is.na(name_map[m2$feature])), nrow(m2), missing_n))

                maaslin_res <- data.frame(
                    metabolite = orig_names,
                    log2FC = m2$coef,
                    pvalue = m2$pval,
                    padj = m2$qval,
                    method = "MaAsLin2",
                    stringsAsFactors = FALSE
                )
            }
        }, error = function(e) {
            cat("  MaAsLin2 error:", conditionMessage(e), "\n")
        })
    }
    
    # Annotation flags (uses feat_names from above; do NOT redeclare to avoid shadowing)
    anno_flag <- build_annotation_flag(colnames(t_patients))
    
    # Build unified output
    out_list <- list(
        data.frame(metabolite = res$metabolite, log2FC = res$logFC,
                   pvalue = res$P.Value, padj = res$adj.P.Val, method = "limma",
                   annotation_confidence = anno_flag, stringsAsFactors = FALSE),
        data.frame(metabolite = colnames(t_patients), log2FC = w_fc,
                   pvalue = w_pvals, padj = w_padj, method = "wilcoxon",
                   annotation_confidence = anno_flag, stringsAsFactors = FALSE)
    )
    if (!is.null(maaslin_res)) {
        # Map annotation flags onto MaAsLin2 result (by restored original metabolite name)
        maaslin_res$annotation_confidence <- anno_flag[match(maaslin_res$metabolite, colnames(t_patients))]
        maaslin_res$annotation_confidence[is.na(maaslin_res$annotation_confidence)] <- "endogenous_likely"
        out_list[[3]] <- maaslin_res
    }
    
    out <- do.call(rbind, out_list)
    
    # Join chemical metadata
    out <- merge(out, maf_lookup, by = "metabolite", all.x = TRUE)
    
    # Stats
    out_limma <- out[out$method == "limma", ]
    n_features <- nrow(out_limma)
    n_sig <- sum(out_limma$padj < 0.05 & abs(out_limma$log2FC) > 0.5)
    n_end <- sum(out_limma$annotation_confidence == "endogenous_likely")
    n_end_sig <- sum(out_limma$padj < 0.05 & abs(out_limma$log2FC) > 0.5 &
                      out_limma$annotation_confidence == "endogenous_likely")
    
    cat(sprintf("  Total features: %d\n", n_features))
    cat(sprintf("  Significant (padj<0.05, |log2FC|>0.5): %d\n", n_sig))
    cat(sprintf("  Endogenous: %d, Endogenous sig: %d\n", n_end, n_end_sig))
    
    fname <- file.path(dir_main, paste0("diff_results_", polarity_label, ".tsv"))
    write.table(out, fname, sep = "\t", row.names = FALSE, quote = TRUE)
    cat(sprintf("  Saved %s (%d rows x %d cols)\n", fname, nrow(out), ncol(out)))
    
    out
}

res_pos <- run_full_analysis(dat$intensity_pos, dat$feature_meta_pos,
                               dat$sample_cols_pos, lookup_pos, "POS")
res_neg <- run_full_analysis(dat$intensity_neg, dat$feature_meta_neg,
                               dat$sample_cols_neg, lookup_neg, "NEG")

# =============================================================================
# 3. Build significant_features_annotated.tsv (endogenous only, combined)
# =============================================================================
cat("\n=== Building significant_features_annotated.tsv ===\n")

limma_pos <- res_pos[res_pos$method == "limma", ]
limma_neg <- res_neg[res_neg$method == "limma", ]

sig_pos_end <- limma_pos[limma_pos$padj < 0.05 & abs(limma_pos$log2FC) > 0.5 &
                           limma_pos$annotation_confidence == "endogenous_likely", ]
sig_neg_end <- limma_neg[limma_neg$padj < 0.05 & abs(limma_neg$log2FC) > 0.5 &
                           limma_neg$annotation_confidence == "endogenous_likely", ]

combined_sig <- rbind(
    cbind(sig_pos_end, polarity = "POS", stringsAsFactors = FALSE),
    cbind(sig_neg_end, polarity = "NEG", stringsAsFactors = FALSE)
)

# Mark cross-mode evidence
hmdb_both <- intersect(
    sig_pos_end$HMDB_ID[!is.na(sig_pos_end$HMDB_ID) & sig_pos_end$HMDB_ID != "" & sig_pos_end$HMDB_ID != "null"],
    sig_neg_end$HMDB_ID[!is.na(sig_neg_end$HMDB_ID) & sig_neg_end$HMDB_ID != "" & sig_neg_end$HMDB_ID != "null"]
)
combined_sig$cross_mode <- combined_sig$HMDB_ID %in% hmdb_both &
    !is.na(combined_sig$HMDB_ID) & combined_sig$HMDB_ID != "" & combined_sig$HMDB_ID != "null"

write.table(combined_sig, file.path(dir_main, "significant_features_annotated.tsv"), sep = "\t",
    row.names = FALSE, quote = TRUE)
cat(sprintf("Saved significant_features_annotated.tsv: %d rows, %d cross-mode\n",
    nrow(combined_sig), sum(combined_sig$cross_mode)))
cat(sprintf("  POS endogenous sig: %d, NEG endogenous sig: %d\n",
    nrow(sig_pos_end), nrow(sig_neg_end)))
cat(sprintf("  Cross-mode HMDB IDs: %d\n", length(hmdb_both)))

# =============================================================================
# 4. Replication check (with GlcCer fix)
# =============================================================================
cat("\n=== Replication check vs Li 2024 Table S4 ===\n")

li2024 <- read_excel("Other_files/Supp info from the Article/mmc3.xlsx", sheet = "TableS4")
hdrs <- as.character(li2024[1, ])
li2024 <- li2024[-1, ]
colnames(li2024) <- hdrs

li_met <- trimws(li2024$metabolite)
li_lfc <- suppressWarnings(as.numeric(li2024$log2FC))
li_alt <- li2024$Alteration
li_group <- li2024$Group
li_kegg <- trimws(li2024[["kegg compound"]])

all_matches <- data.frame(
    our_metabolite = character(), our_log2FC = numeric(), our_padj = numeric(),
    li_metabolite = character(), li_log2FC = numeric(), li_Alteration = character(),
    li_group = character(), polarity = character(), match_type = character(),
    direction_match = logical(), stringsAsFactors = FALSE
)

for (pol in c("POS", "NEG")) {
    our_sig <- if (pol == "POS") sig_pos_end else sig_neg_end
    
    # Name match
    our_names <- tolower(our_sig$metabolite)
    li_names <- tolower(li_met[!is.na(li_met)])
    name_matches <- intersect(our_names, li_names)
    
    for (nm in name_matches) {
        our_row <- our_sig[our_names == nm, ][1, ]
        li_idx <- which(li_names == nm)[1]
        li_fc <- li_lfc[li_idx]
        
        # Determine Li direction
        if (!is.na(li_fc)) li_dir <- sign(li_fc)
        else if (!is.na(li_alt[li_idx]) && li_alt[li_idx] == "UP") li_dir <- 1
        else if (!is.na(li_alt[li_idx]) && li_alt[li_idx] == "DOWN") li_dir <- -1
        else li_dir <- NA
        
        all_matches <- rbind(all_matches, data.frame(
            our_metabolite = our_row$metabolite, our_log2FC = our_row$log2FC,
            our_padj = our_row$padj,
            li_metabolite = li_met[li_idx],
            li_log2FC = if (!is.na(li_fc)) li_fc else NA_real_,
            li_Alteration = if (!is.na(li_alt[li_idx])) li_alt[li_idx] else "",
            li_group = if (!is.na(li_group[li_idx])) li_group[li_idx] else "",
            polarity = pol, match_type = "name",
            direction_match = !is.na(li_dir) && sign(our_row$log2FC) == li_dir,
            stringsAsFactors = FALSE
        ))
    }
    
    # KEGG match (additive — different names, same KEGG ID)
    # FIX 5 (2026-04-29): suppress lipid-class KEGG matches that compare different
    # acyl-chain compositions (e.g., GlcCer 34:3 vs GlcCer 16:0 are different species
    # but share the GlcCer KEGG class ID — direction comparison is meaningless).
    is_lipid_class_collision <- function(our_name, li_name) {
        # Patterns indicating lipid-class names with different chain composition
        # If both names match a lipid-class pattern AND have different "X:Y" specs, suppress.
        lipid_pat <- "(GlcCer|Cer|PC|PE|PS|PI|PG|TG|DG|MG|LPC|LPE|SM|CL)"
        if (!grepl(lipid_pat, our_name, ignore.case=TRUE) ||
            !grepl(lipid_pat, li_name,  ignore.case=TRUE)) return(FALSE)
        # Extract X:Y chain compositions
        get_chains <- function(s) {
            m <- regmatches(s, gregexpr("[0-9]+:[0-9]+", s))[[1]]
            if (length(m) == 0) return("") else return(paste(sort(m), collapse="|"))
        }
        oc <- get_chains(our_name); lc <- get_chains(li_name)
        if (oc == "" || lc == "") return(FALSE)
        return(oc != lc)
    }

    our_kegg <- unique(our_sig$KEGG_ID[!is.na(our_sig$KEGG_ID) & our_sig$KEGG_ID != "" & our_sig$KEGG_ID != "null"])
    li_kegg_clean <- unique(li_kegg[!is.na(li_kegg) & li_kegg != ""])
    kegg_matches <- intersect(our_kegg, li_kegg_clean)

    n_class_collisions <- 0
    for (kid in kegg_matches) {
        our_rows <- our_sig[our_sig$KEGG_ID == kid & !is.na(our_sig$KEGG_ID), ]
        li_rows <- which(li_kegg == kid & !is.na(li_kegg))
        for (oi in seq_len(min(nrow(our_rows), 1))) {
            already_matched <- tolower(our_rows$metabolite[oi]) %in% name_matches
            if (already_matched) next

            our_nm <- our_rows$metabolite[oi]
            li_nm  <- li_met[li_rows[1]]

            # Skip lipid-class collisions
            if (is_lipid_class_collision(our_nm, li_nm)) {
                n_class_collisions <- n_class_collisions + 1
                next
            }

            li_fc <- li_lfc[li_rows[1]]
            if (!is.na(li_fc)) li_dir <- sign(li_fc)
            else if (!is.na(li_alt[li_rows[1]]) && li_alt[li_rows[1]] == "UP") li_dir <- 1
            else if (!is.na(li_alt[li_rows[1]]) && li_alt[li_rows[1]] == "DOWN") li_dir <- -1
            else li_dir <- NA

            all_matches <- rbind(all_matches, data.frame(
                our_metabolite = our_nm,
                our_log2FC = our_rows$log2FC[oi], our_padj = our_rows$padj[oi],
                li_metabolite = li_nm,
                li_log2FC = if (!is.na(li_fc)) li_fc else NA_real_,
                li_Alteration = if (!is.na(li_alt[li_rows[1]])) li_alt[li_rows[1]] else "",
                li_group = if (!is.na(li_group[li_rows[1]])) li_group[li_rows[1]] else "",
                polarity = pol, match_type = "kegg_id",
                direction_match = !is.na(li_dir) && sign(our_rows$log2FC[oi]) == li_dir,
                stringsAsFactors = FALSE
            ))
        }
    }
    if (n_class_collisions > 0) {
        cat(sprintf("  Suppressed %d lipid-class KEGG collisions in %s\n", n_class_collisions, pol))
    }
}

# Report
n_total <- nrow(all_matches)
n_agree <- sum(all_matches$direction_match)
n_disagree <- n_total - n_agree

cat(sprintf("Total matched metabolites: %d\n", n_total))
cat(sprintf("Direction agree: %d (%.1f%%)\n", n_agree, n_agree/n_total*100))
cat(sprintf("Direction disagree: %d\n", n_disagree))

if (n_disagree > 0) {
    cat("\nDisagreeing metabolites:\n")
    disagree_rows <- all_matches[!all_matches$direction_match, ]
    for (i in seq_len(nrow(disagree_rows))) {
        cat(sprintf("  %s (our=%.3f, li=%.3f)\n",
            disagree_rows$our_metabolite[i], disagree_rows$our_log2FC[i],
            disagree_rows$li_log2FC[i]))
    }
}

write.table(all_matches, file.path(dir_supp, "replication_vs_Li2024_combined.tsv"), sep = "\t",
    row.names = FALSE, quote = TRUE)
cat("Saved", file.path(dir_supp, "replication_vs_Li2024_combined.tsv"), "\n")

# Also save per-polarity
for (pol in c("POS", "NEG")) {
    sub <- all_matches[all_matches$polarity == pol, ]
    write.table(sub, file.path(dir_supp, paste0("replication_vs_Li2024_", pol, ".tsv")),
        sep = "\t", row.names = FALSE, quote = TRUE)
}

# =============================================================================
# 5. Volcano plots with annotation coloring
# =============================================================================
cat("\n=== Regenerating volcano plots ===\n")

for (pol in c("POS", "NEG")) {
    df <- if (pol == "POS") limma_pos else limma_neg
    volcano_df <- df
    volcano_df$neg_log10_padj <- -log10(volcano_df$padj)
    volcano_df$neg_log10_padj[!is.finite(volcano_df$neg_log10_padj)] <- max(
        volcano_df$neg_log10_padj[is.finite(volcano_df$neg_log10_padj)], na.rm = TRUE)
    volcano_df$significant <- with(volcano_df, padj < 0.05 & abs(log2FC) > 0.5)
    # Endogenous-only counts for the subtitle (FIXED: previously labeled but not filtered)
    end_sig <- volcano_df$significant & volcano_df$annotation_confidence == "endogenous_likely"
    n_up_end   <- sum(end_sig & volcano_df$log2FC > 0)
    n_down_end <- sum(end_sig & volcano_df$log2FC < 0)
    n_up_all   <- sum(volcano_df$significant & volcano_df$log2FC > 0)
    n_down_all <- sum(volcano_df$significant & volcano_df$log2FC < 0)

    p <- ggplot(volcano_df, aes(x = log2FC, y = neg_log10_padj,
            color = annotation_confidence)) +
        geom_point(size = 0.8, alpha = 0.5) +
        scale_color_manual(values = c(
            "endogenous_likely" = "grey60",
            "plant_exogenous" = "#E69F00",
            "drug_related" = "#D55E00",
            "contaminant_xenobiotic" = "#CC0000"
        )) +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
        labs(title = paste(pol, "- CRC vs CTRL (annotation-filtered volcano)"),
             subtitle = sprintf("Endogenous sig: Up=%d Down=%d  |  All sig: Up=%d Down=%d  |  Red=contaminant Orange=drug Yellow=plant",
                                n_up_end, n_down_end, n_up_all, n_down_all),
             x = "limma log2 Fold Change, CRC / CTRL (age-adjusted)",
             y = "-log10(adjusted p-value)",
             color = "Annotation") +
        theme_minimal(base_size = 16) +
        theme(legend.position = "bottom")

    ggsave(file.path(dir_main, paste0("volcano_", pol, ".png")),
           p, width = 10, height = 8, dpi = 300)
    cat(sprintf("  Saved %s/volcano_%s.png\n", dir_main, pol))
}

# =============================================================================
# 6. Summary table for report
# =============================================================================
cat("\n=== Summary statistics for report ===\n")
cat(sprintf("POS features: %d, sig: %d, endogenous sig: %d\n",
    nrow(limma_pos), sum(limma_pos$padj<0.05 & abs(limma_pos$log2FC)>0.5),
    nrow(sig_pos_end)))
cat(sprintf("NEG features: %d, sig: %d, endogenous sig: %d\n",
    nrow(limma_neg), sum(limma_neg$padj<0.05 & abs(limma_neg$log2FC)>0.5),
    nrow(sig_neg_end)))
cat(sprintf("Cross-mode shared HMDBs: %d\n", length(hmdb_both)))
cat(sprintf("Li2024 replication: %d/%d agree (%.1f%%)", n_agree, n_total, n_agree/n_total*100))
if (n_disagree > 0) {
    cat(sprintf(", %d disagree:", n_disagree))
    for (i in seq_len(nrow(disagree_rows))) {
        cat(sprintf(" %s", disagree_rows$our_metabolite[i]))
    }
}
cat("\n")

# Annotation flag distribution
cat(sprintf("POS annotation: %s\n", 
    paste(names(table(limma_pos$annotation_confidence)), table(limma_pos$annotation_confidence), 
          sep="=", collapse=", ")))
cat(sprintf("NEG annotation: %s\n",
    paste(names(table(limma_neg$annotation_confidence)), table(limma_neg$annotation_confidence),
          sep="=", collapse=", ")))

sink()
writeLines(capture.output(sessionInfo()), sessionInfo_path)
cat("Done.\n")
