# Dissertation-quality plot regeneration
# Consistent aesthetics, journal-style p-values, 300 dpi, large readable fonts.
#
# 2026-04-29 - generated for dissertation use.
# Inputs (all already exist on disk):
#   analysis/intensity_matrices.rds
#   analysis/diff_results_POS.tsv, diff_results_NEG.tsv
#   metadata_merged_619samples.tsv
#   MTBLS10232_inventory/m_*_maf.tsv

sessionInfo_path <- file.path("R_scripts", paste0("sessionInfo_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

suppressMessages({
    library(ggplot2)
    library(limma)
})

# Output routing (2026 reorg): headline figures -> Main, all others -> Supplementary
dir_main <- "analysis/main"
dir_supp <- "analysis/supplementary"
dir.create(dir_main, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_supp, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helpers - journal-style p-value formatting
# ============================================================

# Plain text p-value (for subtitles / captions)
format_p_text <- function(p, prefix = "p") {
    if (is.na(p)) return(paste0(prefix, " = NA"))
    if (p < 0.001) return(paste0(prefix, " < 0.001"))
    if (p < 0.01)  return(sprintf("%s = %.3f", prefix, p))
    if (p < 0.05)  return(sprintf("%s = %.3f", prefix, p))
    if (p < 0.1)   return(sprintf("%s = %.2f", prefix, p))
    return(sprintf("%s = %.2f", prefix, p))
}

# Plotmath p-value (for in-plot annotation)
format_p_math <- function(p) {
    if (is.na(p)) return("italic(p)*' = NA'")
    if (p < 0.001) return("italic(p) < 0.001")
    return(sprintf("italic(p) == %.3f", round(p, 3)))
}

format_lfc <- function(lfc) sprintf("%+.2f", lfc)

# ============================================================
# Theme for dissertation - large, clean, print-ready
# ============================================================
theme_dissertation <- function(base_size = 16) {
    theme_classic(base_size = base_size) +
    theme(
        plot.title       = element_text(face = "bold", size = base_size + 2,
                                        margin = margin(b = 6)),
        plot.subtitle    = element_text(size = base_size - 2, color = "grey25",
                                        margin = margin(b = 10)),
        plot.caption     = element_text(size = base_size - 4, face = "italic",
                                        color = "grey40", hjust = 1,
                                        margin = margin(t = 8)),
        axis.title.x     = element_text(face = "bold", size = base_size,
                                        margin = margin(t = 8)),
        axis.title.y     = element_text(face = "bold", size = base_size,
                                        margin = margin(r = 8)),
        axis.text        = element_text(size = base_size - 1, color = "grey15"),
        axis.line        = element_line(color = "grey20", size = 0.5),
        axis.ticks       = element_line(color = "grey20", size = 0.4),
        legend.title     = element_text(face = "bold", size = base_size - 2),
        legend.text      = element_text(size = base_size - 2),
        legend.position  = "bottom",
        legend.box.margin= margin(t = 6),
        panel.grid.major.y = element_line(color = "grey92", size = 0.3),
        panel.grid.minor   = element_blank(),
        plot.margin      = margin(15, 15, 12, 12)
    )
}

# ============================================================
# Color palettes
# ============================================================
pal_crc      <- c("CTRL" = "#1B7837", "CRC" = "#B2182B")
pal_qc       <- c("Patient" = "#C84B4B", "QC" = "#3CB371")
pal_cohort   <- c("older_CTRL"   = "#7AAED9",
                  "older_CRC"    = "#D44848",
                  "younger_CTRL" = "#4878A6",
                  "younger_CRC"  = "#9C3030",
                  "QC"           = "#3CB371")
pal_location <- c("Control"         = "#4878A6",
                  "RCC"             = "#C84B4B",
                  "LCC"             = "#E0922F",
                  "Rectum"          = "#7B5BA6",
                  "CRC_unspecified" = "#888888",
                  "QC"              = "#3CB371")
pal_anno <- c(
    "endogenous_likely"      = "grey55",
    "plant_exogenous"        = "#E0922F",
    "drug_related"           = "#D55E00",
    "contaminant_xenobiotic" = "#C30000"
)

# ============================================================
# Load common data
# ============================================================
cat("Loading data ...\n")
dat       <- readRDS("analysis/intensity_matrices.rds")
metadata  <- read.table("metadata_merged_619samples.tsv", header = TRUE, sep = "\t",
                        check.names = FALSE, stringsAsFactors = FALSE)
diff_pos <- read.table(file.path(dir_main, "diff_results_POS.tsv"), header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, comment.char = "", quote = "\"", fill = TRUE)
diff_neg <- read.table(file.path(dir_main, "diff_results_NEG.tsv"), header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, comment.char = "", quote = "\"", fill = TRUE)

# Build QN-corrected matrices once (correct axis: features × samples)
build_qn <- function(intensity_mat) {
    imp <- intensity_mat
    imp[imp == 0] <- min(intensity_mat[intensity_mat > 0]) / 2
    normalizeBetweenArrays(log2(imp), method = "quantile")
}
qn_pos <- build_qn(dat$intensity_pos)
qn_neg <- build_qn(dat$intensity_neg)

# ============================================================
# 1. Headline boxplot - Sarcosine CRC vs CTRL (NEG)
# ============================================================
cat("Plot 1: Sarcosine CRC vs CTRL ...\n")
{
    sarc_idx <- which(dat$feature_meta_neg$metabolite_identification == "Sarcosine" &
                       dat$feature_meta_neg$database_identifier == "HMDB0000271")
    sarc_v <- qn_neg[sarc_idx, ]
    sc_meta <- metadata[match(colnames(qn_neg), metadata$MTBLS_Sample_Name), ]
    df <- data.frame(
        log2_intensity = sarc_v,
        group = ifelse(grepl("CRC",  sc_meta$cohort_subset), "CRC",
                ifelse(grepl("CTRL", sc_meta$cohort_subset), "CTRL", NA)),
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$group), ]
    df$group <- factor(df$group, levels = c("CTRL", "CRC"))

    crc_v  <- df$log2_intensity[df$group == "CRC"]
    ctrl_v <- df$log2_intensity[df$group == "CTRL"]
    wp     <- wilcox.test(crc_v, ctrl_v)$p.value
    limma_lfc <- diff_neg$log2FC[diff_neg$method == "limma" & diff_neg$metabolite == "Sarcosine"][1]
    limma_padj <- diff_neg$padj[diff_neg$method == "limma" & diff_neg$metabolite == "Sarcosine"][1]

    y_max <- max(df$log2_intensity)
    bracket_y <- y_max + 0.4
    label_y   <- y_max + 0.7

    # Two-line subtitle so it fits comfortably on a 9-inch wide canvas
    sub_text <- sprintf("CTRL n = %d  |  CRC n = %d  |  NEG mode\nlimma log2FC = %s (age-adjusted), %s",
                        length(ctrl_v), length(crc_v),
                        format_lfc(limma_lfc), format_p_text(limma_padj, "FDR-adjusted p"))

    p <- ggplot(df, aes(x = group, y = log2_intensity, fill = group)) +
        geom_jitter(width = 0.18, alpha = 0.30, size = 0.7, color = "grey20") +
        geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55,
                     color = "grey15", size = 0.5) +
        scale_fill_manual(values = pal_crc, guide = "none") +
        scale_y_continuous(breaks = seq(8, 16, by = 1)) +
        # Significance bracket
        annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, color = "grey15") +
        annotate("segment", x = 1, xend = 1, y = bracket_y - 0.1, yend = bracket_y, color = "grey15") +
        annotate("segment", x = 2, xend = 2, y = bracket_y - 0.1, yend = bracket_y, color = "grey15") +
        annotate("text", x = 1.5, y = label_y,
                 label = sprintf("Wilcoxon %s", format_p_text(wp)),
                 size = 7, fontface = "italic") +
        labs(
            title    = "Sarcosine is elevated in CRC fecal samples",
            subtitle = sub_text,
            caption  = "Putative annotation (MSI Level >= 2); isobaric with D-Alanine (m/z 88.04, RT 6.65 min). MS/MS not deposited.",
            x = NULL,
            y = expression(bold(log[2]*"(QN-normalized intensity)"))
        ) +
        theme_dissertation(17)
    ggsave(file.path(dir_main, "sarcosine_CRC_vs_CTRL.png"), p, width = 9, height = 7.5, dpi = 300)
}

# ============================================================
# 2. Sarcosine by cohort (5-group) and by location
# ============================================================
cat("Plot 2: Sarcosine by cohort and by location ...\n")
{
    sarc_idx <- which(dat$feature_meta_neg$metabolite_identification == "Sarcosine" &
                       dat$feature_meta_neg$database_identifier == "HMDB0000271")
    sarc_v <- qn_neg[sarc_idx, ]
    sc_meta <- metadata[match(colnames(qn_neg), metadata$MTBLS_Sample_Name), ]
    df_full <- data.frame(
        log2_intensity = sarc_v,
        cohort = sc_meta$cohort_subset,
        location = sc_meta$rcc_lcc,
        stringsAsFactors = FALSE
    )

    # By cohort (excluding QC)
    df_c <- df_full[df_full$cohort != "QC", ]
    df_c$cohort <- factor(df_c$cohort,
        levels = c("older_CTRL","younger_CTRL","older_CRC","younger_CRC"))
    p2a <- ggplot(df_c, aes(x = cohort, y = log2_intensity, fill = cohort)) +
        geom_jitter(width = 0.18, alpha = 0.30, size = 0.6, color = "grey20") +
        geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6, color = "grey15", size = 0.5) +
        scale_fill_manual(values = pal_cohort, guide = "none") +
        labs(title    = "Sarcosine by age-cohort and disease group",
             subtitle = "MTBLS10232 fecal NEG-mode metabolomics (n = 557 patients)",
             x = NULL,
             y = expression(bold(log[2]*"(QN-normalized intensity)")),
             caption = "Putative Sarcosine annotation; isobaric with D-Alanine.") +
        theme_dissertation(16) +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))
    ggsave(file.path(dir_supp, "sarcosine_by_cohort.png"), p2a, width = 8, height = 6.5, dpi = 300)

    # By tumor location (CRC subgroups + Control + QC excluded)
    df_l <- df_full[df_full$location %in% c("Control", "RCC", "LCC", "Rectum", "CRC_unspecified"), ]
    df_l$location <- factor(df_l$location,
        levels = c("Control", "RCC", "LCC", "Rectum", "CRC_unspecified"))
    p2b <- ggplot(df_l, aes(x = location, y = log2_intensity, fill = location)) +
        geom_jitter(width = 0.18, alpha = 0.30, size = 0.6, color = "grey20") +
        geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6, color = "grey15", size = 0.5) +
        scale_fill_manual(values = pal_location, guide = "none") +
        labs(title    = "Sarcosine by tumor location",
             subtitle = "RCC = right colon, LCC = left colon",
             x = NULL,
             y = expression(bold(log[2]*"(QN-normalized intensity)")),
             caption = "Putative Sarcosine annotation; isobaric with D-Alanine.") +
        theme_dissertation(16) +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))
    ggsave(file.path(dir_supp, "sarcosine_by_location.png"), p2b, width = 8, height = 6.5, dpi = 300)
}

# ============================================================
# 3. Glycine, Betaine - pathway-neighbor boxplots, two polarities each
# ============================================================
cat("Plot 3: Glycine and Betaine pathway-neighbor boxplots ...\n")
make_neighbor_plot <- function(metabolite_name, polarity, group_var = c("cohort", "location")) {
    group_var <- match.arg(group_var)
    intensity <- if (polarity == "POS") dat$intensity_pos else dat$intensity_neg
    feat_meta <- if (polarity == "POS") dat$feature_meta_pos else dat$feature_meta_neg
    samples   <- if (polarity == "POS") dat$sample_cols_pos else dat$sample_cols_neg

    feat_idx <- which(feat_meta$metabolite_identification == metabolite_name)
    if (length(feat_idx) == 0) {
        cat(sprintf("  [skip] '%s' not detected in %s mode (exact match); no %s plot written.\n",
                    metabolite_name, polarity, group_var))
        return(invisible(NULL))
    }
    raw_v <- intensity[feat_idx[1], ]
    min_nz <- min(intensity[intensity > 0])
    log_v <- log2(pmax(raw_v, min_nz/2))

    sc_meta <- metadata[match(samples, metadata$MTBLS_Sample_Name), ]
    df <- data.frame(
        log2_intensity = log_v,
        cohort = sc_meta$cohort_subset,
        location = sc_meta$rcc_lcc,
        stringsAsFactors = FALSE
    )

    if (group_var == "cohort") {
        df <- df[df$cohort != "QC" & !is.na(df$cohort), ]
        df$grp <- factor(df$cohort,
            levels = c("older_CTRL","younger_CTRL","older_CRC","younger_CRC"))
        pal <- pal_cohort
        suffix <- "cohort"
        xlab <- NULL
    } else {
        df <- df[df$location %in% c("Control", "RCC", "LCC", "Rectum", "CRC_unspecified"), ]
        df$grp <- factor(df$location,
            levels = c("Control","RCC","LCC","Rectum","CRC_unspecified"))
        pal <- pal_location
        suffix <- "location"
        xlab <- NULL
    }

    p <- ggplot(df, aes(x = grp, y = log2_intensity, fill = grp)) +
        geom_jitter(width = 0.18, alpha = 0.30, size = 0.6, color = "grey20") +
        geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6, color = "grey15", size = 0.5) +
        scale_fill_manual(values = pal, guide = "none") +
        labs(title    = sprintf("%s - %s mode", metabolite_name, polarity),
             subtitle = sprintf("Stratification: %s", suffix),
             x = xlab,
             y = expression(bold(log[2]*"(intensity, raw)"))) +
        theme_dissertation(16) +
        theme(axis.text.x = element_text(angle = 25, hjust = 1))

    fname <- file.path(dir_supp, sprintf("%s_%s_%s.png",
                     tolower(metabolite_name), suffix, polarity))
    ggsave(fname, p, width = 7.5, height = 6, dpi = 300)
    invisible(NULL)
}

for (mn in c("Glycine", "Betaine")) {
    for (pol in c("POS", "NEG")) {
        for (gv in c("cohort", "location")) {
            make_neighbor_plot(mn, pol, gv)
        }
    }
}

# ============================================================
# 4. Volcano plots - POS and NEG (annotation-colored, large fonts)
# ============================================================
cat("Plot 4: Volcanoes ...\n")
make_volcano <- function(diff_df, pol, highlight = NULL) {
    v <- diff_df[diff_df$method == "limma", ]
    v$neg_log10_padj <- -log10(v$padj)
    v$neg_log10_padj[!is.finite(v$neg_log10_padj)] <- max(v$neg_log10_padj[is.finite(v$neg_log10_padj)], na.rm=TRUE)
    v$significant <- with(v, padj < 0.05 & abs(log2FC) > 0.5)

    # Grey out ONLY points below BOTH cut-offs (not significant by FDR AND |log2FC| <= 0.5).
    # Any point clearing EITHER cut-off keeps its direction colour, so FDR-significant features
    # with a modest fold change (e.g. Sarcosine) stay coloured.
    # NB: a point is a metabolite feature, not a sample; colour encodes direction of enrichment.
    v$below_both <- (v$padj >= 0.05) & (abs(v$log2FC) <= 0.5)
    v$plot_group <- ifelse(v$below_both, "Below both cut-offs",
                    ifelse(v$log2FC > 0, "Higher in Cancer", "Higher in Healthy"))
    v$plot_group <- factor(v$plot_group,
                           levels = c("Higher in Healthy", "Higher in Cancer", "Below both cut-offs"))
    pal_dir <- c("Higher in Healthy"   = "#1B7837",   # CTRL green
                 "Higher in Cancer"    = "#B2182B",   # CRC red
                 "Below both cut-offs" = "grey78")

    # Significant counts by direction (all features; the annotation filter is no longer applied here)
    n_up_cancer  <- sum(v$significant & v$log2FC > 0)
    n_up_healthy <- sum(v$significant & v$log2FC < 0)

    # Plot the greyed (below-both) points first, coloured points on top for visibility
    v <- v[order(v$below_both, decreasing = TRUE), ]

    p <- ggplot(v, aes(x = log2FC, y = neg_log10_padj, color = plot_group)) +
        geom_point(size = 1.6, alpha = 0.6) +
        scale_color_manual(values = pal_dir, name = NULL) +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey30", size = 0.55) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30", size = 0.55) +
        # Cut-off labels (values also stated in the subtitle)
        annotate("text", x = -Inf, y = -log10(0.05), label = "FDR = 0.05",
                 hjust = -0.08, vjust = -0.6, size = 5, colour = "grey25", fontface = "italic") +
        annotate("text", x = 0.5, y = Inf, label = "|log2FC| = 0.5",
                 hjust = -0.08, vjust = 1.5, size = 5, colour = "grey25", fontface = "italic") +
        labs(
            title    = sprintf("%s mode - Differential fecal metabolites: Cancer vs Healthy", pol),
            subtitle = sprintf("Significant: %d higher in Cancer, %d higher in Healthy  (FDR < 0.05 and |log2FC| > 0.5)",
                               n_up_cancer, n_up_healthy),
            x = expression(bold("log"[2]~"fold change (CRC / CTRL, age-adjusted)")),
            y = expression(bold("-log"[10]~"(FDR-adjusted "*italic(p)*")"))
        ) +
        guides(color = guide_legend(override.aes = list(size = 5, alpha = 1), nrow = 1)) +
        theme_dissertation(19)

    # Optionally designate a metabolite of interest (e.g. Sarcosine) with a ringed point + leader label.
    if (!is.null(highlight)) {
        hl <- v[v$metabolite == highlight, , drop = FALSE]
        if (nrow(hl) > 0) {
            hl <- hl[which.min(hl$padj), , drop = FALSE]   # most-significant row if duplicated names
            p <- p +
                ggrepel::geom_text_repel(
                    data = hl, aes(x = log2FC, y = neg_log10_padj, label = metabolite),
                    inherit.aes = FALSE, size = 7, fontface = "bold", colour = "black",
                    nudge_x = -0.30, nudge_y = 4.2, box.padding = 0.8, point.padding = 0.1,
                    segment.colour = "black", segment.size = 0.5,
                    min.segment.length = 0, max.overlaps = Inf)
        } else {
            message(sprintf("  [volcano] highlight '%s' not found in %s; skipping label.", highlight, pol))
        }
    }
    p
}
ggsave(file.path(dir_main, "volcano_POS.png"), make_volcano(diff_pos, "POS"), width = 11, height = 9, dpi = 300)
ggsave(file.path(dir_main, "volcano_NEG.png"), make_volcano(diff_neg, "NEG", highlight = "Sarcosine"), width = 11, height = 9, dpi = 300)

# ============================================================
# 5. PCA - sample_type (Patient vs QC), POS and NEG
# ============================================================
cat("Plot 5: PCA ...\n")
make_pca <- function(intensity_mat, feature_meta, samples, pol) {
    log_mat <- log2(pmax(intensity_mat, min(intensity_mat[intensity_mat>0])/2))
    var_feat <- apply(log_mat, 1, var)
    log_mat  <- log_mat[var_feat > 0, ]
    pca <- prcomp(t(log_mat), center = TRUE, scale. = TRUE)
    pc_var <- summary(pca)$importance[2, 1:2] * 100
    sc_meta <- metadata[match(samples, metadata$MTBLS_Sample_Name), ]
    df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                     sample_type = sc_meta$sample_type)

    p <- ggplot(df, aes(x = PC1, y = PC2, color = sample_type, fill = sample_type)) +
        geom_point(size = 2, alpha = 0.55, shape = 21, color = "grey15", stroke = 0.2) +
        scale_fill_manual(values = pal_qc, name = "Sample type") +
        labs(
            title = sprintf("%s mode - PCA (raw log2, no QN)", pol),
            subtitle = "QCs (green) cluster tightly in centre - confirms instrument stability",
            x = sprintf("PC1 (%.1f%%)", pc_var[1]),
            y = sprintf("PC2 (%.1f%%)", pc_var[2])
        ) +
        guides(fill = guide_legend(override.aes = list(size = 4, alpha = 0.9))) +
        theme_dissertation(16)
    p
}
ggsave(file.path(dir_supp, "pca_POS_sample_type.png"),
       make_pca(dat$intensity_pos, dat$feature_meta_pos, dat$sample_cols_pos, "POS"),
       width = 8, height = 7, dpi = 300)
ggsave(file.path(dir_supp, "pca_NEG_sample_type.png"),
       make_pca(dat$intensity_neg, dat$feature_meta_neg, dat$sample_cols_neg, "NEG"),
       width = 8, height = 7, dpi = 300)

# ============================================================
# 6. TIS drift
# ============================================================
cat("Plot 6: TIS drift ...\n")
make_tis <- function(intensity_mat, samples, pol) {
    tis <- colSums(intensity_mat, na.rm = TRUE)
    sc_meta <- metadata[match(samples, metadata$MTBLS_Sample_Name), ]
    df <- data.frame(order = seq_along(tis), tis = tis,
                     sample_type = sc_meta$sample_type)
    p <- ggplot(df, aes(x = order, y = tis, color = sample_type)) +
        geom_point(size = 1.2, alpha = 0.55) +
        geom_smooth(aes(group = 1), method = "loess", se = TRUE,
                    color = "grey15", fill = "grey80", size = 0.7) +
        scale_color_manual(values = pal_qc, name = "Sample type") +
        scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
        labs(
            title    = sprintf("%s mode - Total Ion Sum vs Sample Order", pol),
            subtitle = "QC samples (green) are tight around the running mean - minimal acquisition drift",
            x = "Sample order in assay file", y = "Total ion sum"
        ) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.9))) +
        theme_dissertation(16)
    p
}
ggsave(file.path(dir_supp, "tis_drift_POS.png"),
       make_tis(dat$intensity_pos, dat$sample_cols_pos, "POS"),
       width = 10, height = 6, dpi = 300)
ggsave(file.path(dir_supp, "tis_drift_NEG.png"),
       make_tis(dat$intensity_neg, dat$sample_cols_neg, "NEG"),
       width = 10, height = 6, dpi = 300)

# ============================================================
# 7. Cohort heterogeneity scatter (full vs older / younger)
# ============================================================
cat("Plot 7: Cohort heterogeneity scatter ...\n")
make_cohort_scatter <- function(intensity_mat, feature_meta, samples, pol, age_filter, age_label) {
    qn <- build_qn(intensity_mat)
    sc_meta <- metadata[match(samples, metadata$MTBLS_Sample_Name), ]
    is_pat <- sc_meta$sample_type == "Patient"
    grp <- ifelse(grepl("CRC", sc_meta$cohort_subset[is_pat]), "CRC", "CTRL")
    age <- ifelse(grepl("^older_", sc_meta$cohort_subset[is_pat]), "older", "younger")
    qn_pat <- qn[, is_pat]

    fc_full <- rowMeans(qn_pat[, grp=="CRC"]) - rowMeans(qn_pat[, grp=="CTRL"])

    keep <- (age == age_filter)
    fc_sub <- rowMeans(qn_pat[, grp=="CRC"  & keep]) - rowMeans(qn_pat[, grp=="CTRL" & keep])

    df <- data.frame(full = fc_full, sub = fc_sub)
    r  <- cor(df$full, df$sub, method = "spearman")

    p <- ggplot(df, aes(x = full, y = sub)) +
        geom_point(size = 0.8, alpha = 0.4, color = "grey25") +
        geom_abline(slope = 1, intercept = 0, color = "#C84B4B", linetype = "dashed", size = 0.7) +
        labs(
            title    = sprintf("%s mode: effect-size concordance", pol),
            subtitle = sprintf("Full cohort vs %s   |   Spearman r = %.3f   |   n = %d features",
                               age_label, r, nrow(df)),
            x = expression(bold("log"[2]~"FC (full cohort)")),
            y = bquote(bold("log"[2]~"FC ("*.(age_label)*")"))
        ) +
        theme_dissertation(16)
    p
}
for (pol in c("POS", "NEG")) {
    intensity <- if (pol == "POS") dat$intensity_pos else dat$intensity_neg
    feat_meta <- if (pol == "POS") dat$feature_meta_pos else dat$feature_meta_neg
    samples   <- if (pol == "POS") dat$sample_cols_pos else dat$sample_cols_neg
    ggsave(file.path(dir_supp, sprintf("cohort_fc_full_vs_older_%s.png", pol)),
           make_cohort_scatter(intensity, feat_meta, samples, pol, "older",   "older subset"),
           width = 8.5, height = 7.5, dpi = 300)
    ggsave(file.path(dir_supp, sprintf("cohort_fc_full_vs_younger_%s.png", pol)),
           make_cohort_scatter(intensity, feat_meta, samples, pol, "younger", "younger subset"),
           width = 8.5, height = 7.5, dpi = 300)
}

cat("All dissertation plots regenerated to analysis/*.png at 300 dpi.\n")

writeLines(capture.output(sessionInfo()), sessionInfo_path)
