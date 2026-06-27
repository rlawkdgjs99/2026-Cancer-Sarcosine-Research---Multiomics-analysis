#!/usr/bin/env Rscript
# ============================================================
#  03_downstream_analysis.R
#    16S genus-level downstream analysis — Healthy(control) vs CRC(colorectal cancer)
#    입력:
#      - pe_pipeline/genus_table.tsv         (genus x sample, counts; DADA2 산출)
#      - raw_download/pe_samples.tsv          (run, layout, phenotype, sample_name)
#    분석:
#      (1) 알파다양성 (Observed/Shannon/Simpson) + Wilcoxon
#      (2) 베타다양성 (Bray-Curtis PCoA) + PERMANOVA(adonis2)
#      (3) 군집 조성 (top genera, group-mean relative abundance)
#      (4) 차등존재비: ANCOM-BC2 (compositional) + Wilcoxon(rel.abund)+BH 교차검증
#    출력: pe_pipeline/downstream/  *.png(300dpi), *.tsv, summary.txt, sessionInfo
#    재현성: set.seed(42), sessionInfo 저장. plot: PNG, 주제목 강조.
# ============================================================
suppressMessages({ library(phyloseq); library(vegan); library(ggplot2) })
set.seed(42)

## ---------- 경로 / 파라미터 ----------
BASE <- "16S_데이터_분석/raw_download"
PP   <- file.path(BASE, "pe_pipeline")
OUT  <- file.path(PP, "downstream")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

MIN_DEPTH <- 1000     # 샘플 최소 read 수(이하 제외)
PREV      <- 0.10     # genus 최소 출현율(10% 미만 제외)
TOPN      <- 15       # 조성 그림에 표시할 top genus 수
COLZ      <- c(Healthy = "#1B7837", CRC = "#B2182B")   # colorblind-safe

theme_pub <- theme_bw(base_size = 16) +
  theme(plot.title    = element_text(face = "bold", size = 21, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 15, colour = "grey30", margin = margin(b = 8)),
        axis.title    = element_text(face = "bold"),
        strip.text    = element_text(face = "bold", size = 16),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.title  = element_text(face = "bold"),
        plot.margin   = margin(12, 14, 10, 12))
savepng <- function(p, file, w = 8, h = 6)  # ragg fails on this path (colon) -> base png
  ggsave(file.path(OUT, file), p, width = w, height = h, dpi = 300, bg = "white", device = grDevices::png)

## ---------- 데이터 적재 + phyloseq ----------
gt <- read.delim(file.path(PP, "genus_table.tsv"), row.names = 1, check.names = FALSE)
## 비세균 오염(숙주 mitochondria / chloroplast / eukaryota) feature 제거 -> bacteria-only 후 상대존재비 재계산
nonbact <- grepl("Mitochondria|Chloroplast|Eukaryota", rownames(gt), ignore.case = TRUE)
cat(sprintf("[filter] non-bacterial genera removed: %d {%s}\n", sum(nonbact), paste(rownames(gt)[nonbact], collapse = ", ")))
gt <- gt[!nonbact, , drop = FALSE]
meta <- read.delim(file.path(BASE, "pe_samples.tsv"), header = FALSE,
                   col.names = c("run", "layout", "phenotype", "sample_name"),
                   stringsAsFactors = FALSE)
meta$group <- ifelse(meta$phenotype == "Colorectal Neoplasms", "CRC",
               ifelse(meta$phenotype == "Health", "Healthy", NA))
rownames(meta) <- meta$run
common <- intersect(colnames(gt), rownames(meta))
gt   <- gt[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]
meta <- meta[!is.na(meta$group), , drop = FALSE]
gt   <- gt[, rownames(meta), drop = FALSE]
meta$group <- factor(meta$group, levels = c("Healthy", "CRC"))   # Healthy = reference

ps0 <- phyloseq(otu_table(as.matrix(gt), taxa_are_rows = TRUE), sample_data(meta))
cat(sprintf("[load] %d samples, %d genera (pre-filter)\n", nsamples(ps0), ntaxa(ps0)))

## ---------- QC 필터 ----------
ps <- prune_samples(sample_sums(ps0) >= MIN_DEPTH, ps0)
ps <- filter_taxa(ps, function(x) sum(x > 0) >= PREV * nsamples(ps), prune = TRUE)
grp <- factor(sample_data(ps)$group, levels = c("Healthy", "CRC"))
tab <- table(grp); nH <- tab[["Healthy"]]; nC <- tab[["CRC"]]
cat(sprintf("[filter] %d samples (Healthy %d / CRC %d), %d genera; depth>=%d, prev>=%.0f%%\n",
            nsamples(ps), nH, nC, ntaxa(ps), MIN_DEPTH, 100 * PREV))

## ---------- (1) 알파다양성 ----------
ad <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
ad$group <- grp
metrics <- c("Observed", "Shannon", "Simpson")
pa <- sapply(metrics, function(m) wilcox.test(ad[[m]] ~ ad$group)$p.value)
write.table(data.frame(sample = rownames(ad), ad), file.path(OUT, "alpha_diversity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

labs <- setNames(sprintf("%s  (p = %s)", metrics,
                         formatC(pa, format = "g", digits = 2)), metrics)
adl <- do.call(rbind, lapply(metrics, function(m)
  data.frame(metric = factor(labs[m], levels = labs[metrics]),
             value = ad[[m]], group = ad$group)))
g1 <- ggplot(adl, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .75, width = .6) +
  geom_jitter(width = .14, size = .7, alpha = .35) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = COLZ) +
  labs(title = "Fecal microbiome alpha diversity: colorectal cancer vs. healthy controls",
       subtitle = sprintf("Genus-level 16S rRNA (V3-V4), n = %d (Healthy = %d, CRC = %d); Wilcoxon rank-sum test",
                          nsamples(ps), nH, nC),
       x = NULL, y = "Diversity index value") +
  theme_pub + theme(legend.position = "none")
savepng(g1, "01_alpha_diversity.png", 11, 5)

## ---------- (2) 베타다양성 (Bray-Curtis PCoA + PERMANOVA) ----------
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
bc  <- phyloseq::distance(ps.rel, method = "bray")
set.seed(42)
pn  <- adonis2(bc ~ group, data = data.frame(sample_data(ps.rel)), permutations = 999)
write.table(as.data.frame(pn), file.path(OUT, "permanova_braycurtis.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
ord <- ordinate(ps.rel, method = "PCoA", distance = bc)
ev  <- round(ord$values$Relative_eig * 100, 1)
g2 <- plot_ordination(ps.rel, ord, color = "group") +
  stat_ellipse(aes(group = group), linewidth = .9, type = "norm") +
  geom_point(size = 2.3, alpha = .85) +
  scale_color_manual(values = COLZ) +
  labs(title = "Gut microbiota structure: CRC vs. healthy",
       subtitle = sprintf("Bray-Curtis PCoA (genus-level); PERMANOVA R2 = %.3f, p = %.3f",
                          pn$R2[1], pn$`Pr(>F)`[1]),
       x = sprintf("PCoA axis 1 (%.1f%%)", ev[1]),
       y = sprintf("PCoA axis 2 (%.1f%%)", ev[2]), color = "Group") +
  theme_pub
savepng(g2, "02_beta_diversity_pcoa.png", 8, 6.5)

## ---------- (3) 군집 조성 (top genera, group-mean) ----------
relmat <- as.matrix(otu_table(ps.rel))     # genus x sample
relmean <- sort(rowMeans(relmat), decreasing = TRUE)
top <- names(relmean)[seq_len(min(TOPN, length(relmean)))]
comp <- do.call(rbind, lapply(levels(grp), function(g) {
  v   <- rowMeans(relmat[, grp == g, drop = FALSE])
  lab <- ifelse(names(v) %in% top, names(v), "Other")
  a   <- aggregate(v, by = list(genus = lab), FUN = sum)
  data.frame(group = g, genus = a$genus, value = a$x)
}))
comp$genus <- factor(comp$genus, levels = rev(c(top, "Other")))
pal <- setNames(c(colorRampPalette(
  c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d","#666666"))(length(top)),
  "grey80"), c(top, "Other"))
g3 <- ggplot(comp, aes(group, value, fill = genus)) +
  geom_col(width = .62, colour = "white", linewidth = .1) +
  scale_fill_manual(values = pal, breaks = c(top, "Other")) +
  scale_y_continuous(expand = expansion(mult = c(0, .02))) +
  labs(title = sprintf("Top %d fecal genera: CRC vs. healthy", TOPN),
       subtitle = "Genus-level 16S rRNA (V3-V4); bars show group-mean relative abundance",
       x = NULL, y = "Mean relative abundance", fill = "Genus") +
  theme_pub
savepng(g3, "03_composition_top_genera.png", 8.5, 7)

## ---------- (4) 차등존재비 ----------
# (4a) Wilcoxon on relative abundance + BH  (강건; 항상 실행)
wil <- data.frame(genus = rownames(relmat),
                  p = apply(relmat, 1, function(x)
                    tryCatch(wilcox.test(x ~ grp)$p.value, error = function(e) NA_real_)),
                  median_Healthy = apply(relmat[, grp == "Healthy", drop = FALSE], 1, median),
                  median_CRC     = apply(relmat[, grp == "CRC", drop = FALSE], 1, median),
                  stringsAsFactors = FALSE)
wil$q <- p.adjust(wil$p, "BH")
wil$log2FC_CRC_vs_Healthy <- log2((wil$median_CRC + 1e-6) / (wil$median_Healthy + 1e-6))
wil <- wil[order(wil$q), c("genus","median_Healthy","median_CRC","log2FC_CRC_vs_Healthy","p","q")]
write.table(wil, file.path(OUT, "DA_wilcoxon_relabund.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# (4b) ANCOM-BC2 (compositional; best-effort)
da.plot <- NULL; da.method <- NULL
if (requireNamespace("ANCOMBC", quietly = TRUE)) {
  res <- tryCatch({
    set.seed(42)
    ANCOMBC::ancombc2(data = ps, fix_formula = "group", group = "group",
                      p_adj_method = "BH", prv_cut = PREV, lib_cut = 0,
                      struc_zero = TRUE, neg_lb = FALSE, alpha = 0.05,
                      n_cl = 4, verbose = FALSE)
  }, error = function(e) { cat("[ANCOM-BC2] 실패:", conditionMessage(e), "\n"); NULL })
  if (!is.null(res)) {
    o  <- res$res
    lc <- grep("^lfc_group",  names(o), value = TRUE)[1]
    qc <- grep("^q_group",    names(o), value = TRUE)[1]
    dc <- grep("^diff_group", names(o), value = TRUE)[1]
    anc <- data.frame(genus = o$taxon, lfc = o[[lc]], q = o[[qc]], diff = o[[dc]])
    anc <- anc[order(anc$q), ]
    write.table(anc, file.path(OUT, "DA_ancombc2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    sig <- anc[which(anc$diff & anc$q < 0.05), ]
    if (nrow(sig) > 0) {
      da.plot <- data.frame(genus = sig$genus, lfc = sig$lfc, q = sig$q)
      da.method <- "ANCOM-BC2"
    }
  }
}
if (is.null(da.plot)) {                       # fallback: Wilcoxon
  sw <- wil[which(wil$q < 0.05), ]
  if (nrow(sw) > 0) {
    da.plot <- data.frame(genus = sw$genus, lfc = sw$log2FC_CRC_vs_Healthy, q = sw$q)
    da.method <- "Wilcoxon (relative abundance) + BH"
  }
}
if (!is.null(da.plot) && nrow(da.plot) > 0) {
  ntot <- nrow(da.plot)                                                    # 전체 유의 genus 수
  d <- da.plot[order(da.plot$lfc), ]
  if (nrow(d) > 30) d <- d[order(abs(d$lfc), decreasing = TRUE)[1:30], ]   # 너무 많으면 효과크기 top30
  d <- d[order(d$lfc), ]
  d$genus <- factor(d$genus, levels = d$genus)
  d$dir <- ifelse(d$lfc > 0, "Enriched in CRC", "Enriched in Healthy")
  subt <- if (ntot > nrow(d))
    sprintf("%s; q < 0.05; top %d of %d significant (by |effect|); LFC > 0 = higher in CRC", da.method, nrow(d), ntot)
  else sprintf("%s; q < 0.05; %d significant genera; LFC > 0 = higher in CRC", da.method, nrow(d))
  g4 <- ggplot(d, aes(lfc, genus, fill = dir)) +
    geom_col() + geom_vline(xintercept = 0, linewidth = .3) +
    scale_fill_manual(values = c("Enriched in CRC" = "#B2182B", "Enriched in Healthy" = "#1B7837")) +
    labs(title = "Differentially abundant fecal genera: CRC vs. healthy controls",
         subtitle = subt,
         x = "Log fold change (CRC vs. Healthy)", y = NULL, fill = NULL) +
    theme_pub + theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14))
  savepng(g4, "04_differential_abundance.png", 11, max(4, 0.32 * nrow(d) + 2))
  nsig <- ntot                                  # 요약엔 그림 표시수(top30)가 아니라 실제 유의 총수
} else { nsig <- 0; da.method <- "none significant" }

## ---------- 요약 + sessionInfo ----------
dHe <- median(sample_sums(ps)[grp == "Healthy"]); dCr <- median(sample_sums(ps)[grp == "CRC"])
sink(file.path(OUT, "summary.txt"))
cat("=== 16S downstream summary — Healthy vs CRC ===\n")
cat(sprintf("Samples (post-QC): %d  (Healthy %d / CRC %d)\n", nsamples(ps), nH, nC))
cat(sprintf("Genera (post-filter, prev>=%.0f%%): %d\n", 100*PREV, ntaxa(ps)))
cat(sprintf("Median read depth: Healthy %.0f / CRC %.0f\n", dHe, dCr))
cat(sprintf("Alpha (Wilcoxon p): Observed %.3g | Shannon %.3g | Simpson %.3g\n",
            pa["Observed"], pa["Shannon"], pa["Simpson"]))
cat(sprintf("Beta PERMANOVA (Bray-Curtis): R2 = %.4f, p = %.3f\n", pn$R2[1], pn$`Pr(>F)`[1]))
cat(sprintf("Differential abundance: method=%s, significant genera=%d (q<0.05)\n", da.method, nsig))
cat("\nOutputs: 01_alpha_diversity.png, 02_beta_diversity_pcoa.png,\n",
    "         03_composition_top_genera.png, 04_differential_abundance.png,\n",
    "         alpha_diversity.tsv, permanova_braycurtis.tsv, DA_wilcoxon_relabund.tsv, DA_ancombc2.tsv\n")
sink()
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_downstream.txt"))
cat("DOWNSTREAM_DONE\n")
