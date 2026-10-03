#!/usr/bin/env Rscript
# ============================================================
#  04_sarcosine_analysis.R  —  CORE: gut-microbial sarcosine metabolism & CRC
#  가설: 분해-연관 속(genus)이 많을수록 분변 sarcosine 낮음; 분해속은 암에서 감소.
#  분석(faithful 5단계 a-d):
#   (a) 분변 sarcosine: Healthy vs CRC
#   (b) 분해속 abundance ~ sarcosine 상관 (가설: 음의 상관)   [핵심]
#   (c) 분해속(+생성속) abundance: Healthy vs CRC
#   (d) within-patient 삼각: 분해속 low <-> sarcosine high <-> Cancer 요약
#  입력: genus_table.tsv(16S), pe_samples.tsv, neg-mode MAF(sarcosine)
#  출력: pe_pipeline/sarcosine/  *.png(300dpi), *.tsv, summary.txt, sessionInfo
#  ⚠️ 한계: sarcosine 이성질체(alanine 등) 분리 미확인·reliability 공란;
#     16S=genus 수준; 샘플 매칭=이름 기반 가정; Clostridium=sensu stricto 1(C.symbiosum 매핑 주의)
#  NOTE: 플롯 텍스트는 ASCII만 사용(grDevices::png가 한글/화살표/그리스문자 미렌더 + ragg는 이 경로에 못 씀).
# ============================================================
suppressMessages(library(ggplot2))
set.seed(42)

BASE <- "."   # analysis-folder root (working dir)
PP   <- file.path(BASE, "16S_데이터_분석/raw_download/pe_pipeline")
PES  <- file.path(BASE, "16S_데이터_분석/raw_download/pe_samples.tsv")
NEG  <- file.path(BASE, "Metabolomics_Data_MetaboLights/Result_files",
                  "m_MTBLS10232_LC-MS_negative_reverse-phase_metabolite_profiling_v2_maf.tsv")
OUT  <- file.path(PP, "sarcosine"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

COLZ <- c(Healthy = "#1B7837", CRC = "#B2182B")
theme_pub <- theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 20, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 14, colour = "grey30", margin = margin(b = 8)),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 15),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.title = element_text(face = "bold"), plot.margin = margin(12, 16, 10, 12))
# ragg(agg_png) fails on this path (colon in 'Liang:Ma'); base grDevices::png works.
savepng <- function(p, f, w = 8, h = 6) ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 300, bg = "white", device = grDevices::png)

## ---------- 데이터 적재 ----------
cnt <- as.matrix(read.delim(file.path(PP, "genus_table.tsv"), row.names = 1, check.names = FALSE))
nonbact <- grepl("Mitochondria|Chloroplast|Eukaryota", rownames(cnt), ignore.case = TRUE)  # 비세균 오염 제거
cat(sprintf("[filter] non-bacterial removed: %d {%s}\n", sum(nonbact), paste(rownames(cnt)[nonbact], collapse = ", ")))
cnt <- cnt[!nonbact, , drop = FALSE]
relab <- sweep(cnt, 2, colSums(cnt), "/")          # genus x sample 상대존재비 (bacteria-only 재계산)
meta <- read.delim(PES, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("run", "layout", "phenotype", "sample_name"))
meta$group <- ifelse(meta$phenotype == "Colorectal Neoplasms", "CRC",
               ifelse(meta$phenotype == "Health", "Healthy", NA))

maf <- read.delim(NEG, check.names = FALSE, stringsAsFactors = FALSE)
srow <- which(maf$metabolite_identification == "Sarcosine")
stopifnot(length(srow) == 1)
scols <- 24:ncol(maf)
sarc <- as.numeric(maf[srow, scols]); names(sarc) <- colnames(maf)[scols]
cat(sprintf("[load] genera=%d, 16S samples=%d, sarcosine samples=%d\n",
            nrow(relab), ncol(relab), length(sarc)))

## ---------- 타깃 속 정의 + 페어드 데이터 구성 ----------
degraders <- c("Lachnospira", "Faecalibacterium", "Coprococcus",
               "Fusicatenibacter", "Roseburia", "Blautia")
producers <- c("Hungatella")   # Clostridium s.s.1 제거: C.symbiosum 아님(strict-sense Clostridium)
eligens   <- "[Eubacterium] eligens group"          # L.eligens 특이(SILVA); 개별 패널에만 추가, degrader_sum 제외
panel_genera <- c(degraders, eligens, producers)    # 개별 표시(04c/prevalence)용
getg <- function(g) if (g %in% rownames(relab)) relab[g, ] else setNames(rep(0, ncol(relab)), colnames(relab))

runs <- colnames(relab)
snN  <- gsub("Control", "CTRL", meta$sample_name[match(runs, meta$run)])   # 이름 정규화
df <- data.frame(run = runs, sample = snN,
                 group = meta$group[match(runs, meta$run)],
                 sarcosine = as.numeric(sarc[snN]), stringsAsFactors = FALSE)
for (g in panel_genera) df[[g]] <- as.numeric(getg(g))
df$degrader_sum <- rowSums(df[, degraders, drop = FALSE])   # 원래 6 분해속만 (eligens/Blautia 합산 영향 없음)
df <- df[!is.na(df$sarcosine) & !is.na(df$group), ]   # 페어드만
df$group <- factor(df$group, levels = c("Healthy", "CRC"))
df$log10_sarcosine <- log10(df$sarcosine)
nH <- sum(df$group == "Healthy"); nC <- sum(df$group == "CRC")
cat(sprintf("[paired] %d samples (Healthy %d / CRC %d)\n", nrow(df), nH, nC))
write.table(df, file.path(OUT, "merged_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

prev <- data.frame(genus = panel_genera,
  role = ifelse(panel_genera %in% producers, "producer", "degrader"),
  prevalence_pct = sapply(panel_genera, function(g) round(100 * mean(getg(g)[df$run] > 0), 1)),
  mean_relab_pct = sapply(panel_genera, function(g) round(100 * mean(getg(g)[df$run]), 3)))
write.table(prev, file.path(OUT, "target_genus_prevalence.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- (a) sarcosine: Healthy vs CRC ----------
pa <- wilcox.test(sarcosine ~ group, data = df)$p.value
medH <- median(df$sarcosine[df$group == "Healthy"]); medC <- median(df$sarcosine[df$group == "CRC"])
ga <- ggplot(df, aes(group, sarcosine, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .75, width = .6) +
  geom_jitter(width = .15, size = .8, alpha = .4) +
  scale_y_log10() + scale_fill_manual(values = COLZ) +
  labs(title = "Fecal sarcosine in colorectal cancer vs. healthy controls",
       subtitle = sprintf("Untargeted LC-MS, negative mode (HMDB0000271); n = %d (Healthy %d, CRC %d)\nWilcoxon p = %.3g; median Healthy %.0f vs CRC %.0f (a.u.)\nCaveat: isomer (alanine / b-alanine) separation not confirmed",
                          nrow(df), nH, nC, pa, medH, medC),
       x = NULL, y = "Sarcosine intensity (a.u., log scale)") +
  theme_pub + theme(legend.position = "none")
savepng(ga, "04a_sarcosine_by_group.png", 8.5, 6.5)

## ---------- (b) 분해속 ~ sarcosine 상관 (핵심) ----------
vars <- c("degrader_sum", degraders)
cor.tab <- do.call(rbind, lapply(vars, function(v) {
  ct <- suppressWarnings(cor.test(df[[v]], df$sarcosine, method = "spearman"))
  data.frame(feature = v, rho = unname(ct$estimate), p = ct$p.value)
}))
cor.tab$q <- p.adjust(cor.tab$p, "BH")
cor.tab <- cor.tab[order(cor.tab$rho), ]
write.table(cor.tab, file.path(OUT, "DEGRADER_sarcosine_correlation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# within-group 상관 (암 상태 교란 점검): pooled 음의 상관이 그룹 내에서도 유지되는가?
wg <- do.call(rbind, lapply(vars, function(v)
  do.call(rbind, lapply(c("Healthy", "CRC"), function(gg) {
    d <- df[df$group == gg, ]
    ct <- suppressWarnings(cor.test(d[[v]], d$sarcosine, method = "spearman"))
    data.frame(feature = v, group = gg, n = nrow(d), rho = unname(ct$estimate), p = ct$p.value)
  }))))
write.table(wg, file.path(OUT, "DEGRADER_sarcosine_correlation_withingroup.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
wgH <- wg$p[wg$feature == "degrader_sum" & wg$group == "Healthy"]
wgC <- wg$p[wg$feature == "degrader_sum" & wg$group == "CRC"]

lab.b <- setNames(sprintf("%s\n(rho = %.2f, p = %.2g)", cor.tab$feature, cor.tab$rho, cor.tab$p),
                  cor.tab$feature)
dlong <- do.call(rbind, lapply(vars, function(v)
  data.frame(feature = factor(lab.b[v], levels = lab.b[vars]),
             relab = df[[v]], log10_sarcosine = df$log10_sarcosine, group = df$group)))
gb <- ggplot(dlong, aes(relab, log10_sarcosine)) +
  geom_point(aes(color = group), size = .9, alpha = .5) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = .6) +
  facet_wrap(~feature, scales = "free_x", ncol = 4) +
  scale_color_manual(values = COLZ) +
  labs(title = "Sarcosine-degrading commensals vs. fecal sarcosine (hypothesis: negative)",
       subtitle = sprintf("%d paired samples; Spearman correlation; line = linear fit. Negative rho supports 'more degraders -> less sarcosine'", nrow(df)),
       x = "Genus relative abundance", y = expression(bold(log[10]~"sarcosine intensity")), color = "Group") +
  theme_pub + theme(strip.text = element_text(size = 12))
savepng(gb, "04b_degrader_vs_sarcosine.png", 14, 7)

## ---------- (c) 분해속/생성속 abundance: Healthy vs CRC ----------
# 패널 순서: clean 분해속 5 -> [Eubacterium] eligens(L.eligens 특이) -> Hungatella(생성) -> Blautia(맨 끝, 이질적)
allg <- c("Lachnospira", "Faecalibacterium", "Coprococcus", "Fusicatenibacter", "Roseburia",
          "[Eubacterium] eligens group", "Hungatella", "Blautia")
grp.tab <- do.call(rbind, lapply(allg, function(g) {
  p <- suppressWarnings(wilcox.test(df[[g]] ~ df$group)$p.value)
  data.frame(genus = g, role = ifelse(g %in% producers, "producer", ifelse(g == "Blautia", "mixed", "degrader")),
             median_Healthy = median(df[[g]][df$group == "Healthy"]),
             median_CRC = median(df[[g]][df$group == "CRC"]),
             log2FC_CRC_vs_H = log2((median(df[[g]][df$group == "CRC"]) + 1e-6) /
                                    (median(df[[g]][df$group == "Healthy"]) + 1e-6)), p = p)
}))
grp.tab$q <- p.adjust(grp.tab$p, "BH")
write.table(grp.tab, file.path(OUT, "DEGRADER_by_group.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

disp <- allg; disp[disp == "[Eubacterium] eligens group"] <- "[Eub.] eligens group"   # 표시 단축
lab.c <- setNames(sprintf("%s\n(p = %.2g)", disp, grp.tab$p[match(allg, grp.tab$genus)]), allg)
clong <- do.call(rbind, lapply(allg, function(g)
  data.frame(genus = factor(lab.c[g], levels = lab.c[allg]),
             relab = df[[g]] + 1e-5, group = df$group)))
gc <- ggplot(clong, aes(group, relab, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .75, width = .6) +
  geom_jitter(width = .14, size = .5, alpha = .3) +
  facet_wrap(~genus, scales = "free_y", ncol = 4) + scale_y_log10() +
  scale_fill_manual(values = COLZ) +
  labs(title = "Sarcosine-degrading and -producing genera: CRC vs. healthy controls",
       subtitle = sprintf("Genus-level 16S (V3-V4); n = %d (Healthy %d, CRC %d); Wilcoxon p per genus (+1e-5 pseudocount for log scale)",
                          nrow(df), nH, nC),
       caption = "Blautia (last) = genus-level heterogeneous, increases in CRC -- not a clean degrader.   [Eub.] eligens group = SILVA label for L. eligens.",
       x = NULL, y = "Relative abundance (log scale)", fill = "Group") +
  theme_pub + theme(legend.position = "none", strip.text = element_text(size = 12),
                    plot.caption = element_text(hjust = 0, size = 11, colour = "grey30"))
savepng(gc, "04c_genera_by_group.png", 14, 7.4)

## ---------- (d) within-patient 삼각 요약 ----------
legA <- sprintf("Sarcosine in CRC: median CRC %.0f vs Healthy %.0f, p=%.3g (CRC %s)",
                medC, medH, pa, ifelse(medC > medH, "higher", "lower"))
rho.sum <- cor.tab$rho[cor.tab$feature == "degrader_sum"]; p.sum <- cor.tab$p[cor.tab$feature == "degrader_sum"]
legB <- sprintf("Degrader_sum vs sarcosine: pooled rho=%.2f (p=%.3g) BUT within-group NS (Healthy p=%.2f, CRC p=%.2f) -> confounded by cancer status",
                rho.sum, p.sum, wgH, wgC)
ds.p <- suppressWarnings(wilcox.test(df$degrader_sum ~ df$group)$p.value)
ds.mH <- median(df$degrader_sum[df$group == "Healthy"]); ds.mC <- median(df$degrader_sum[df$group == "CRC"])
legC <- sprintf("Degrader_sum in CRC: median CRC %.3f vs Healthy %.3f, p=%.3g (%s)",
                ds.mC, ds.mH, ds.p, ifelse(ds.mC < ds.mH, "CRC lower: supports H1", "CRC higher"))
tri <- data.frame(leg = c("(a) Cancer -> Sarcosine", "(b) Degraders -> Sarcosine", "(c) Cancer -> Degraders"),
                  stat = c(legA, legB, legC))
gd <- ggplot(tri, aes(x = 1, y = factor(leg, levels = rev(leg)))) +
  geom_text(aes(label = paste0(leg, "\n", stat)), hjust = 0, x = 0, size = 5.2) +
  xlim(0, 1) +
  labs(title = "The degrader-sarcosine-cancer triangle (within paired CRC cohort)",
       subtitle = sprintf("%d paired samples. Hypothesis holds if: degraders lower in CRC AND degrader-sarcosine negative AND sarcosine higher in CRC",
                          nrow(df))) +
  theme_void(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 14, colour = "grey30"),
        axis.text = element_blank(), plot.margin = margin(14, 16, 14, 16))
savepng(gd, "04d_triangle_summary.png", 15, 4.6)

## ---------- 요약 ----------
sink(file.path(OUT, "summary.txt"))
cat("=== Sarcosine-degrader-cancer 분석 요약 (Liang/Ma CRC 16S+metabolome) ===\n")
cat(sprintf("페어드 샘플: %d (Healthy %d / CRC %d), genus-level 16S\n\n", nrow(df), nH, nC))
cat("[타깃 속 prevalence]\n"); print(prev, row.names = FALSE)
cat("\n[(a) sarcosine Healthy vs CRC]\n", legA, "\n")
cat("\n[(b) 분해속 ~ sarcosine Spearman (pooled)]\n"); print(cor.tab, row.names = FALSE)
cat("\n[(b') within-group 상관 (암 상태 교란 점검)]\n")
print(wg[wg$feature %in% c("degrader_sum", "Faecalibacterium"), ], row.names = FALSE)
cat("\n[(c) 속 abundance Healthy vs CRC]\n"); print(grp.tab, row.names = FALSE)
cat("\n[(d) 삼각 요약]\n", legA, "\n", legB, "\n", legC, "\n")
cat(sprintf("\n[★중요한계] leg(b) 분해속-sarcosine 음의 상관은 그룹 내에서 비유의(Healthy p=%.2f, CRC p=%.2f) -> 대부분 암 상태 교란. 직접 within-individual 관계 미확정.\n", wgH, wgC))
cat("[기타한계] sarcosine 이성질체(alanine 등) 분리 미확인; reliability 공란; 16S=genus 수준;",
    "샘플 매칭=이름 기반 가정; Clostridium=sensu stricto 1(C.symbiosum 매핑 불확실); young/old 혼재.\n")
sink()
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_sarcosine.txt"))
cat("SARCOSINE_ANALYSIS_DONE\n")
