#!/usr/bin/env Rscript
# ============================================================================
#  04b_sarcosine_analysis_eligens.R
#  --------------------------------------------------------------------------
#  변형본(variant) — 원본 04_sarcosine_analysis.R 은 수정하지 않는다.
#
#  왜 만들었나
#    원본의 degrader_sum 은
#        Lachnospira, Faecalibacterium, Coprococcus, Fusicatenibacter,
#        Roseburia, Blautia          (eligens 제외 / Blautia 포함)
#    으로 계산되는데, 원고 Methods 와 Sup Fig 15a legend 는
#        Lachnospira, Faecalibacterium, Coprococcus, Fusicatenibacter,
#        Roseburia, [Eubacterium] eligens group   (Blautia 는 "별도 취급,
#        degradation-associated set 에 포함하지 않음")
#    이라고 서술한다. 서술과 계산이 어긋난다.
#
#    Blautia 는 CRC 에서 오히려 증가하는 속(log2FC +0.88)이므로 "분해연관 합"에
#    넣는 것은 정의상 모순이고, eligens 는 L. eligens 의 SILVA 대응 taxon 이라
#    Fig 2/3 의 core degradation 종과 직접 이어진다. 따라서 이 변형본은
#    원고 서술을 정본으로 보고 그대로 계산한다.
#
#  원본과 다른 점은 정확히 세 가지뿐
#    (1) BASE 경로를 현재 폴더 위치로 수정 (원본은 옛 폴더명이 하드코딩되어
#        지금은 실행되지 않는다)
#    (2) degraders 벡터: Blautia -> [Eubacterium] eligens group
#    (3) 출력 폴더: sarcosine/ -> sarcosine_eligens/   (원본 산출물 무손상)
#  그 외 필터·통계·플롯 설정은 원본과 동일하게 유지하여 두 결과가 직접
#  비교 가능하도록 했다.
#
#  검증
#    - 5개 공통 속의 합이 원본 merged_per_sample.tsv 와 일치하는지 assert
#    - 샘플 수 308 (Healthy 121 / CRC 187) assert
#    - 원본 값과의 대조표를 콘솔과 summary.txt 에 출력
#
#  NOTE: 이 경로는 'Liang:Ma' 의 콜론 + 한글 상위폴더 때문에 ragg(agg_png) 가
#        실패한다. base grDevices::png 만 동작하므로 savepng 이 이를 고정한다.
#        (원본 스크립트에도 동일한 주석이 있다.)
# ============================================================================
suppressMessages(library(ggplot2))
set.seed(42)

BASE <- "."   # analysis-folder root (working dir)
PP   <- file.path(BASE, "16S_데이터_분석/raw_download/pe_pipeline")
PES  <- file.path(BASE, "16S_데이터_분석/raw_download/pe_samples.tsv")
NEG  <- file.path(BASE, "Metabolomics_Data_MetaboLights/Result_files",
                  "m_MTBLS10232_LC-MS_negative_reverse-phase_metabolite_profiling_v2_maf.tsv")
OUT  <- file.path(PP, "sarcosine_eligens"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ORIG <- file.path(PP, "sarcosine")   # 원본 산출물 (읽기 전용, 대조용)

COLZ <- c(Healthy = "#1B7837", CRC = "#B2182B")
theme_pub <- theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 20, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 14, colour = "grey30", margin = margin(b = 8)),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 15),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.title = element_text(face = "bold"), plot.margin = margin(12, 16, 10, 12))
savepng <- function(p, f, w = 8, h = 6) {
  ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 300, bg = "white",
         device = grDevices::png)
  if (!file.exists(file.path(OUT, f)) || file.size(file.path(OUT, f)) < 10000)
    stop("PNG 저장 실패: ", f)          # ragg 처럼 조용히 실패하는 것을 막는다
}

## ---------- 데이터 적재 (원본과 동일) ----------
cnt <- as.matrix(read.delim(file.path(PP, "genus_table.tsv"), row.names = 1, check.names = FALSE))
nonbact <- grepl("Mitochondria|Chloroplast|Eukaryota", rownames(cnt), ignore.case = TRUE)
cat(sprintf("[filter] non-bacterial removed: %d {%s}\n", sum(nonbact), paste(rownames(cnt)[nonbact], collapse = ", ")))
cnt <- cnt[!nonbact, , drop = FALSE]
relab <- sweep(cnt, 2, colSums(cnt), "/")
meta <- read.delim(PES, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("run", "layout", "phenotype", "sample_name"))
meta$group <- ifelse(meta$phenotype == "Colorectal Neoplasms", "CRC",
               ifelse(meta$phenotype == "Health", "Healthy", NA))

maf <- read.delim(NEG, check.names = FALSE, stringsAsFactors = FALSE)
srow <- which(maf$metabolite_identification == "Sarcosine")
stopifnot(length(srow) == 1)
scols <- 24:ncol(maf)
stopifnot(!grepl("^smallmolecule|database|^keggdatabase|^hmdbdatabase$", colnames(maf)[24]))  # 24번이 샘플 컬럼인지 확인
sarc <- as.numeric(maf[srow, scols]); names(sarc) <- colnames(maf)[scols]
cat(sprintf("[load] genera=%d, 16S samples=%d, sarcosine samples=%d\n",
            nrow(relab), ncol(relab), length(sarc)))

## ---------- 타깃 속 정의  ★★ 원본과 다른 유일한 지점 ★★ ----------
BASE5     <- c("Lachnospira", "Faecalibacterium", "Coprococcus", "Fusicatenibacter", "Roseburia")
eligens   <- "[Eubacterium] eligens group"
degraders <- c(BASE5, eligens)      # ← 원본은 c(BASE5, "Blautia")
producers <- c("Hungatella")
blautia   <- "Blautia"              # 개별 패널에만 표시(합계 제외)
panel_genera <- c(degraders, producers, blautia)
stopifnot(all(panel_genera %in% rownames(relab)))   # 조용한 0 벡터 방지
getg <- function(g) relab[g, ]

runs <- colnames(relab)
snN  <- gsub("Control", "CTRL", meta$sample_name[match(runs, meta$run)])
df <- data.frame(run = runs, sample = snN,
                 group = meta$group[match(runs, meta$run)],
                 sarcosine = as.numeric(sarc[snN]), stringsAsFactors = FALSE)
for (g in panel_genera) df[[g]] <- as.numeric(getg(g))
df$degrader_sum <- rowSums(df[, degraders, drop = FALSE])
df <- df[!is.na(df$sarcosine) & !is.na(df$group), ]
df$group <- factor(df$group, levels = c("Healthy", "CRC"))
df$log10_sarcosine <- log10(df$sarcosine)
nH <- sum(df$group == "Healthy"); nC <- sum(df$group == "CRC")
cat(sprintf("[paired] %d samples (Healthy %d / CRC %d)\n", nrow(df), nH, nC))
stopifnot(nrow(df) == 308, nH == 121, nC == 187)

## ---------- 원본 대비 무결성 검증 ----------
o <- read.delim(file.path(ORIG, "merged_per_sample.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
o <- o[match(df$run, o$run), ]
stopifnot(identical(o$run, df$run))
d5_new <- rowSums(df[, BASE5, drop = FALSE]); d5_old <- rowSums(o[, BASE5, drop = FALSE])
cat(sprintf("[check] 공통 5속 합 최대 절대차 (원본 대비): %.3e\n", max(abs(d5_new - d5_old))))
cat(sprintf("[check] sarcosine 값 최대 절대차            : %.3e\n", max(abs(df$sarcosine - o$sarcosine))))
stopifnot(max(abs(d5_new - d5_old)) < 1e-12, max(abs(df$sarcosine - o$sarcosine)) < 1e-9)
write.table(df, file.path(OUT, "merged_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- (b) 분해속 ~ sarcosine 상관 ----------
vars <- c("degrader_sum", degraders)
cor.tab <- do.call(rbind, lapply(vars, function(v) {
  ct <- suppressWarnings(cor.test(df[[v]], df$sarcosine, method = "spearman"))
  data.frame(feature = v, rho = unname(ct$estimate), p = ct$p.value)
}))
cor.tab$q <- p.adjust(cor.tab$p, "BH")
cor.tab <- cor.tab[order(cor.tab$rho), ]
write.table(cor.tab, file.path(OUT, "DEGRADER_sarcosine_correlation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

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

## ---------- (c) 속 abundance: Healthy vs CRC ----------
allg <- c(BASE5, eligens, producers, blautia)
grp.tab <- do.call(rbind, lapply(allg, function(g) {
  p <- suppressWarnings(wilcox.test(df[[g]] ~ df$group)$p.value)
  data.frame(genus = g,
             role = ifelse(g %in% producers, "producer", ifelse(g == blautia, "mixed", "degrader")),
             median_Healthy = median(df[[g]][df$group == "Healthy"]),
             median_CRC = median(df[[g]][df$group == "CRC"]),
             log2FC_CRC_vs_H = log2((median(df[[g]][df$group == "CRC"]) + 1e-6) /
                                    (median(df[[g]][df$group == "Healthy"]) + 1e-6)), p = p)
}))
grp.tab$q <- p.adjust(grp.tab$p, "BH")
write.table(grp.tab, file.path(OUT, "DEGRADER_by_group.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- 원본 대비 핵심 수치 대조 ----------
ds.p  <- suppressWarnings(wilcox.test(df$degrader_sum ~ df$group)$p.value)
ds.mH <- median(df$degrader_sum[df$group == "Healthy"]); ds.mC <- median(df$degrader_sum[df$group == "CRC"])
ds.l2 <- log2(ds.mC / ds.mH)
rho.sum <- cor.tab$rho[cor.tab$feature == "degrader_sum"]; p.sum <- cor.tab$p[cor.tab$feature == "degrader_sum"]

o.cor <- read.delim(file.path(ORIG, "DEGRADER_sarcosine_correlation.tsv"))
o.rho <- o.cor$rho[o.cor$feature == "degrader_sum"]; o.p <- o.cor$p[o.cor$feature == "degrader_sum"]
o.sum <- rowSums(o[, c(BASE5, "Blautia"), drop = FALSE])
o.gp  <- suppressWarnings(wilcox.test(o.sum ~ df$group)$p.value)
o.mH  <- median(o.sum[df$group == "Healthy"]); o.mC <- median(o.sum[df$group == "CRC"])

cmp <- data.frame(
  metric = c("구성", "median Healthy", "median CRC", "log2FC (CRC/H)",
             "Wilcoxon p (H vs CRC)", "Spearman rho (vs sarcosine)", "Spearman p",
             "within-group p (Healthy)", "within-group p (CRC)"),
  ORIGINAL_5plusBlautia = c("5속 + Blautia", sprintf("%.4f", o.mH), sprintf("%.4f", o.mC),
                            sprintf("%.2f", log2(o.mC/o.mH)), sprintf("%.2e", o.gp),
                            sprintf("%.3f", o.rho), sprintf("%.2e", o.p),
                            sprintf("%.2f", wg$p[wg$feature=="degrader_sum" & wg$group=="Healthy"]),
                            sprintf("%.3f", wg$p[wg$feature=="degrader_sum" & wg$group=="CRC"])),
  VARIANT_5plusEligens  = c("5속 + eligens", sprintf("%.4f", ds.mH), sprintf("%.4f", ds.mC),
                            sprintf("%.2f", ds.l2), sprintf("%.2e", ds.p),
                            sprintf("%.3f", rho.sum), sprintf("%.2e", p.sum),
                            sprintf("%.2f", wgH), sprintf("%.3f", wgC)),
  stringsAsFactors = FALSE)
# 원본의 within-group 값은 원본 파일에서 읽어와 대체
o.wg <- read.delim(file.path(ORIG, "DEGRADER_sarcosine_correlation_withingroup.tsv"))
cmp$ORIGINAL_5plusBlautia[8] <- sprintf("%.2f", o.wg$p[o.wg$feature=="degrader_sum" & o.wg$group=="Healthy"])
cmp$ORIGINAL_5plusBlautia[9] <- sprintf("%.3f", o.wg$p[o.wg$feature=="degrader_sum" & o.wg$group=="CRC"])
write.table(cmp, file.path(OUT, "COMPARISON_vs_original.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=========== 원본 vs 변형본 (degrader_sum) ===========\n")
print(cmp, row.names = FALSE)

## ---------- 요약 ----------
sink(file.path(OUT, "summary.txt"))
cat("=== 04b VARIANT: degrader_sum = 5속 + [Eubacterium] eligens group (Blautia 제외) ===\n")
cat("원본 04_sarcosine_analysis.R 은 5속 + Blautia (eligens 제외) 로 계산한다.\n")
cat("이 변형본은 원고 Methods / Sup Fig 15a legend 의 서술을 정본으로 삼는다.\n\n")
cat(sprintf("페어드 샘플: %d (Healthy %d / CRC %d)\n\n", nrow(df), nH, nC))
cat("[원본 대비 대조]\n"); print(cmp, row.names = FALSE)
cat("\n[(b) 분해속 ~ sarcosine Spearman (pooled)]\n"); print(cor.tab, row.names = FALSE)
cat("\n[(b') within-group 상관]\n"); print(wg[wg$feature == "degrader_sum", ], row.names = FALSE)
cat("\n[(c) 속 abundance Healthy vs CRC]\n"); print(grp.tab, row.names = FALSE)
cat(sprintf("\n[결론] degrader_sum: median CRC %.4f vs Healthy %.4f, log2FC %.2f, p=%.3g\n",
            ds.mC, ds.mH, ds.l2, ds.p))
cat(sprintf("[한계 — 원본과 동일] pooled 음의 상관(rho=%.2f, p=%.3g)은 그룹 내에서 비유의(Healthy p=%.2f, CRC p=%.2f) -> 암 상태 교란. within-individual 관계 미확정.\n",
            rho.sum, p.sum, wgH, wgC))
sink()

cat("\n산출물 ->", OUT, "\n")
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))
