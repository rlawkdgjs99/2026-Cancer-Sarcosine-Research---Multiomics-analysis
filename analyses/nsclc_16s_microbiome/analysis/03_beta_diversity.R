## 03_beta_diversity.R
## Between-sample (beta) diversity: Health vs NSCLC.
##  - Primary distance : Bray-Curtis on genus TSS relative abundance (HGMT paper).
##  - Sensitivity      : Aitchison (Euclidean on CLR) — compositionally aware.
##  - Test : PERMANOVA (adonis2, 999 perm), group-only AND covariate-adjusted.
##  - betadisper/permutest checks whether a PERMANOVA signal could be a dispersion artifact.

source("analysis/00_setup.R")
log_con <- file(file.path(LOG_DIR, "03_beta.log"), open = "wt"); sink(log_con, split = TRUE)
cat("==== 03_beta_diversity ====", as.character(Sys.time()), "\n")

prep <- readRDS(file.path(DERIVED_DIR, "prepared_16S.rds"))
meta <- prep$meta; tss <- prep$genus_tss; raw <- prep$genus_raw
stopifnot(all(meta$run_id == rownames(tss)))
md <- as.data.frame(meta)                              # adonis2 needs a data.frame

## ============ Primary: Bray-Curtis ============
bc <- vegdist(tss, method = "bray")

## PCoA (classical MDS)
pcoa <- cmdscale(bc, k = 2, eig = TRUE)
pos_eig <- pcoa$eig[pcoa$eig > 0]
var_expl <- round(100 * pcoa$eig[1:2] / sum(pos_eig), 1)
coords <- data.table(run_id = rownames(tss), group = meta$group,
                     PCo1 = pcoa$points[, 1], PCo2 = pcoa$points[, 2])
fwrite(coords, file.path(RESULTS_DIR, "beta_pcoa_braycurtis_coords.csv"))

## PERMANOVA — group only
set.seed(42)
perm_grp <- adonis2(bc ~ group, data = md, permutations = N_PERM, by = "terms")
## PERMANOVA — covariate-adjusted (marginal effect of group given age/sex/BMI)
set.seed(42)
perm_adj <- adonis2(bc ~ age + sex + bmi + group, data = md, permutations = N_PERM, by = "margin")
cat("\n[PERMANOVA Bray-Curtis ~ group]\n"); print(perm_grp)
cat("\n[PERMANOVA Bray-Curtis ~ age+sex+bmi+group (marginal)]\n"); print(perm_adj)

## dispersion check (betadisper)
set.seed(42)
bd <- betadisper(bc, meta$group)
bd_test <- permutest(bd, permutations = N_PERM)
cat("\n[betadisper permutest — dispersion homogeneity]\n"); print(bd_test)

## tidy PERMANOVA tables to CSV
save_adonis <- function(a, file) {
  d <- as.data.frame(a); d$term <- rownames(d); fwrite(as.data.table(d), file)
}
save_adonis(perm_grp, file.path(RESULTS_DIR, "beta_permanova_bray_grouponly.csv"))
save_adonis(perm_adj, file.path(RESULTS_DIR, "beta_permanova_bray_adjusted.csv"))
bd_df <- as.data.frame(bd_test$tab); bd_df$term <- rownames(bd_df)
fwrite(as.data.table(bd_df), file.path(RESULTS_DIR, "beta_betadisper_bray.csv"))

## ============ Sensitivity: Aitchison (CLR + Euclidean) ============
## CLR is scale-invariant per sample; add a documented pseudocount for zeros.
pseudo <- min(tss[tss > 0]) / 2                        # half smallest non-zero proportion
clr <- t(apply(tss + pseudo, 1, function(z) { lz <- log(z); lz - mean(lz) }))
ait <- dist(clr, method = "euclidean")
set.seed(42)
perm_ait <- adonis2(ait ~ group, data = md, permutations = N_PERM, by = "terms")
cat("\n[PERMANOVA Aitchison ~ group] (pseudocount =", signif(pseudo, 3), ")\n"); print(perm_ait)
save_adonis(perm_ait, file.path(RESULTS_DIR, "beta_permanova_aitchison_grouponly.csv"))

## ============ PCoA plot (Bray-Curtis) ============
r2  <- perm_grp$R2[1]; pval <- perm_grp$`Pr(>F)`[1]
p <- ggplot(coords, aes(PCo1, PCo2, color = group)) +
  geom_point(size = 1.8, alpha = 0.75) +
  stat_ellipse(level = 0.68, linewidth = 0.6) +
  scale_color_manual(values = GROUP_COLORS) +
  labs(x = sprintf("PCo1 (%.1f%%)", var_expl[1]),
       y = sprintf("PCo2 (%.1f%%)", var_expl[2]),
       title = "PCoA (Bray-Curtis)",
       subtitle = sprintf("PERMANOVA R2 = %.3f, P = %.3g", r2, pval)) +
  theme_hvc
save_png(p, file.path(FIG_DIR, "beta_pcoa_braycurtis.png"), 5, 4)
ggsave(file.path(FIG_DIR, "beta_pcoa_braycurtis.pdf"), p, width = 5, height = 4)

cat("\n==== done ====\n"); sink(); close(log_con)
