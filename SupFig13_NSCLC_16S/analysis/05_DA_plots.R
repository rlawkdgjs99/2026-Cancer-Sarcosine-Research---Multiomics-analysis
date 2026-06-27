## 05_DA_plots.R
## Visualize the top differentially abundant genera (NSCLC vs Health) as a
## diverging lollipop plot (primary) and a horizontal bar chart (alternative).
## Source: primary MaAsLin2 result (TSS+LOG, age-adjusted). Concordance marked vs CLR.

source("analysis/00_setup.R")
log_con <- file(file.path(LOG_DIR, "05_DA_plots.log"), open = "wt"); sink(log_con, split = TRUE)
cat("==== 05_DA_plots ====", as.character(Sys.time()), "\n")

TOPN <- 15                                   # top genera per direction to display

da     <- fread(file.path(RESULTS_DIR, "DA_maaslin2_TSS_LOG.csv"))
da_clr <- fread(file.path(RESULTS_DIR, "DA_maaslin2_CLR.csv"))
stopifnot(all(c("coef", "stderr", "qval", "genus", "lineage") %in% names(da)))

## concordance: also significant (q<0.05) under CLR normalization
da[, concordant := feature %in% da_clr[qval < 0.05, feature]]

## readable labels: strip GTDB genome-id suffix; for unclassified genera use family name
da[, label := gsub("_[0-9]+$", "", genus)]
ix <- grepl("^unclassified_genus", da$genus)
if (any(ix)) {
  fam <- sub("\\|g__.*$", "", da$lineage[ix]); fam <- sub(".*\\|f__", "", fam)
  fam[fam == "" | is.na(fam)] <- "unclassified"
  da$label[ix] <- paste0(gsub("_[0-9]+$", "", fam), " (g.uncl)")
}

## select significant genera, top N per direction by |effect size|
sig <- da[qval < ALPHA_SIG]
sig[, direction := ifelse(coef > 0, "Higher in NSCLC", "Higher in Health")]
sel <- sig[order(-abs(coef)), head(.SD, TOPN), by = direction]
sel[, label2 := ifelse(concordant, paste0(label, " *"), label)]
sel[, label2 := make.unique(label2)]                       # guarantee unique y labels
sel[, label2 := factor(label2, levels = sel[order(coef), label2])]  # order by effect
cat("significant (q<0.05):", nrow(sig), "| plotted:", nrow(sel),
    "| concordant w/ CLR among plotted:", sum(sel$concordant), "\n")

fwrite(sel[order(-coef), .(genus, label, direction, coef, stderr, pval, qval, concordant)],
       file.path(RESULTS_DIR, "DA_top_genera_plotted.csv"))

dir_cols <- c("Higher in Health" = unname(GROUP_COLORS["Health"]),
              "Higher in NSCLC"  = unname(GROUP_COLORS["NSCLC"]))
ttl  <- "Top DA genera: NSCLC vs Health"
subtitle <- sprintf("MaAsLin2 (TSS+LOG, age-adj.); q<0.05; top %d/dir. of %d\n*  = also significant under CLR",
                    TOPN, nrow(sig))
xlab <- "MaAsLin2 coefficient (effect size;  >0 = higher in NSCLC)"

## ---- lollipop (primary) ----
p_lolli <- ggplot(sel, aes(coef, label2, color = direction)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_segment(aes(x = coef - 1.96 * stderr, xend = coef + 1.96 * stderr,
                   y = label2, yend = label2), linewidth = 0.4, alpha = 0.5) +
  geom_segment(aes(x = 0, xend = coef, yend = label2), linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(size = -log10(qval))) +
  scale_color_manual(values = dir_cols, name = NULL) +
  scale_size_continuous(range = c(2, 5.5), name = expression(-log[10](q))) +
  labs(x = xlab, y = NULL, title = ttl, subtitle = subtitle) +
  theme_hvc + theme(legend.position = "right")
save_png(p_lolli, file.path(FIG_DIR, "DA_lollipop_top_genera.png"), 8.5, 8.6)
ggsave(file.path(FIG_DIR, "DA_lollipop_top_genera.pdf"), p_lolli, width = 8.5, height = 8.6)

## ---- horizontal bar (alternative) ----
p_bar <- ggplot(sel, aes(coef, label2, fill = direction)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_vline(xintercept = 0, color = "grey40") +
  scale_fill_manual(values = dir_cols, name = NULL) +
  labs(x = xlab, y = NULL, title = ttl, subtitle = subtitle) +
  theme_hvc + theme(legend.position = "right")
save_png(p_bar, file.path(FIG_DIR, "DA_bar_top_genera.png"), 8.5, 8.6)
ggsave(file.path(FIG_DIR, "DA_bar_top_genera.pdf"), p_bar, width = 8.5, height = 8.6)

cat("\n[plotted genera, ordered by effect]\n")
print(sel[order(-coef), .(label, direction, coef = round(coef, 2),
                          qval = signif(qval, 2), concordant)])
cat("\n==== done ====\n"); sink(); close(log_con)
