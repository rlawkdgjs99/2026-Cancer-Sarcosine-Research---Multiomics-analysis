#!/usr/bin/env Rscript
# ============================================================================
# TAXA <-> FUNCTION INTEGRATION (R vs NR): do Responder-enriched species associate
# with sarcosine DEGRADATION and Non-responder-enriched species with PRODUCTION?
#   For each DIFFERENTIALLY-ABUNDANT species (meta DA q<0.05), plot its meta Spearman
#   correlation with the DEGRADATION score (x) vs the PRODUCTION score (y), coloured by
#   the ICI-response group it is enriched in (R / NR). Tests the double dissociation with
#   a Wilcoxon of each correlation against enrichment direction.
# READ-ONLY on verified inputs (sarcosine_species_assoc_meta.csv = per-species meta rho/q;
#   pooled_species_meta.csv = DA). No upstream statistic recomputed. seed 42.
# Run: Rscript sarcosine_taxa_function_scatter.R
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
# Run with the R working directory set to the analysis-folder root; paths below are relative to it.
main_dir <- "."; res <- file.path(main_dir, "pooled_analysis", "results")
save_png <- function(p, out, w, h) { tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp) }

am <- fread(file.path(res, "sarcosine_species_assoc_meta.csv"))
w  <- dcast(am, species + da_dir + da_q ~ score, value.var = c("pooled_rho", "qval"))
setnames(w, c("pooled_rho_degradation","pooled_rho_production","qval_degradation","qval_production"),
            c("deg_rho","prod_rho","deg_q","prod_q"))
ds <- w[!is.na(da_q) & da_q < 0.05 & is.finite(deg_rho) & is.finite(prod_rho)]
ds[, grp := factor(fifelse(da_dir == "higher in R", "Enriched in Responder (R)", "Enriched in Non-responder (NR)"),
                   levels = c("Enriched in Responder (R)", "Enriched in Non-responder (NR)"))]

## double-dissociation tests (species-level Wilcoxon of rho by enrichment direction)
wd <- wilcox.test(deg_rho ~ da_dir, data = ds); wp <- wilcox.test(prod_rho ~ da_dir, data = ds)
fmtp <- function(p) ifelse(p < 1e-3, sprintf("%.1e", p), sprintf("%.3f", p))
cat(sprintf("DA-sig species: %d (R=%d, NR=%d)\n", nrow(ds), sum(ds$da_dir=="higher in R"), sum(ds$da_dir=="higher in NR")))
cat(sprintf("degradation-rho R-vs-NR Wilcoxon p=%s ; production-rho p=%s\n", fmtp(wd$p.value), fmtp(wp$p.value)))

## NOTE: H. hathewayi is intentionally NOT plotted here -- its differential abundance is NOT
## significant (q=0.55), so it is not part of the DA-significant species set. (Plotting a non-sig
## point on a figure of significant species would be misleading.) Its sarcosine-production
## association is shown in its own dedicated figure (Figure_Hungatella_hathewayi_main.png).

## labels: strongest R degraders + notable NR producers
labR <- head(ds[da_dir == "higher in R"][order(-deg_rho)], 6)
labNR <- ds[da_dir == "higher in NR"][order(deg_rho)][seq_len(min(6, .N))]
lab <- rbind(labR, labNR); lab[, sp := gsub("_", " ", species)]

COL <- c("Enriched in Responder (R)" = "#1B7837", "Enriched in Non-responder (NR)" = "#B2182B")
xlo <- min(ds$deg_rho) - 0.03; xhi <- max(ds$deg_rho) + 0.03
ylo <- min(ds$prod_rho) - 0.03; yhi <- max(ds$prod_rho) + 0.065

p <- ggplot(ds, aes(deg_rho, prod_rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(aes(color = grp), size = 3.1, alpha = 0.85) +
  ggrepel::geom_text_repel(data = lab, aes(label = sp, color = grp), size = 4.8, fontface = "italic",
                           max.overlaps = 20, min.segment.length = 0, segment.color = "grey75", show.legend = FALSE) +
  annotate("text", x = xhi - 0.005, y = ylo + 0.01, hjust = 1, vjust = 0, size = 5.7, fontface = "bold", color = "#1B7837",
           label = "Responder-enriched ->\nsarcosine DEGRADERS") +
  annotate("text", x = xlo + 0.005, y = yhi - 0.005, hjust = 0, vjust = 1, size = 5.7, fontface = "bold", color = "#B2182B",
           label = "Non-responder-enriched ->\nsarcosine PRODUCERS") +
  scale_color_manual(values = COL, name = "ICI response (DA direction)") +
  coord_cartesian(xlim = c(xlo, xhi), ylim = c(ylo, yhi)) +
  labs(x = expression("Correlation with sarcosine DEGRADATION score (meta Spearman " * rho * ")"),
       y = expression("Correlation with PRODUCTION score (meta " * rho * ")"),
       title = "Responder vs Non-responder bacteria split by sarcosine function",
       caption = sprintf("DA-significant species (meta q<0.05, n=%d) across 3 discovery cohorts.\nDouble dissociation Wilcoxon p: degradation-rho %s, production-rho %s (association, not causation).", nrow(ds), fmtp(wd$p.value), fmtp(wp$p.value))) +
  theme_classic(base_size = 16) +
  theme(legend.position = "top", legend.title = element_text(size = 16), legend.text = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 22),
        plot.caption = element_text(size = 12, hjust = 0, color = "grey40"))
save_png(p, file.path(res, "sarcosine_taxa_function_integration.png"), 10.5, 7.4)
fwrite(ds[order(da_dir, -deg_rho), .(species, da_dir, da_q, deg_rho, deg_q, prod_rho, prod_q)],
       file.path(res, "sarcosine_taxa_function_integration.csv"))
cat("wrote sarcosine_taxa_function_integration.png/.csv\n")
writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_taxa_function_scatter.txt"))
