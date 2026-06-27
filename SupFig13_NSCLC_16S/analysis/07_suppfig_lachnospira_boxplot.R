## 07_suppfig_lachnospira_boxplot.R
## Supplementary figure: Lachnospira genus relative abundance (16S, TSS%), Health vs NSCLC.
## Self-contained on purpose: it does NOT source 00_setup.R, whose BASE_DIR points to a
## different machine/user. Edit BASE below if the project is relocated.
## Stats are NOT recomputed for the annotation — the MaAsLin2 q-values are read straight
## from the existing pipeline results (results/DA_maaslin2_*.csv) for integrity.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE     <- "."
ANALYSIS <- file.path(BASE, "analysis")
DERIVED  <- file.path(ANALYSIS, "data_derived")
RESULTS  <- file.path(ANALYSIS, "results")
FIG      <- file.path(ANALYSIS, "figures")

## ---- house style (copied from 00_setup.R) ----
GROUP_COLORS <- c(Health = "#1B7837", NSCLC = "#B2182B")
theme_hvc <- theme_bw(base_size = 16) + theme(panel.grid.minor = element_blank())
## macOS-native quartz backend (cairo/X11 not available on this machine; ragg breaks on the
## non-ASCII project path). quartz handles both the Unicode path and Unicode glyphs cleanly.
save_png <- function(plot, file, width, height, dpi = 300) {
  grDevices::png(file, width = width, height = height, units = "in", res = dpi, type = "quartz")
  print(plot); grDevices::dev.off(); invisible(file)
}

## ---- data: Lachnospira TSS relative abundance (%) ----
prep <- readRDS(file.path(DERIVED, "prepared_16S.rds"))
tss  <- prep$genus_tss; meta <- prep$meta
stopifnot("Lachnospira" %in% colnames(tss), all(meta$run_id == rownames(tss)))

dt <- data.table(run_id = rownames(tss),
                 group  = meta$group,
                 Lachnospira_pct = tss[, "Lachnospira"] * 100)   # TSS proportion -> %

## ---- authoritative stats from the pipeline (read, do not recompute) ----
da   <- fread(file.path(RESULTS, "DA_maaslin2_TSS_LOG.csv"))
dac  <- fread(file.path(RESULTS, "DA_maaslin2_CLR.csv"))
qval <- da [genus == "Lachnospira", qval]
pval <- da [genus == "Lachnospira", pval]   # raw p from the SAME age-adjusted model (q is its BH-adjusted form)
qclr <- dac[genus == "Lachnospira", qval]
stopifnot(length(qval) == 1, length(pval) == 1, length(qclr) == 1)

## ---- group summary + source data ----
summ <- dt[, .(n = .N, detected = sum(Lachnospira_pct > 0),
               median_pct = round(median(Lachnospira_pct), 4),
               mean_pct   = round(mean(Lachnospira_pct), 4)), by = group]
cat("[group summary]\n"); print(summ)
fwrite(dt,   file.path(RESULTS, "suppfig_Lachnospira_TSS_per_sample.csv"))
fwrite(summ, file.path(RESULTS, "suppfig_Lachnospira_summary.csv"))

## ---- plot builder (linear or log10) ----
build <- function(scale = c("linear", "log10")) {
  scale <- match.arg(scale)
  d <- copy(dt)
  if (scale == "log10") {
    pseudo <- min(d$Lachnospira_pct[d$Lachnospira_pct > 0]) / 2   # documented: half min non-zero
    d[, yval := Lachnospira_pct + pseudo]
    ylab <- sprintf("Lachnospira rel. abundance (%%)\n[log10 scale, +%.1g pseudocount]", pseudo)
    sc   <- scale_y_log10(); ypos <- max(d$yval) * 1.6
  } else {
    d[, yval := Lachnospira_pct]
    ylab <- "Lachnospira relative abundance (%)"
    sc   <- NULL; ypos <- max(d$yval) * 1.03
  }
  ggplot(d, aes(group, yval, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.45) +
    sc +
    scale_fill_manual(values = GROUP_COLORS) +
    annotate("text", x = 1.5, y = ypos,
             label = sprintf("q = %.2g  (p = %.2g)", qval, pval), size = 5.3) +
    labs(x = NULL, y = ylab,
         title = "Lachnospira: Health vs NSCLC") +
    theme_hvc + theme(legend.position = "none")
}

p_lin <- build("linear")
p_log <- build("log10")
save_png(p_lin, file.path(FIG, "suppfig_Lachnospira_boxplot_linear.png"), 5.2, 4.6)
ggsave(file.path(FIG, "suppfig_Lachnospira_boxplot_linear.pdf"), p_lin, width = 5.2, height = 4.6)
save_png(p_log, file.path(FIG, "suppfig_Lachnospira_boxplot_log10.png"), 5.2, 4.6)
ggsave(file.path(FIG, "suppfig_Lachnospira_boxplot_log10.pdf"), p_log, width = 5.2, height = 4.6)

cat("\n[done] wrote linear + log10 PNG/PDF and 2 source CSVs to figures/ and results/\n")
cat(sprintf("MaAsLin2 age-adj q = %.4g | CLR q = %.4g\n", qval, qclr))
