#!/usr/bin/env Rscript
# ============================================================================
# MAIN FIGURE — Hungatella hathewayi in NSCLC (CRC<->NSCLC sarcosine thread).
#   (A) Sarcosine-PRODUCTION correlation forest  [ROBUST: 4/4 cohorts +, meta q<1e-4, I2=0]
#   (B) Relative abundance, ICI Responder vs Non-responder  [pooled discovery, 3 cohorts]
# Honest framing: the BACKBONE is the production-association (A) + the CRC literature.
#   Reads verified inputs; recomputes only display panels. cowplot. seed 42.
#   (The earlier overall-survival KM panel [C] was removed per user request, 2026-06-09.)
# Run: Rscript figure_hungatella_main.R
# ============================================================================
set.seed(42); Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
# Run with the R working directory set to the analysis-folder root; paths below are relative to it.
main_dir <- "."; res <- file.path(main_dir, "pooled_analysis", "results")
save_png <- function(p, out, w, h) { tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp) }
read_species <- function(cdir) { f <- list.files(cdir, "^Bacteria_.*\\.txt$", full.names = TRUE); stopifnot(length(f) == 1L)
  b <- fread(f, sep = "\t", header = TRUE, quote = "", showProgress = FALSE); stopifnot(ncol(b) == 3L); setnames(b, c("Taxa","RunID","Abundance"))
  sp <- b[grepl("|s__", Taxa, fixed = TRUE) & !grepl("|t__", Taxa, fixed = TRUE)]
  w <- dcast(sp, RunID ~ Taxa, value.var = "Abundance", fill = 0, fun.aggregate = sum)
  mm <- as.matrix(w[, -1L, with = FALSE]); rownames(mm) <- w$RunID; colnames(mm) <- make.unique(sub(".*\\|s__", "", colnames(mm))); mm }
HH <- "Hungatella_hathewayi"

## ---- Panel A: production-correlation forest (read verified CSV) -------------
fo <- fread(file.path(res, "suppl_Hungatella_production_forest.csv"))
ord <- fo$label; fo[, label := factor(label, levels = rev(ord))]
fo[, role := factor(role, levels = c("discovery","validation","pooled"))]
xlabpos <- max(fo$hi) + 0.06
pA <- ggplot(fo, aes(rho, label, color = role)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.6) +
  geom_point(size = 2.9) +
  geom_text(aes(x = xlabpos, label = stat_label), hjust = 0, size = 3.4, color = "grey25") +
  scale_color_manual(values = c(discovery = "#0072B2", validation = "#999999", pooled = "#117733"), guide = "none") +
  scale_x_continuous(limits = c(min(fo$lo) - 0.03, xlabpos + 0.78)) +
  labs(x = expression("Spearman " * rho * "  (abundance vs sarcosine production score)"), y = NULL,
       title = expression(bold("Sarcosine-production correlation")~" (q<1e-4, "*I^2*"=0)")) +
  theme_classic(base_size = 16) + theme(plot.title = element_text(size = 14), axis.text.y = element_text(size = 12))

## ---- Panel B: abundance R vs NR (pooled discovery, 3 cohorts) ---------------
DISC <- c(PRJNA751792 = "NSCLC_PRJNA751792", PRJNA1023797 = "NSCLC_PRJNA1023797", PRJEB22863 = "NSCLC_RCC_PRJEB22863")
bx <- rbindlist(lapply(names(DISC), function(cn) {
  cd <- file.path(main_dir, DISC[cn]); sc <- fread(file.path(cd, "results", paste0("sarcosine_scores_", cn, ".csv")))
  spm <- read_species(cd); ab <- if (HH %in% colnames(spm)) spm[sc$RunID, HH] else rep(0, nrow(sc))
  data.table(group = as.character(sc$group), hh = as.numeric(ab)) }))
bx[, group := factor(group, levels = c("R","NR"))]
pseudo <- min(bx$hh[bx$hh > 0]) / 2; bx[, logab := log10(hh + pseudo)]
da <- fread(file.path(res, "pooled_species_meta.csv"))[species == HH]
annB <- sprintf("meta DA: %s, q=%.2f", sub("higher in ", "", da$dir), da$qval)
pB <- ggplot(bx, aes(group, logab, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55, linewidth = 0.4) +
  geom_jitter(width = 0.15, height = 0, size = 0.8, alpha = 0.35, color = "grey30") +
  annotate("text", x = 1.5, y = max(bx$logab) + 0.10 * diff(range(bx$logab)), label = annB, size = 3.9, fontface = "bold") +
  scale_fill_manual(values = c(R = "#1B7837", NR = "#B2182B"), guide = "none") +
  scale_x_discrete(limits = c("R", "NR")) +
  coord_cartesian(ylim = c(min(bx$logab), max(bx$logab) + 0.18 * diff(range(bx$logab)))) +
  labs(x = "ICI response", y = expression(log[10]~"(relative abundance %)"),
       title = expression(bold("Abundance: Responder vs Non-responder"))) +
  theme_bw(base_size = 16) + theme(plot.title = element_text(size = 14), panel.grid.minor = element_blank())

## ---- assemble (panels A + B; the overall-survival KM panel C was removed) ---
fig <- plot_grid(pA, pB, ncol = 1, labels = c("A", "B"), rel_heights = c(0.85, 1), label_size = 18)
ttl <- ggdraw() + draw_label("Hungatella hathewayi in NSCLC gut microbiome (ICI cohorts)", fontface = "bold.italic", size = 18, x = 0.01, hjust = 0)
fig2 <- plot_grid(ttl, fig, ncol = 1, rel_heights = c(0.06, 1))
save_png(fig2, file.path(res, "Figure_Hungatella_hathewayi_main.png"), 8.5, 8)
cat("wrote Figure_Hungatella_hathewayi_main.png | panel B meta DA:", annB, "\n")
writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_figure_hungatella_main.txt"))
