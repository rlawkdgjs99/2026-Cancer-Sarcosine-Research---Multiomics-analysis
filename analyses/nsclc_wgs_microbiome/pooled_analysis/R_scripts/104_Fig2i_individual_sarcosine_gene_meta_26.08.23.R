#!/usr/bin/env Rscript
# Re-render the existing Sup. Fig. 11a individual-KO meta-analysis as a
# compact Main-Figure candidate. This script is visualization-only: it reads
# the locked upstream DerSimonian-Laird meta-analysis estimates and BH q values
# and does not refit a model or alter any statistic.
#
# Statistical design inherited from pooled_sarcosine.R:
# - unit displayed: one sarcosine-associated KEGG orthologue (KO)
# - discovery design: three independent NSCLC ICI WGS cohorts
# - upstream cohort effects: MaAsLin2 R-versus-NR coefficients
# - upstream synthesis: inverse-variance random-effects meta-analysis (DL)
# - upstream inclusion: estimable in >=2 discovery cohorts
# - upstream multiplicity: Benjamini-Hochberg across the 7 included KOs
#
# Display rule requested on 2026-08-23:
# - y axis contains gene/enzyme name plus functional role only
# - KO identifiers are retained in source data/provenance, not on the axis
# - official KEGG symbols are used where present; for K08687 and K08688 KEGG
#   provides enzyme-number-style SYMBOL fields, so their verified enzyme names
#   are used instead of inventing gene symbols.

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(ggplot2)
  library(ggtext)
  library(here)
  library(patchwork)
  library(ragg)
  library(svglite)
})

stopifnot(file.exists(here::here("PROJECT_HANDOFF.md")))

nsc_root <- here::here(
  "공공_Metabolomics&Metagenomics_분석모음",
  "HGMT_NSCLC_ICI_RvsNR_WGS"
)
pooled_dir <- file.path(nsc_root, "pooled_analysis")
result_dir <- file.path(pooled_dir, "results")
out_dir <- file.path(result_dir, "Fig2i_individual_sarcosine_gene_meta_26.08.23")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_meta <- file.path(result_dir, "pooled_sarcosine_KO_meta.csv")
input_koset <- file.path(nsc_root, "sarcosine_KO_set.csv")
input_upstream_script <- file.path(pooled_dir, "R_scripts", "pooled_sarcosine.R")
input_selected_source <- here::here("Figure_Panel_Source_Data_26.07.23", "SupFig9b.csv")
input_files <- c(input_meta, input_koset, input_upstream_script, input_selected_source)
stopifnot(all(file.exists(input_files)))

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

meta <- fread(input_meta)
koset <- fread(input_koset)
selected <- fread(input_selected_source)

# Input inspection and integrity checks before plotting.
required_meta <- c(
  "KO", "k_cohorts", "pooled_coef", "ci_lb", "ci_ub", "pval",
  "I2", "tau2", "dir", "qval", "valid_coef", "valid_q"
)
required_selected <- c("KO", "pooled_coef", "ci_lb", "ci_ub", "role_label")
stopifnot(
  identical(names(meta), required_meta),
  all(required_selected %in% names(selected)),
  nrow(meta) == 7L,
  !anyDuplicated(meta$KO),
  !anyNA(meta),
  all(vapply(meta[, .(pooled_coef, ci_lb, ci_ub, pval, I2, tau2, qval)],
             function(x) all(is.finite(x)), logical(1))),
  all(meta$ci_lb <= meta$pooled_coef),
  all(meta$pooled_coef <= meta$ci_ub),
  all(meta$qval >= 0 & meta$qval <= 1),
  all(meta$k_cohorts >= 2L)
)

selected_check <- merge(
  meta[, .(KO, pooled_coef, ci_lb, ci_ub)],
  selected[, .(KO, pooled_coef_selected = pooled_coef,
               ci_lb_selected = ci_lb, ci_ub_selected = ci_ub)],
  by = "KO", all = TRUE
)
stopifnot(
  nrow(selected_check) == 7L,
  !anyNA(selected_check),
  max(abs(selected_check$pooled_coef - selected_check$pooled_coef_selected)) < 1e-12,
  max(abs(selected_check$ci_lb - selected_check$ci_lb_selected)) < 1e-12,
  max(abs(selected_check$ci_ub - selected_check$ci_ub_selected)) < 1e-12
)

# Frozen annotation snapshot checked against official KEGG KO entries on
# 2026-08-23. KEGG URLs are retained for traceability; the plot does not make a
# live web request, so reruns are deterministic.
annotation <- data.table(
  KO = c("K00302", "K00303", "K08687", "K08688", "K18897", "K21833", "K21834"),
  kegg_symbol = c("soxA", "soxB", "E3.5.1.59", "E3.5.3.3", "sdmt", "dgcA, ddhC", "dgcB"),
  display_gene = c(
    "soxA", "soxB", "N-carbamoylsarcosine amidase", "creatinase",
    "sdmt", "dgcA/ddhC", "dgcB"
  ),
  display_label_type = c(
    "KEGG gene symbol", "KEGG gene symbol", "verified KEGG enzyme name",
    "verified KEGG enzyme name", "KEGG gene symbol", "KEGG gene symbols",
    "KEGG gene symbol"
  ),
  role = c(
    "degradation", "degradation", "production", "production",
    "degradation", "degradation;production", "degradation;production"
  ),
  role_display = c(
    "Degradation", "Degradation", "Production", "Production",
    "Degradation", "Bidirectional", "Bidirectional"
  ),
  kegg_entry_url = paste0("https://www.kegg.jp/entry/", c(
    "K00302", "K00303", "K08687", "K08688", "K18897", "K21833", "K21834"
  )),
  annotation_checked_date = "2026-08-23"
)

annotation[, gene_axis_html := c(
  "<i>soxA</i>",
  "<i>soxB</i>",
  "N-carbamoylsarcosine amidase",
  "creatinase",
  "<i>sdmt</i>",
  "<i>dgcA</i>/<i>ddhC</i>",
  "<i>dgcB</i>"
)]

# Confirm functional roles against the locally verified pathway definition.
annotation_check <- merge(
  annotation,
  koset[, .(KO, KO_name, role_local = role)],
  by = "KO", all.x = TRUE
)
stopifnot(
  nrow(annotation_check) == 7L,
  !anyNA(annotation_check$KO_name),
  identical(annotation_check$role, annotation_check$role_local)
)

plot_data <- merge(meta, annotation, by = "KO", all.x = TRUE)
stopifnot(nrow(plot_data) == 7L, !anyNA(plot_data$display_gene))

# Biological ordering is fixed a priori by functional role, not by effect size.
display_order <- c("K00302", "K00303", "K18897", "K08687", "K08688", "K21833", "K21834")
plot_data[, display_rank := match(KO, display_order)]
stopifnot(!anyNA(plot_data$display_rank))
setorder(plot_data, display_rank)
plot_data[, y := rev(seq_len(.N))]
plot_data[, axis_label := paste0(
  gene_axis_html,
  " <span style='color:#666666'>(", role_display, ")</span>"
)]
plot_data[, direction := factor(
  ifelse(pooled_coef >= 0, "R higher", "NR higher"),
  levels = c("R higher", "NR higher")
)]
plot_data[, significant := qval < 0.05]
plot_data[, q_label := ifelse(qval < 0.001, "<0.001", sprintf("%.3f", qval))]
plot_data[, I2_label := ifelse(I2 < 0.05, "0", sprintf("%.1f", I2))]

RESP_COLS <- c("R higher" = "#2E5F8A", "NR higher" = "#C47B3B")
ROW_SEPARATORS <- c(4.5, 2.5)

forest <- ggplot(plot_data, aes(y = y)) +
  geom_hline(yintercept = ROW_SEPARATORS, colour = "#D9D9D9", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "#8C8C8C", linewidth = 0.45,
             linetype = "22") +
  geom_segment(
    aes(x = ci_lb, xend = ci_ub, yend = y, colour = direction),
    linewidth = 0.9, lineend = "round", show.legend = FALSE
  ) +
  geom_segment(
    aes(x = ci_lb, xend = ci_lb, y = y - 0.085, yend = y + 0.085,
        colour = direction),
    linewidth = 0.7, show.legend = FALSE
  ) +
  geom_segment(
    aes(x = ci_ub, xend = ci_ub, y = y - 0.085, yend = y + 0.085,
        colour = direction),
    linewidth = 0.7, show.legend = FALSE
  ) +
  geom_point(
    aes(x = pooled_coef, fill = direction, colour = direction, alpha = significant),
    shape = 21, size = 3.0, stroke = 0.45, show.legend = FALSE
  ) +
  scale_colour_manual(values = RESP_COLS) +
  scale_fill_manual(values = RESP_COLS) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.70)) +
  scale_x_continuous(
    limits = c(-0.78, 1.28),
    breaks = c(-0.5, 0, 0.5, 1.0),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0.5, 8.0),
    breaks = plot_data$y,
    labels = plot_data$axis_label,
    expand = expansion(mult = 0)
  ) +
  labs(
    x = "Pooled MaAsLin2 coefficient (R − NR)\n← higher in NR                         higher in R →",
    y = NULL
  ) +
  theme_classic(base_family = "Arial", base_size = 7.0) +
  theme(
    axis.title.x = element_text(size = 7.0, margin = margin(t = 5)),
    axis.text.x = element_text(size = 6.5, colour = "black"),
    axis.text.y = ggtext::element_markdown(
      size = 6.7, colour = "black", hjust = 1, margin = margin(r = 5)
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(2, 2, 3, 4)
  )

stats_panel <- ggplot(plot_data, aes(y = y)) +
  geom_hline(yintercept = ROW_SEPARATORS, colour = "#D9D9D9", linewidth = 0.35) +
  geom_text(
    aes(x = 1, label = q_label, fontface = ifelse(significant, "bold", "plain")),
    family = "Arial", size = 6.5 / ggplot2::.pt, colour = "black"
  ) +
  geom_text(
    aes(x = 2, label = I2_label),
    family = "Arial", size = 6.5 / ggplot2::.pt, colour = "#4D4D4D"
  ) +
  annotate(
    "text", x = 1, y = 7.78, label = "q", family = "Arial",
    fontface = "bold", size = 6.5 / ggplot2::.pt
  ) +
  annotate(
    "text", x = 2, y = 7.78, label = "I² (%)", family = "Arial",
    fontface = "bold", size = 6.5 / ggplot2::.pt
  ) +
  scale_x_continuous(limits = c(0.45, 2.45), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0.5, 8.0), expand = expansion(mult = 0)) +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(2, 4, 3, 0))

combined <- (forest | stats_panel) +
  plot_layout(widths = c(4.9, 1.15)) +
  plot_annotation(
    title = "Individual microbial sarcosine genes",
    subtitle = "Random-effects meta-analysis across 3 NSCLC ICI discovery cohorts",
    theme = theme(
      plot.title = element_text(
        family = "Arial", face = "bold", size = 8.5,
        hjust = 0, margin = margin(b = 1)
      ),
      plot.subtitle = element_text(
        family = "Arial", size = 6.5, colour = "#4D4D4D",
        hjust = 0, margin = margin(b = 3)
      ),
      plot.margin = margin(4, 5, 4, 4)
    )
  )

png_ppt <- file.path(out_dir, "Proposed_Fig2i_individual_sarcosine_gene_meta_PPT.png")
png_600 <- file.path(out_dir, "Proposed_Fig2i_individual_sarcosine_gene_meta_600dpi.png")
svg_path <- file.path(out_dir, "Proposed_Fig2i_individual_sarcosine_gene_meta.svg")

save_via_ascii_temp <- function(path, device, dpi = NULL) {
  extension <- paste0(".", tools::file_ext(path))
  tmp <- tempfile(pattern = "fig2i_gene_meta_", fileext = extension)
  args <- list(
    filename = tmp, plot = combined, device = device,
    width = 5.9, height = 3.25, units = "in", bg = "white"
  )
  if (!is.null(dpi)) args$dpi <- dpi
  do.call(ggsave, args)
  stopifnot(file.exists(tmp), file.info(tmp)$size > 0)
  stopifnot(file.copy(tmp, path, overwrite = TRUE))
  unlink(tmp)
}

save_via_ascii_temp(png_ppt, ragg::agg_png, dpi = 300)
save_via_ascii_temp(png_600, ragg::agg_png, dpi = 600)
save_via_ascii_temp(svg_path, svglite::svglite)

source_data_path <- file.path(out_dir, "Fig2i_individual_sarcosine_gene_meta_source_data.csv")
annotation_path <- file.path(out_dir, "KEGG_gene_annotation_snapshot_2026-08-23.csv")
input_manifest_path <- file.path(out_dir, "input_manifest_sha256.tsv")
session_path <- file.path(out_dir, "sessionInfo.txt")
log_path <- file.path(out_dir, "analysis_log.txt")
readme_path <- file.path(out_dir, "README.txt")

fwrite(
  plot_data[, .(
    display_order = display_rank, KO, display_gene, role_display, k_cohorts,
    pooled_coef, ci_lb, ci_ub, pval, qval, I2, tau2, direction,
    validation_coef = valid_coef, validation_q = valid_q
  )],
  source_data_path
)
fwrite(annotation_check, annotation_path)
fwrite(
  data.table(
    file = normalizePath(input_files, mustWork = TRUE),
    sha256 = vapply(input_files, sha256_file, character(1))
  ),
  input_manifest_path,
  sep = "\t"
)

writeLines(capture.output(sessionInfo()), session_path)

writeLines(c(
  "Figure purpose",
  "--------------",
  "Compact redraw of the former Supplementary Fig. 11a individual sarcosine-KO meta-analysis for proposed placement as Main Fig. 2i.",
  "",
  "Statistical integrity",
  "---------------------",
  "No model was refitted. Point estimates, 95% CIs, BH q values and I2 values are copied from pooled_analysis/results/pooled_sarcosine_KO_meta.csv.",
  "The upstream analysis combined MaAsLin2 coefficients from the 3 discovery NSCLC ICI cohorts using a DerSimonian-Laird random-effects model.",
  "The upstream inclusion criterion was estimability in at least 2 discovery cohorts. All 7 qualifying KOs are displayed; none was removed for appearance or significance.",
  "Positive coefficients mean higher abundance in responders (R); negative coefficients mean higher abundance in non-responders (NR).",
  "",
  "Axis-label policy",
  "-----------------",
  "The y axis displays only gene/enzyme name and functional role, as requested; KO identifiers remain in the source-data CSV.",
  "Official KEGG gene symbols are used for soxA, soxB, sdmt, dgcA/ddhC and dgcB.",
  "K08687 and K08688 have enzyme-number-style KEGG SYMBOL entries rather than conventional gene symbols, so the verified enzyme names N-carbamoylsarcosine amidase and creatinase are shown instead of inventing symbols.",
  "KEGG annotation was checked on 2026-08-23; exact entry URLs are stored in KEGG_gene_annotation_snapshot_2026-08-23.csv.",
  "",
  "Display files",
  "-------------",
  "Use Proposed_Fig2i_individual_sarcosine_gene_meta_PPT.png for PowerPoint placement.",
  "Use the 600-dpi PNG or SVG for publication assembly/export.",
  "The proposed Fig. 2i identity remains provisional until the author inserts the panel into the deck."
), readme_path)

writeLines(c(
  paste0("Generated: ", format(Sys.time(), tz = "Asia/Seoul", usetz = TRUE)),
  paste0("Rows: ", nrow(plot_data)),
  paste0("KOs with q < 0.05: ", paste(plot_data[qval < 0.05, KO], collapse = ", ")),
  paste0("Input meta SHA256: ", sha256_file(input_meta)),
  paste0("Selected source-data match max |delta|: ",
         format(max(abs(selected_check$pooled_coef - selected_check$pooled_coef_selected)), scientific = TRUE)),
  "No statistical re-estimation performed."
), log_path)

output_files <- c(
  png_ppt, png_600, svg_path, source_data_path, annotation_path,
  input_manifest_path, session_path, log_path, readme_path
)
stopifnot(all(file.exists(output_files)), all(file.info(output_files)$size > 0))

output_manifest_path <- file.path(out_dir, "output_manifest_sha256.tsv")
fwrite(
  data.table(
    file = basename(output_files),
    sha256 = vapply(output_files, sha256_file, character(1)),
    bytes = file.info(output_files)$size
  ),
  output_manifest_path,
  sep = "\t"
)

message("Created: ", normalizePath(png_ppt))
message("All statistics matched the locked source table; no model was refitted.")
