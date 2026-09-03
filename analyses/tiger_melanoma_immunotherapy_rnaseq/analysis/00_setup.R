# =============================================================================
# 00_setup.R  -- shared colors, theme, and helpers (sourced by 02/03)
# User preferences (standing):
#   - Group colors: ICI Responder (R) = GREEN ; Non-responder (NR) = RED
#   - KM survival : expression Low = GREEN ; High = RED
#   - Fonts large for readability; TITLE noticeably larger; titles short
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(survival)
  library(survminer)
})

# ---- Color constants --------------------------------------------------------
COL_GREEN <- "#1F9E5A"   # Responder / Low expression
COL_RED   <- "#D6291E"   # Non-responder / High expression

RESP_COLORS <- c(R = COL_GREEN, NR = COL_RED)        # R vs NR boxplots
EXPR_COLORS <- c(Low = COL_GREEN, High = COL_RED)    # KM Low vs High

# ---- Publication theme (large fonts, larger short title) --------------------
theme_pub <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title   = element_text(size = 22, face = "bold", hjust = 0.5,
                                  margin = margin(b = 8)),
      plot.subtitle= element_text(size = 15, hjust = 0.5, margin = margin(b = 6)),
      axis.title   = element_text(size = 17),
      axis.text    = element_text(size = 15, color = "black"),
      axis.line    = element_line(linewidth = 0.6),
      axis.ticks   = element_line(linewidth = 0.6, color = "black"),
      legend.title = element_text(size = 15),
      legend.text  = element_text(size = 14),
      legend.position = "none",
      plot.margin  = margin(10, 14, 10, 10)
    )
}

# ---- Save helper: write PNG (300 dpi) only — for PowerPoint figure planning --
# NOTE: the project path contains non-ASCII (Korean) characters. The default
# ragg/agg PNG device FAILS on such paths ("agg could not write"); cairo_pdf and
# cairo-type png() both work. We therefore force device = png, type = "cairo".
#
# Figures are auto-routed into results/figures/<analysis>/<treatment>/ by stem:
#   treatment: stem "sensitivity_*" -> outlier_removed ; "tiger_*" -> TIGER_all_samples
#              ; otherwise -> primary_PRE_only
#   analysis : stem containing "box" -> 01_response_RvsNR ; "km" -> 02_survival_KM
.fig_dir <- function(stem) {
  treat <- if (grepl("^sensitivity_", stem)) "outlier_removed"
           else if (grepl("^tiger_", stem)) "TIGER_all_samples"
           else "primary_PRE_only"
  body  <- sub("^(sensitivity|tiger)_", "", stem)
  atype <- if (grepl("box", body)) "01_response_RvsNR"
           else if (grepl("km", body)) "02_survival_KM"
           else "99_misc"
  d <- file.path("results/figures", atype, treat)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}
save_plot <- function(plot, file_stem, width = 5, height = 5.4) {
  d <- .fig_dir(file_stem)
  # PNG only (for PowerPoint figure planning). To also emit vector PDF for final
  # publication, add: ggsave(file.path(d, paste0(file_stem,".pdf")), plot,
  #   width=width, height=height, device=cairo_pdf)
  ggsave(file.path(d, paste0(file_stem, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white",
         device = png, type = "cairo")
  invisible(plot)
}

# ---- Pretty p-value formatter ----------------------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 1e-4) return("P < 0.0001")
  paste0("P = ", formatC(p, format = "f", digits = ifelse(p < 0.001, 4, 3)))
}
