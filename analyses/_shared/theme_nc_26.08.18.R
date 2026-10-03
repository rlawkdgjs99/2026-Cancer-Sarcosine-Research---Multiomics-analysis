# Shared Nature-Communications figure theme (2026-08-18)
#
# WHY THIS EXISTS
# The 26.07.23 pipelines export each panel onto a 6.2-14 inch canvas and the panel is then
# placed into the figure at 1.5-4.3 inches. The resulting 60-70% shrink turns 15 pt tick
# text into 4.5-5.5 pt, below the legibility floor for print. This theme instead renders
# every panel AT ITS FINAL PLACED SIZE, so a point size written here is the point size the
# reader sees. Panels rendered with it must be placed at their declared footprint and NOT
# rescaled in PowerPoint.
#
# Type scale (final, on-page):
NC_TICK_PT   <- 6.0     # axis tick labels
NC_TITLE_PT  <- 7.0     # axis titles
NC_STRIP_PT  <- 6.5     # facet strip labels
NC_ANNOT_PT  <- 6.0     # in-panel statistics
NC_LEGEND_PT <- 6.0

# ggplot2 linewidth is in mm; 1 pt = 0.3528 mm.
pt_lw <- function(pt) pt * 0.3528
# Nature Communications: "The thinnest lines in the final figure should be no smaller
# than one point wide." Panels are rendered at final size, so 1 pt here IS 1 pt on page.
NC_AXIS_LW  <- pt_lw(1.0)
NC_BOX_LW   <- pt_lw(1.0)
NC_GRID_LW  <- pt_lw(1.0)

# ggplot2 geom_text size is in mm of cap height: pt / .pt
mm_text <- function(pt) pt / ggplot2::.pt

# Species short forms follow the manuscript's own Results wording:
#   "R. faecis, L. eligens, KLE1615, Coprococcus eutactus and
#    Lachnospiraceae bacterium AM48_27BH"
# Binomials are italic; bare strain designations are roman.
SPECIES_SHORT <- c(
  "Roseburia faecis"                    = "*R. faecis*",
  "Lachnospira eligens"                 = "*L. eligens*",
  "Clostridiales bacterium KLE1615"     = "KLE1615",
  "Coprococcus eutactus"                = "*C. eutactus*",
  "Lachnospiraceae bacterium AM48 27BH" = "AM48 27BH"
)
short_species <- function(x) unname(ifelse(x %in% names(SPECIES_SHORT),
                                           SPECIES_SHORT[x], x))

COL_HEALTHY <- "#1B9E8F"
COL_CANCER  <- "#C43C3C"

# Three group palettes, one per kind of contrast. All share the tonal register of the
# Healthy/Cancer pair (S 50-54%, L 36-50%) so the figures read as one system; only the
# hue changes, so no two contrasts are ever encoded by the same colour.
#   disease status   Healthy / Cancer   green  / red     H146 / H0
#   ICI response     R / NR             blue   / amber   H208 / H28
#   marker split     Low / High         violet, light / dark (ordinal, one hue)
COL_R    <- "#2E5F8A"   # H208 S50 L36 - matches COL_HEALTHY's tone
COL_NR   <- "#C47B3B"   # H28  S54 L50 - matches COL_CANCER's tone
COL_LOW  <- "#9868B1"   # H280 S32 L55
COL_HIGH <- "#642F7F"   # H280 S46 L34
GROUP_COLORS    <- c(Healthy = COL_HEALTHY, Cancer = COL_CANCER)
RESPONSE_COLORS <- c(R = COL_R, NR = COL_NR)
SPLIT_COLORS    <- c(Low = COL_LOW, High = COL_HIGH)
COL_POINT   <- "#202020"   # neutral dark: gives contrast over ANY group fill

# No in-panel titles: panel letters and titles are set in the figure layout, not by ggplot.
# No grey strip boxes: the strip label alone carries the facet identity.
theme_nc <- function(base_pt = NC_TICK_PT, grid_y = FALSE) {
  th <- ggplot2::theme_classic(base_family = "Arial", base_size = base_pt) +
    ggplot2::theme(
      plot.title       = ggplot2::element_blank(),
      plot.subtitle    = ggplot2::element_blank(),
      plot.caption     = ggplot2::element_blank(),
      axis.title       = ggplot2::element_text(family = "Arial", size = NC_TITLE_PT,
                                               colour = "black"),
      axis.text        = ggplot2::element_text(family = "Arial", size = base_pt,
                                               colour = "black"),
      axis.line        = ggplot2::element_line(linewidth = NC_AXIS_LW, colour = "black"),
      axis.ticks       = ggplot2::element_line(linewidth = NC_AXIS_LW, colour = "black"),
      axis.ticks.length = ggplot2::unit(1.4, "pt"),
      legend.title     = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(family = "Arial", size = NC_LEGEND_PT),
      legend.key.size  = ggplot2::unit(7, "pt"),
      legend.margin    = ggplot2::margin(0, 0, 0, 0),
      strip.text       = ggplot2::element_text(family = "Arial", size = NC_STRIP_PT,
                                               colour = "black",
                                               margin = ggplot2::margin(b = 1.5, t = 0.5)),
      strip.background = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing    = ggplot2::unit(3.5, "pt"),
      plot.margin      = ggplot2::margin(2, 3, 1, 1)
    )
  if (grid_y) {
    th <- th + ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = NC_GRID_LW),
      panel.grid.major.x = ggplot2::element_blank()
    )
  } else {
    th <- th + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  }
  th
}

# Box + points with the summary drawn ON TOP of the data.
# Layer 1: filled box, no outline emphasis, so the group colour reads as an area.
# Layer 2: neutral dark points -- a same-colour point over a same-colour fill has no
#          contrast no matter what alpha is used, which is why the points are not
#          mapped to Group.
# Layer 3: the same box redrawn with fill = NA, putting median/hinges/whiskers back
#          above the point cloud without hiding it.
box_points_layers <- function(width = 0.58, pt_size = 0.45, pt_alpha = 0.30,
                              jitter_w = 0.13, seed = 260723) {
  list(
    # First layer is fill only: colour = NA. Drawing the outline here as well would
    # render it twice and darken it (caught in review, 2026-08-19).
    ggplot2::geom_boxplot(width = width, outlier.shape = NA, alpha = 1,
                          linewidth = NC_BOX_LW, colour = NA),
    ggplot2::geom_point(shape = 16, colour = COL_POINT, size = pt_size, alpha = pt_alpha,
                        position = ggplot2::position_jitter(width = jitter_w, height = 0,
                                                            seed = seed)),
    ggplot2::geom_boxplot(width = width, outlier.shape = NA, fill = NA,
                          linewidth = NC_BOX_LW, colour = "#2D2D2D")
  )
}

# Render at the FINAL placed size. width/height are the on-page inches.
save_nc <- function(plot, path, width_in, height_in, dpi = 600) {
  if (file.exists(path)) unlink(path)
  ggplot2::ggsave(path, plot, width = width_in, height = height_in,
                  dpi = dpi, bg = "white", device = ragg::agg_png)
  stopifnot(file.exists(path), file.info(path)$size > 5000)
  invisible(path)
}
