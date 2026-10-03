# Figure 6e — CD8와 conventional DC를 각각 별도 그림으로
# Display only: NES/q는 frozen 표에서 그대로. 재계산 없음.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr); library(grid)})
ROOT <- "."   # workspace root laid out as described in LAYOUT.md (working dir)
source(file.path(ROOT,"공공_Metabolomics&Metagenomics_분석모음/_shared/theme_nc_26.08.18.R"))
SRC <- file.path(ROOT,"2024_Drug_Res_Updates_NSCLC/RNA-seq공공데이터_GSE207422/analysis_sarcosine_FINAL_26.08.25",
                 "15_AllCell_Degradation_Pathway_Overview_26.09.14/results/tables/Fig6_scRNA_GSEA_heatmap.csv")
OUTDIR <- commandArgs(TRUE)[1]
d <- read_csv(SRC, show_col_types = FALSE)

FILL_SIG <- "#C0574C"; FILL_NS <- "#FFFFFF"; EDGE <- "#C0574C"

panel <- function(comp, progs, title, sub, foot, w, h, file) {
  x <- d %>% filter(compartment==comp, program %in% progs) %>%
    mutate(program=factor(program, levels=rev(progs)),
           sig = BH_q < 0.05,
           lab = ifelse(sig, sprintf("%.2f*", NES), sprintf("%.2f", NES)),
           qlab= ifelse(BH_q < 1e-3, sprintf("q = %.1e", BH_q), sprintf("q = %.3f", BH_q)))
  stopifnot(nrow(x)==length(progs))
  WRAP <- c("TCR signaling"="TCR signaling",
            "TNFR2-related noncanonical NF-\u03BAB"="TNFR2-related\nnoncanonical NF-\u03BAB",
            "IFN-\u03B3 response"="IFN-\u03B3 response",
            "Antigen cross-presentation (MHC-I)"="Antigen cross-\npresentation (MHC-I)")
  x <- x %>% mutate(plab = factor(WRAP[as.character(program)],
                                 levels=unname(WRAP[rev(progs)])))
  XMAX <- 4.30
  p <- ggplot(x, aes(y=plab, x=NES)) +
    geom_col(aes(fill=sig, colour=sig), width=.56, linewidth=pt_lw(0.9)) +
    geom_text(aes(x=NES+0.10, label=lab), hjust=0, family="Arial",
              size=mm_text(6.6), fontface="bold", colour="#101010") +
    geom_text(aes(x=NES+0.10, label=qlab), hjust=0, vjust=2.35, family="Arial",
              size=mm_text(5.0), colour="#5A5A5A") +
    scale_fill_manual(values=c(`TRUE`=FILL_SIG, `FALSE`=FILL_NS), guide="none") +
    scale_colour_manual(values=c(`TRUE`=FILL_SIG, `FALSE`=EDGE), guide="none") +
    scale_x_continuous(limits=c(0,XMAX), breaks=c(0,1,2,3), expand=c(0,0)) +
    scale_y_discrete(expand=expansion(add=.52)) +
    labs(x="NES (Degradation-High vs Low)", y=NULL) +
    theme_nc() +
    theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(),
          axis.text.y=element_text(size=6.4, colour="black", hjust=1, lineheight=1.02),
          axis.text.x=element_text(size=5.8), axis.title.x=element_text(size=6.2, margin=margin(t=1.4)),
          panel.grid.major.x=element_line(colour="#EDEDED", linewidth=pt_lw(0.5)),
          plot.margin=margin(t=9.5, r=2.0, b=5.2, l=2.0, unit="mm"))
  ragg::agg_png(file.path(OUTDIR,file), width=w, height=h, units="mm", res=600, background="white")
  grid.newpage(); grid.draw(ggplotGrob(p))
  grid.text(title, x=unit(2.0,"mm"), y=unit(1,"npc")-unit(2.4,"mm"), hjust=0, vjust=1,
            gp=gpar(fontfamily="Arial", fontsize=7.4, fontface="bold"))
  grid.text(sub, x=unit(2.0,"mm"), y=unit(1,"npc")-unit(5.5,"mm"), hjust=0, vjust=1,
            gp=gpar(fontfamily="Arial", fontsize=5.3, col="#5A5A5A"))
  grid.text(foot, x=unit(2.0,"mm"), y=unit(1.5,"mm"), hjust=0, vjust=0,
            gp=gpar(fontfamily="Arial", fontsize=4.6, col="#7A7A7A"))
  invisible(dev.off())
  cat(file, "|", nrow(x), "programmes\n"); print(x %>% select(program, NES, BH_q))
}

panel("CD8",
      c("TCR signaling","TNFR2-related noncanonical NF-κB","IFN-γ response"),
      "CD8⁺ T cells",
      "Patient pseudobulk \u00B7 12 patients (High 7 / Low 5) \u00B7 16,229 cells",
      "* BH q < 0.05 within the 5,959 CD8 gene sets tested. Open bar: not significant.",
      74, 46, "CD8.png")

panel("cDC",
      c("Antigen cross-presentation (MHC-I)"),
      "Conventional DCs",
      "Patient pseudobulk \u00B7 8 patients (High 5 / Low 3) \u00B7 772 cells",
      "* BH q < 0.05 within the 3,797 conventional-DC gene sets tested.",
      74, 28, "DC.png")

panel("Whole_tumour",
      c("Antigen cross-presentation (MHC-I)","TCR signaling",
        "TNFR2-related noncanonical NF-κB","IFN-γ response"),
      "Whole-tumour pseudobulk",
      "Patient pseudobulk · 12 patients (High 7 / Low 5) · 78,192 cells",
      "* BH q < 0.05 within the 4,958 whole-tumour gene sets tested.",
      74, 56, "Pseudobulk.png")
