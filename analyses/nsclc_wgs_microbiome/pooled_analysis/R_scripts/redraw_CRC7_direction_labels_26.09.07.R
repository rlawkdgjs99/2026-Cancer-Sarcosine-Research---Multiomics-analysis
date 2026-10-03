#!/usr/bin/env Rscript
# DISPLAY ONLY: use frozen tables and the hash-locked original plotting block.
# Never source the original analysis script or rerun its statistical procedures.
suppressPackageStartupMessages({library(ggplot2); library(ggtext); library(ggrepel);
  library(patchwork); library(ragg); library(digest)})
options(stringsAsFactors=FALSE)
arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
script <- normalizePath(sub('^--file=', '', arg))
pool <- dirname(dirname(script))
old <- file.path(pool, 'results', 'CRC_Method_7KO_NSCLC_26.09.07')
out <- file.path(old, 'Direction_Labels')
original <- file.path(pool, 'R_scripts', 'CRC_method_seven_KO_NSCLC_26.09.07.R')
hash <- function(p) digest(file=p, algo='sha256', serialize=FALSE)
stopifnot(hash(original)=='06664a11a0ebbf0d810262365cee5e23ed81b06f3c4608d26e3f7ba3965e1a2b')
read <- function(p) read.csv(p, check.names=FALSE, fileEncoding='UTF-8-BOM')
inputs <- read(file.path(old, 'qa', 'input_sha256.csv'))
protected <- unique(c(original, inputs$path,
  list.files(file.path(old, 'tables'), full.names=TRUE),
  list.files(file.path(old, 'figures'), full.names=TRUE)))
stopifnot(all(file.exists(protected)))
before <- vapply(protected, hash, '')
stopifnot(identical(unname(vapply(inputs$path, hash, '')), inputs$SHA256))
dir.create(out, showWarnings=FALSE)
for(d in c('figures','qa')) dir.create(file.path(out,d), showWarnings=FALSE)
ko <- read(file.path(old,'tables','NSCLC_seven_KO_pooled_statistics.csv'))
pa <- read(file.path(old,'tables','NSCLC_pathway_pooled_statistics.csv'))
cross <- read(file.path(old,'tables','CRC_NSCLC_pooled_KO_concordance.csv'))
concordance <- read(file.path(old,'tables','concordance_summary.csv'))
target <- c('K00301','K00302','K00303','K00305','K00306','K00315','K08688')
gene <- c('Sarcosine oxidase','soxA','soxB','soxG','PIPOX','DMGDH','Creatinase')
role <- c(rep('Degradation',5),rep('Production',2))
bg <- cross[cross$genomewide_background,]
stopifnot(identical(ko$ID,target), nrow(pa)==3L, nrow(cross)==6982L,
  nrow(bg)==6978L, sum(cross$targeted)==7L,
  all(ko$n_R==432L), all(ko$n_NR==392L))

# Parse the known file. Evaluate ONLY the expression block from `cols <- ...`
# through the DEFINITION of `save_png <- function(...)`. The analysis, its
# input-reading/statistical calls, and its output-writing calls are not evaluated.
expr <- parse(file=original, keep.source=FALSE)
lhs <- vapply(expr, function(e) {
  if(is.call(e) && identical(e[[1]],as.name('<-')) && is.symbol(e[[2]]))
    as.character(e[[2]]) else ''
}, '')
start <- which(lhs=='cols'); end <- which(lhs=='save_png')
stopifnot(length(start)==1L, length(end)==1L, start<end)
set.seed(42)
for(i in seq.int(start,end)) eval(expr[[i]])

# Capture every original ggplot layer, including points and CIs.
baseline <- lapply(list(KO=forest, Scores=pp, Concordance=fig3),
  function(p) ggplot_build(p)$data)
forest <- forest +
  annotate('text', x=-xf*.97, y=7.65, label='← Higher in NR',
    hjust=0, family='Arial', fontface='bold', size=3, colour=cols[['NR higher']]) +
  annotate('text', x=xf*.97, y=7.65, label='Higher in R →',
    hjust=1, family='Arial', fontface='bold', size=3, colour=cols[['R higher']])
pp <- pp +
  annotate('text', x=-xx*.97, y=3.65, label='← Higher in NR',
    hjust=0, family='Arial', fontface='bold', size=3, colour=cols[['NR higher']]) +
  annotate('text', x=xx*.97, y=3.65, label='Higher in R →',
    hjust=1, family='Arial', fontface='bold', size=3, colour=cols[['R higher']])
fig1 <- forest+right+plot_layout(widths=c(4,1.25))+plot_annotation(
  title='Individual microbial sarcosine genes',
  subtitle='Three NSCLC cohorts pooled · 432 R / 392 NR Run records',
  caption='95% CI: 5,000 bootstrap resamples. BH correction across 7 KOs.\nOpen circles: pooled prevalence <10%. Unadjusted Run-level comparison.\nR: responders. NR: non-responders.', theme=theme_plot())
fig2 <- pp+pr+plot_layout(widths=c(4,1.5))+plot_annotation(
  title='Sarcosine pathway balance',
  subtitle='Three NSCLC cohorts pooled',
  caption='95% bootstrap CI. Two-sided Wilcoxon P. BH correction across 3 scores.\nRatio = log₂[(production + 10⁻⁸)/(degradation + 10⁻⁸)]. Unadjusted Run-level comparison.\nR: responders. NR: non-responders.', theme=theme_plot())

# Secondary axes are identity copies for direction labels, NOT new quantities.
# Neutral label colour is intentional: scatter point colours encode KO role.
# Give the seven LABELS deterministic non-overlapping positions. Original KO
# data points stay fixed. Replace only the original text-repel annotation layer.
stopifnot(inherits(fig3$layers[[6]]$geom,'GeomTextRepel'))
fig3$layers <- fig3$layers[-6]
hi$label_x <- c(.06,-.12,.34,-.26,.22,-.22,-.29)
hi$label_y <- c(.32,.36,.13,.08,-.07,.27,.19)
hi$label_hjust <- c(0,1,0,1,0,1,1)
fig3 <- fig3 +
  geom_segment(data=hi,aes(x=CRC_effect,y=NSCLC_effect,xend=label_x,yend=label_y),
    inherit.aes=FALSE,colour='#71808A',linewidth=.3) +
  geom_text(data=hi,aes(x=label_x,y=label_y,label=paste0(gene,' (',KO,')'),hjust=label_hjust),
    inherit.aes=FALSE,family='Arial',size=2.9,colour='#2D3D47')
fig3 <- fig3 +
  scale_x_continuous(sec.axis=dup_axis(breaks=c(-.52*li,.52*li),
    labels=c('← Higher in Cancer','Higher in Healthy →'), name=NULL)) +
  scale_y_continuous(sec.axis=dup_axis(breaks=c(-.52*li,.52*li),
    labels=c('Higher in NR ↓','Higher in R ↑'), name=NULL)) +
  theme(axis.text.x.top=element_text(face='bold',size=8,colour='#293944',margin=margin(b=6)),
    axis.text.y.right=element_text(face='bold',size=8,angle=0,colour='#293944',margin=margin(l=7)),
    axis.line.x.top=element_blank(),axis.line.y.right=element_blank(),
    axis.ticks.x.top=element_blank(),axis.ticks.y.right=element_blank()) +
  labs(subtitle='Seven sarcosine KOs highlighted',
    caption=paste0('Blue: degradation. Orange: production. Open circles: prevalence <10% in either disease.\n',
    'Background and concordance: KOs with pooled prevalence ≥10% in both diseases.\n',
    'Unadjusted Run-level comparisons. Low-prevalence targets retained separately.\n',
    'R: responders. NR: non-responders. Group-direction labels do not encode KO role.'))

# Provenance stays in analysis records, not the NSCLC display subtitles.
stopifnot(!grepl('CRC',fig1$patches$annotation$subtitle,fixed=TRUE),
  !grepl('CRC',fig2$patches$annotation$subtitle,fixed=TRUE),
  !grepl('CRC',fig3$labels$subtitle,fixed=TRUE),
  identical(fig3$labels$x,'CRC pooled rank-biserial effect (Healthy − Cancer)'))

current <- lapply(list(KO=forest,Scores=pp,Concordance=fig3), function(p) ggplot_build(p)$data)
checks <- list()
for(n in names(baseline)) for(i in seq_along(baseline[[n]])) {
  # Only the seven text-label positions and their leader lines were changed.
  # All other original layers must be identical, not merely numerically close.
  if(n=='Concordance' && i==6L) next
  j <- if(n=='Concordance' && i==7L) 6L else i
  ok <- identical(baseline[[n]][[i]],current[[n]][[j]])
  stopifnot(ok)
  checks[[length(checks)+1L]] <- data.frame(plot=n,original_layer=i,
    rows=nrow(baseline[[n]][[i]]),entire_layer_identical=ok)
}
stopifnot(identical(current$Concordance[[7]]$x,hi$CRC_effect),
  identical(current$Concordance[[7]]$y,hi$NSCLC_effect),
  identical(current$Concordance[[8]]$label,paste0(hi$gene,' (',hi$KO,')')))
save_png(fig1,'NSCLC_CRC7_Individual_KOs.png',6.6,4.1)
save_png(fig2,'NSCLC_CRC7_Pathway_Balance.png',6.6,3.05)
save_png(fig3,'CRC_NSCLC_CRC7_Concordance.png',7.1,6.6)
after <- vapply(protected,hash,'')
stopifnot(identical(before,after))
write.csv(data.frame(path=protected,sha256_before=before,sha256_after=after,
  unchanged=before==after),file.path(out,'qa','protected_file_hashes.csv'),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(out,'qa','plot_layer_validation.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'qa','sessionInfo.txt'))
writeLines(c('# Direction-label redraw — display only',
  'Frozen results loaded from the parent tables directory; no statistical analysis rerun.',
  'Both forest plots: negative/left/orange = Higher in NR; positive/right/blue = Higher in R.',
  'Scatter: CRC negative/left = Higher in Cancer; positive/right = Higher in Healthy.',
  'Scatter: NSCLC negative/bottom = Higher in NR; positive/top = Higher in R.',
  'Scatter colours remain KO roles (blue degradation, orange production), not response groups.',
  'All 13 retained original built plot layers, including every point/CI coordinate, are identical.',
  'Only the scatter text-repel layer was replaced with seven manually placed labels and leader lines.',
  'Every leader line starts at its original KO data coordinate; label strings are unchanged.',
  'Canvas/annotation layout changed; physical pixel positions are not claimed identical.',
  'These are still unadjusted pooled Run-level exploratory candidates, NOT adjusted meta-analysis.',
  'The known Sample-name overlap and manuscript adoption questions are unchanged.',
  'No original inputs, tables, plots, manuscript, or author PPT were overwritten.'),
  file.path(out,'DISPLAY_NOTES.md'))
cat('PASS:',length(checks),'original plot layers unchanged;',length(protected),
  'protected files unchanged. Three PNGs rendered in',out,'\n')
