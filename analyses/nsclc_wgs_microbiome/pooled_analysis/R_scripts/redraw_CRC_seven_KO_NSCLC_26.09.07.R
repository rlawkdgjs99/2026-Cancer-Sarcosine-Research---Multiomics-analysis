#!/usr/bin/env Rscript
# AUTHOR REJECTED this display-only interpretation on 2026-09-07.
# Historical draft only. Do not rerun or install as the requested final figure.
# Display-only change to the AUTHOR-selected CRC seven-KO target list.
# Preserve all original estimates, CI, BH families, filters and observation units.
# Missing eligible meta-estimates remain NA, never a plotted zero.
suppressPackageStartupMessages({
  library(ggplot2); library(ggtext); library(patchwork); library(ggrepel)
  library(digest); library(metafor)
})
options(stringsAsFactors = FALSE)
set.seed(42)
args <- commandArgs(FALSE)
script <- normalizePath(sub('^--file=', '', grep('^--file=', args, value=TRUE)))
pooled <- dirname(dirname(script)); base <- dirname(pooled)
out <- file.path(pooled, 'results', 'CRC_Seven_KO_Redraw_26.09.07')
dir.create(out, showWarnings=FALSE)
oldforest <- file.path(pooled, 'results', 'Fig2i_individual_sarcosine_gene_meta_26.08.23')
oldcross <- file.path(pooled, 'results', 'Fig2_genomewide_KO_cross_disease_concordance_CANDIDATE_26.08.23')
cohort_dirs <- c(PRJNA751792='NSCLC_PRJNA751792', PRJNA1023797='NSCLC_PRJNA1023797', PRJEB22863='NSCLC_RCC_PRJEB22863')
da_paths <- file.path(base, cohort_dirs, 'results', paste0('sarcosine_KO_DA_',names(cohort_dirs),'.csv'))
inputs <- unique(c(file.path(pooled,'results','pooled_sarcosine_KO_meta.csv'),
  file.path(base,'sarcosine_KO_set.csv'), da_paths,
  list.files(oldforest,full.names=TRUE), list.files(oldcross,full.names=TRUE),
  file.path(base,cohort_dirs,'R_scripts','run_sarcosine.R'),
  file.path(pooled,'R_scripts',c('pooled_sarcosine.R','104_Fig2i_individual_sarcosine_gene_meta_26.08.23.R','106_candidate_Fig2_genomewide_KO_cross_disease_concordance_26.08.23.R'))))
stopifnot(all(file.exists(inputs)))
hash <- function(x) digest(file=x,algo='sha256',serialize=FALSE)
before <- vapply(inputs,hash,'')
write.csv(data.frame(path=inputs,SHA256=before), file.path(out,'input_sha256.csv'),row.names=FALSE)
read <- function(p) read.csv(p,check.names=FALSE)
meta <- read(inputs[1]); annotation <- read(file.path(base,'sarcosine_KO_set.csv'))
cross <- read(file.path(oldcross,'cross_disease_common_KOs_primary.csv'))
stats <- read(file.path(oldcross,'global_concordance_statistics.csv'))
stats <- stats[stats$mode=='biological_sample_primary',,drop=FALSE]
target <- c('K00301','K00302','K00303','K00305','K00306','K00315','K08688')
names_gene <- c('Sarcosine oxidase','soxA','soxB','soxG','PIPOX','DMGDH','Creatinase')
roles <- c(rep('Degradation',5),rep('Production',2))
stopifnot(nrow(meta)==7,nrow(cross)==7400,!anyDuplicated(cross$KO),nrow(stats)==1,
  all(target %in% annotation$KO),
  identical(tolower(roles),annotation$role[match(target,annotation$KO)]))
das <- lapply(da_paths,read)
eligible_counts <- vapply(target,function(k) sum(vapply(das,function(d)
  sum(d$feature==k & is.finite(d$coef) & is.finite(d$stderr) & d$stderr>0),integer(1))),integer(1))
stopifnot(identical(unname(eligible_counts),c(0L,3L,3L,0L,1L,0L,3L)))
f <- data.frame(KO=target,gene=names_gene,role=roles,k_available=eligible_counts)
idx <- match(target,meta$KO)
for (column in c('pooled_coef','ci_lb','ci_ub','pval','qval','I2','tau2')) f[[column]] <- meta[[column]][idx]
f$status <- ifelse(is.na(idx),'Not estimable under original criteria','Original meta-estimate retained')
f$reason <- ifelse(f$k_available==0,'No cohort passes original 10% prevalence gate',
  ifelse(f$k_available==1,'Only one eligible cohort; meta-analysis requires at least two',''))
f$BH_family <- 'Original 7 meta-analyzable NSCLC sarcosine KOs (unchanged)'
stopifnot(identical(f$KO[is.finite(f$pooled_coef)],c('K00302','K00303','K08688')))
# Numerical verification only: reproduce the frozen forest from stored cohort estimates.
meta_checks <- lapply(seq_len(nrow(meta)),function(i) {
  d <- do.call(rbind,lapply(das,function(x) x[x$feature==meta$KO[i],,drop=FALSE]))
  m <- rma(yi=d$coef,sei=d$stderr,method='DL')
  difference <- max(abs(c(as.numeric(m$beta)-meta$pooled_coef[i], m$ci.lb-meta$ci_lb[i],
    m$ci.ub-meta$ci_ub[i],m$pval-meta$pval[i],m$I2-meta$I2[i])))
  data.frame(KO=meta$KO[i],max_abs_difference=difference)
})
meta_checks <- do.call(rbind,meta_checks)
stopifnot(max(meta_checks$max_abs_difference)<1e-10,
  max(abs(p.adjust(meta$pval,'BH')-meta$qval))<1e-12,
  abs(cor(cross$crc_effect,cross$nsclc_effect,method='spearman')-stats$spearman_rho)<1e-12,
  sum(sign(cross$crc_effect)==sign(cross$nsclc_effect))==stats$direction_agree_n)
write.csv(meta_checks,file.path(out,'meta_reconstruction_checks.csv'),row.names=FALSE)
write.csv(f,file.path(out,'NSCLC_CRC_seven_KO_forest_source.csv'),row.names=FALSE,na='')
availability <- f[c('KO','gene','role','k_available','status','reason')]
availability$cross_disease_estimable <- target %in% cross$KO
stopifnot(identical(availability$KO[availability$cross_disease_estimable],c('K00302','K00303','K08688')))
write.csv(availability,file.path(out,'seven_KO_availability.csv'),row.names=FALSE)
cross$CRC_seven_KO_target <- cross$KO %in% target
cross$display_gene <- names_gene[match(cross$KO,target)]
cross$display_role <- roles[match(cross$KO,target)]
write.csv(cross,file.path(out,'CRC_NSCLC_concordance_source.csv'),row.names=FALSE,na='')

COL <- c('R higher'='#2E5F8A','NR higher'='#C47B3B')
f$y <- 7:1
f$direction <- ifelse(f$pooled_coef>=0,'R higher','NR higher')
f$axis_label <- paste0(c('Sarcosine oxidase','<i>soxA</i>','<i>soxB</i>','<i>soxG</i>','PIPOX','DMGDH','Creatinase'),
  " <span style='color:#707070'>(",target,')</span>')
f$q_label <- ifelse(is.na(f$qval),'—',sprintf('%.3f',f$qval))
f$I2_label <- ifelse(is.na(f$I2),'—',sprintf('%.1f',f$I2))
yes <- f[is.finite(f$pooled_coef),]; no <- f[!is.finite(f$pooled_coef),]
forest <- ggplot(f,aes(y=y)) +
  geom_vline(xintercept=0,linetype='22',colour='#A3A8AC',linewidth=.4) +
  geom_hline(yintercept=2.5,colour='#DADDE0',linewidth=.35) +
  geom_segment(data=yes,aes(x=ci_lb,xend=ci_ub,yend=y,colour=direction),linewidth=.85) +
  geom_point(data=yes,aes(x=pooled_coef,colour=direction),size=2.6) +
  geom_text(data=no,aes(x=.12,label='Not estimable'),colour='#858B92',size=2.8) +
  scale_colour_manual(values=COL) +
  scale_x_continuous(breaks=c(-.5,0,.5),limits=c(-.6,.85),expand=expansion(mult=0)) +
  scale_y_continuous(breaks=f$y,labels=f$axis_label,limits=c(.5,7.95),expand=expansion(mult=0)) +
  labs(x='Pooled MaAsLin2 coefficient (R − NR)\nNR higher                         R higher',y=NULL) +
  theme_classic(base_family='Arial',base_size=8) +
  theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),
    axis.text.y=ggtext::element_markdown(size=8,colour='#25323A'),axis.text.x=element_text(colour='#25323A'),
    axis.title.x=element_text(size=8,margin=margin(t=6)),plot.margin=margin(4,2,3,0))
right <- ggplot(f,aes(y=y)) +
  geom_hline(yintercept=2.5,colour='#DADDE0',linewidth=.35) +
  geom_text(aes(x=.4,label=k_available),family='Arial',size=2.8) +
  geom_text(aes(x=1.35,label=q_label),family='Arial',size=2.8) +
  geom_text(aes(x=2.25,label=I2_label),family='Arial',size=2.8) +
  annotate('text',x=c(.4,1.35,2.25),y=7.65,label=c('k','q','I² (%)'),family='Arial',size=2.8,fontface='bold') +
  scale_y_continuous(limits=c(.5,7.95),expand=expansion(mult=0)) +
  scale_x_continuous(limits=c(0,2.8),expand=expansion(mult=0)) + theme_void() +
  theme(plot.margin=margin(4,0,3,1))
fp <- forest+right+plot_layout(widths=c(3.8,1.45)) +
  plot_annotation(title='Individual microbial sarcosine genes',
    subtitle='CRC seven-KO set · Three NSCLC discovery cohorts',
    caption=paste('Bars: 95% CI. k: eligible cohorts. Meta-analysis requires k ≥ 2.',
      'q: BH adjustment across the original seven estimable NSCLC KOs, unchanged.',
      'Degradation: K00301–K00306 shown above. Production: K00315 and K08688.',sep='\n'),
    theme=theme(plot.title=element_text(family='Arial',face='bold',size=12,colour='#25323A'),
      plot.subtitle=element_text(family='Arial',size=8.5,colour='#59656D'),
      plot.caption=element_text(family='Arial',size=7,hjust=0,colour='#59656D',lineheight=1.15),
      plot.margin=margin(9,10,9,10)))

lim <- max(.25,ceiling(max(abs(c(cross$crc_effect,cross$nsclc_effect)))*20)/20)
high <- cross[cross$CRC_seven_KO_target,]
scatter <- ggplot(cross,aes(x=crc_effect,y=nsclc_effect)) +
  geom_hline(yintercept=0,colour='#BEC4C8',linewidth=.35,linetype='22') +
  geom_vline(xintercept=0,colour='#BEC4C8',linewidth=.35,linetype='22') +
  geom_abline(slope=1,intercept=0,colour='#D8DDE0',linewidth=.3,linetype='dotted') +
  geom_point(data=cross[!cross$CRC_seven_KO_target,],colour='#8B979F',alpha=.25,size=.65,stroke=0) +
  geom_point(data=high,aes(fill=display_role),shape=21,size=3,colour='white',stroke=.45) +
  geom_text_repel(data=high,aes(label=paste0(display_gene,' (',KO,')')),
    seed=42,nudge_x=c(-.17,.10,-.06),nudge_y=c(-.17,.10,.13),
    colour='#25323A',family='Arial',size=3,segment.colour='#697780',
    min.segment.length=0,box.padding=.4,point.padding=.3,max.overlaps=Inf,show.legend=FALSE) +
  annotate('text',x=-.94*lim,y=.94*lim,hjust=0,vjust=1,family='Arial',size=2.7,
    colour='#34434C',lineheight=1.2,
    label=sprintf('All 7,400 shared KOs\nSpearman ρ = %.3f\nStratified permutation P < 0.001\nDirection agreement = %.1f%%',stats$spearman_rho,stats$direction_agreement_percent)) +
  scale_fill_manual(values=c(Degradation='#2E5F8A',Production='#C47B3B')) +
  coord_equal(xlim=c(-lim,lim),ylim=c(-lim,lim)) +
  labs(title='Cross-disease concordance of microbial KO effects',
    subtitle='CRC seven-KO set: 3 jointly estimable KOs highlighted',
    x='CRC meta rank-biserial effect (Healthy − Cancer)',
    y='NSCLC meta rank-biserial effect (R − NR)',fill=NULL,
    caption=paste('K00301, K00305, K00306 and K00315 have no eligible paired meta-estimates.',
      'Original gates: ≥10% cohort prevalence and ≥2 estimable cohorts per disease.',
      'R, responder; NR, non-responder. Global statistics use all 7,400 KOs.',sep='\n')) +
  theme_classic(base_family='Arial',base_size=9) +
  theme(plot.title=element_text(face='bold',size=11.2,colour='#25323A'),
    plot.subtitle=element_text(size=8.5,colour='#59656D'),
    axis.text=element_text(colour='#34434C'),axis.title=element_text(size=8.5),
    legend.position='top',legend.justification='left',legend.text=element_text(size=8),
    plot.caption=element_text(hjust=0,size=7,colour='#59656D',lineheight=1.15),
    plot.margin=margin(9,10,9,10))

safe_png <- function(plot,name,w,h) {
  seed <- .Random.seed
  temp <- tempfile(fileext='.png')
  ggsave(temp,plot,width=w,height=h,units='in',dpi=600,bg='white',device=ragg::agg_png)
  .Random.seed <<- seed
  dest <- file.path(out,name)
  stopifnot(file.copy(temp,dest,overwrite=TRUE),file.exists(dest),file.info(dest)$size>10000)
  unlink(temp)
}
safe_png(fp,'NSCLC_CRC7_KO_Forest.png',6.4,3.85)
safe_png(scatter,'CRC_NSCLC_CRC7_KO_Concordance.png',5.5,5.7)
stopifnot(identical(before,vapply(inputs,hash,'')))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(c('DISPLAY-ONLY: no new cohort model, threshold or BH family.',
  'Target IDs: K00301 K00302 K00303 K00305 K00306 K00315 K08688',
  'Forest retains original MaAsLin2 coefficients; concordance retains original rank-biserial effects.',
  'Forest eligible cohort counts: 0 3 3 0 1 0 3.',
  'No invented zero points for four missing targets.',
  'Concordance retains all 7400 original coordinates; only highlight membership changed (3 targets).',
  'All original seven forest DL meta-estimates/CI/P/I2 reconstructed within 1e-10; original BH q reproduced.',
  'Global Spearman rho and direction agreement independently reproduced; stored permutation P retained, not rerun.',
  'All input file hashes unchanged. PPT and manuscript were not opened or modified.'),file.path(out,'verification.txt'))
cat('OUTPUT:',out,'\n')
