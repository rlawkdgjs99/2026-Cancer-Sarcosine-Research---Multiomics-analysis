#!/usr/bin/env Rscript
# Render the saved per-cohort statistics. Fixed label layout never moves dots.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg);library(jsonlite)})
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
OUT<-dirname(dirname(normalizePath(sub('^--file=','',arg[1]))))
QA<-file.path(OUT,'qa');ss<-fread(file.path(OUT,'cohort_summary.csv'))
COL<-c('R higher'='#1B9E8F','NR higher'='#C43C3C','Below display thresholds'='#BFC2C4')
pub<-theme_classic(base_family='Arial',base_size=12)+theme(
  plot.title=element_text(size=17,face='bold',margin=margin(b=7)),
  plot.subtitle=element_text(size=10.5,colour='#50565A',margin=margin(b=13)),
  axis.title=element_text(size=13,colour='black'),axis.text=element_text(colour='black',size=11),
  axis.line=element_line(linewidth=.45),axis.ticks=element_line(linewidth=.4),
  legend.position='top',legend.justification='left',legend.title=element_blank(),
  legend.text=element_text(size=10),legend.margin=margin(0,0,4,0),
  plot.caption=element_text(size=9,colour='#50565A',hjust=0,margin=margin(t=12)),
  plot.margin=margin(12,15,12,12))
saveplot<-function(g,dest,stem,width,height) {
  for(ext in c('png','pdf')) {
    tmp<-tempfile(fileext=paste0('.',ext));target<-file.path(dest,paste0(stem,'.',ext))
    if(ext=='png')ggsave(tmp,g,width=width,height=height,units='in',dpi=600,bg='white',device=ragg::agg_png)
    else ggsave(tmp,g,width=width,height=height,units='in',bg='white',device=cairo_pdf)
    stopifnot(file.copy(tmp,target,overwrite=TRUE),file.info(target)$size>1000);unlink(tmp)
  }
}
checks<-list();positions<-list()
for(i in seq_len(nrow(ss))) {
  co<-ss$Cohort[i];s<-ss[i];dest<-file.path(OUT,co)
  d<-fread(file.path(dest,paste0(co,'_species_DA_results.csv')))
  b<-fread(file.path(dest,paste0(co,'_species_DA_bar_data.csv')))
  d[,Display_class:=factor(Display_class,levels=names(COL))]
  b[,Display_class:=factor(Display_class,levels=names(COL))]
  b[,Species_display:=factor(Species_display,levels=Species_display)]
  stopifnot(nrow(d)==s$tested_species,nrow(b)==s$bar_R+s$bar_NR,
    all(b$BH_q<.05 & abs(b$log2FC_NR_vs_R)>1))
  labels<-d[volcano_label==TRUE]
  labels[,label:=gsub('_',' ',Species)]
  labels[,c('label_x','label_y','hjust'):=list(NA_real_,NA_real_,NA_real_)]
  # Explicit taxonomy keys prevent q-order changes from silently remapping labels.
  if(co=='PRJNA751792') {
    stopifnot(identical(labels$Species,'Roseburia_intestinalis'))
    labels[,c('label_x','label_y','hjust'):=list(-3.0,1.65,1)]
  }
  if(co=='PRJNA1023797') {
    layout<-data.table(
      Species=c('Lachnospiraceae_bacterium','GGB9619_SGB15067','Faecalibacillus_faecis',
        'Eubacterium_ventriosum','Faecalibacillus_intestinalis','GGB3005_SGB3996',
        'Methylobacterium_SGB15164','Streptococcus_anginosus','Streptococcus_parasanguinis',
        'Streptococcus_gordonii'),
      label_x=c(rep(-3.2,7),4,4,2.5),
      label_y=c(2.8,2.55,1.05,2.30,2.05,1.80,1.55,2.65,1.75,1.4),
      hjust=c(rep(1,7),0,0,0))
    stopifnot(setequal(labels$Species,layout$Species))
    jj<-match(labels$Species,layout$Species)
    labels[,c('label_x','label_y','hjust'):=list(layout$label_x[jj],layout$label_y[jj],layout$hjust[jj])]
  }
  if(co=='PRJEB22863')stopifnot(nrow(labels)==0L)
  stopifnot(all(is.finite(labels$label_x)),all(is.finite(labels$label_y)))
  labels[,line_end_x:=label_x+ifelse(hjust==1,.04,-.04)]
  positions[[co]]<-labels[,.(Cohort,Species_full,log2FC_NR_vs_R,neg_log10_BH_q,
    label,label_x,label_y,hjust)]
  # Keep common y limits. PRJEB22863 requires a wider x range to retain extreme
  # mean ratios involving zeros; no statistical points are clipped.
  xlims<-if(co=='PRJEB22863')c(-17.5,9.5) else c(-6.8,8.2)
  xbreaks<-if(co=='PRJEB22863')seq(-15,5,5) else seq(-6,6,2)
  stopifnot(all(d$log2FC_NR_vs_R>xlims[1] & d$log2FC_NR_vs_R<xlims[2]),
    all(d$neg_log10_BH_q>=0 & d$neg_log10_BH_q<3))
  subtitle<-sprintf('%s · %d R / %d NR Run records',co,s$n_R,s$n_NR)
  caption<-sprintf('Within-cohort Wilcoxon test; BH across %d species. Colored: q < 0.05 and |log₂FC| > 1.',s$tested_species)
  v<-ggplot(d,aes(log2FC_NR_vs_R,neg_log10_BH_q))+
    geom_vline(xintercept=c(-1,1),linetype='dashed',colour='#A5AAAD',linewidth=.35)+
    geom_hline(yintercept=-log10(.05),linetype='dashed',colour='#A5AAAD',linewidth=.35)+
    geom_segment(data=labels,aes(xend=line_end_x,yend=label_y),colour='#9BA2A5',linewidth=.3)+
    geom_point(data=d[Display_class=='Below display thresholds'],aes(colour=Display_class),alpha=.75,size=1.5,show.legend=FALSE)+
    geom_point(data=d[Display_class!='Below display thresholds'],aes(colour=Display_class),alpha=.9,size=2.05,show.legend=TRUE)+
    geom_text(data=labels,aes(x=label_x,y=label_y,label=label,colour=Display_class,hjust=hjust),
      family='Arial',fontface='italic',size=3.1,show.legend=FALSE)+
    scale_colour_manual(values=COL,drop=FALSE,labels=c('R higher','NR higher','Below cutoffs'))+
    scale_x_continuous(breaks=xbreaks,expand=expansion(mult=c(.025,.025)))+
    scale_y_continuous(breaks=seq(0,3,.5),limits=c(0,3),expand=expansion(mult=c(0,.03)))+
    coord_cartesian(xlim=xlims,clip='off')+
    labs(title='NSCLC species differential abundance',subtitle=subtitle,
      x=expression(log[2]~'fold change (NR/R)'),y=expression(-log[10]~'(BH q)'),caption=caption)+pub
  if(co=='PRJEB22863')v<-v+annotate('text',x=-4,y=1.8,label='No species with BH q < 0.05',
    family='Arial',size=3.6,colour='#666D70')
  barcap<-sprintf('q < 0.05 and |log₂FC| > 1. Up to 15 per direction by q; %d R-higher and %d NR-higher qualify.',
    s$R_display_hits,s$NR_display_hits)
  if(nrow(b)>0L) {
    p<-ggplot(b,aes(log2FC_NR_vs_R,Species_display,fill=Display_class))+
      geom_vline(xintercept=0,colour='#535A5E',linewidth=.4)+
      geom_col(width=if(nrow(b)==1L).5 else .76,colour=NA,show.legend=TRUE)+
      scale_fill_manual(values=COL[1:2],drop=FALSE)+
      scale_y_discrete(expand=expansion(add=if(nrow(b)==1L)1.2 else .55))+
      theme(axis.text.y=element_text(face='italic',size=11.5))
  } else {
    p<-ggplot()+
      annotate('text',x=.5,y=1,label='No species meet\nq < 0.05 and |log₂FC| > 1',
        family='Arial',size=4.1,colour='#50565A')+
      scale_y_continuous(limits=c(.5,1.5),breaks=NULL)
  }
  p<-p+scale_x_continuous(breaks=-3:4,limits=c(-3,4),expand=expansion(mult=c(.02,.02)))+
    labs(title='NSCLC differential species',subtitle=subtitle,
      x=expression(log[2]~'fold change (NR/R)'),y=NULL,caption=barcap)+pub+
    theme(axis.text.y=element_text(face='italic',size=11.5),
      axis.line.y=element_blank(),axis.ticks.y=element_blank())
  if(nrow(b)==0L)p<-p+theme(legend.position='none')
  height<-if(nrow(b)<=1L)4 else 6
  saveplot(v,dest,paste0(co,'_species_DA_volcano'),10,6.5)
  saveplot(p,dest,paste0(co,'_species_DA_bar'),9,height)
  vb<-ggplot_build(v);bb<-ggplot_build(p)
  stopifnot(nrow(vb$data[[4]])==sum(d$Display_class=='Below display thresholds'),
    nrow(vb$data[[5]])==sum(d$Display_class!='Below display thresholds'),
    nrow(vb$data[[6]])==nrow(labels))
  vd<-rbindlist(lapply(vb$data[c(4,5)],function(z)if(nrow(z))as.data.table(z[,c('x','y')]) else NULL))
  observed<-data.table(x=vd$x,y=vd$y);expected<-d[,.(x=log2FC_NR_vs_R,y=neg_log10_BH_q)]
  setorder(observed,x,y);setorder(expected,x,y)
  stopifnot(max(abs(as.matrix(observed)-as.matrix(expected)))<1e-12)
  if(nrow(b))stopifnot(nrow(bb$data[[2]])==nrow(b),
    max(abs(sort(bb$data[[2]]$x)-sort(b$log2FC_NR_vs_R)))<1e-12)
  checks[[co]]<-list(tested=nrow(d),gray=sum(d$Display_class=='Below display thresholds'),
    colored=sum(d$Display_class!='Below display thresholds'),labels=nrow(labels),
    bars=nrow(b),all_plotted_coordinates_match_CSV=TRUE,x_limits=xlims,y_limits=c(0,3),
    png_volcano_pixels=c(6000,3900),png_bar_pixels=c(5400,as.integer(height*600)),
    dpi=600,bar_empty_annotation=nrow(b)==0L)
}
fwrite(rbindlist(positions),file.path(QA,'volcano_label_positions.csv'))
write_json(checks,file.path(QA,'plot_verification.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Saved all six plots as PNG and vector PDF; plot-data checks passed.\n')

