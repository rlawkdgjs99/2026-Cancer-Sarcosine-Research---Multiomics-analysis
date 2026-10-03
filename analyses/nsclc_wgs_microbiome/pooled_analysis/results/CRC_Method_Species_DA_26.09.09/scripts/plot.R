#!/usr/bin/env Rscript
# Draw the author-requested volcano and horizontal DA bars from frozen results.
# No statistics or feature-selection rules are changed here.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg);library(jsonlite)})
set.seed(42)
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
OUT<-dirname(dirname(normalizePath(sub('^--file=','',arg[1]))))
d<-fread(file.path(OUT,'NSCLC_species_DA_results.csv'))
b<-fread(file.path(OUT,'NSCLC_species_DA_bar_data.csv'))
stopifnot(nrow(d)==541L,all(is.finite(d$log2FC_NR_vs_R)),all(is.finite(d$neg_log10_BH_q)),
  all(b$BH_q<.05 & abs(b$log2FC_NR_vs_R)>1),!anyDuplicated(b$Species_full))
COL<-c('R higher'='#1B9E8F','NR higher'='#C43C3C','Below display thresholds'='#BFC2C4')
d[,Display_class:=factor(Display_class,levels=names(COL))]
b[,Display_class:=factor(Display_class,levels=names(COL))]
b[,Species_display:=factor(Species_display,levels=Species_display)]
labels<-d[volcano_label==TRUE]
labels[,label:=gsub('_',' ',Species)]
# Wrap long taxonomy labels at spaces without abbreviating the scientific ID.
labels[,label:=vapply(label,function(s)paste(strwrap(s,width=30),collapse='\n'),character(1))]
# Fixed label coordinates keep text/point separation identical in PNG and PDF.
# Only text is positioned; all dot coordinates remain the saved statistics.
left<-labels[log2FC_NR_vs_R<0]
right<-labels[log2FC_NR_vs_R>0]
stopifnot(nrow(left)==8L,nrow(right)==7L)
left[,`:=`(label_x=-2.9,label_y=c(4.08,3.7,3.22,2.82,2.51,2.2,1.88,1.52),hjust=1)]
right[,`:=`(label_x=c(3.1,3.75,2.1,1.9,2.25,3.1,3.3),
  label_y=c(4.05,2.65,2.95,3.5,1.95,1.6,2.23),hjust=0)]
labels<-rbind(left,right)
labels[,line_end_x:=label_x+ifelse(hjust==1,.04,-.04)]
pub<-theme_classic(base_family='Arial',base_size=12)+theme(
  plot.title=element_text(size=17,face='bold',margin=margin(b=7)),
  plot.subtitle=element_text(size=10.5,colour='#50565A',margin=margin(b=13)),
  axis.title=element_text(size=13,colour='black'),axis.text=element_text(colour='black',size=11),
  axis.line=element_line(linewidth=.45),axis.ticks=element_line(linewidth=.4),
  legend.position='top',legend.justification='left',legend.title=element_blank(),
  legend.text=element_text(size=10),legend.margin=margin(0,0,4,0),
  plot.caption=element_text(size=9,colour='#50565A',hjust=0,margin=margin(t=12)),
  plot.margin=margin(12,15,12,12))
v<-ggplot(d,aes(log2FC_NR_vs_R,neg_log10_BH_q))+
  geom_vline(xintercept=c(-1,1),linetype='dashed',colour='#A5AAAD',linewidth=.35)+
  geom_hline(yintercept=-log10(.05),linetype='dashed',colour='#A5AAAD',linewidth=.35)+
  geom_segment(data=labels,aes(xend=line_end_x,yend=label_y),
    colour='#9BA2A5',linewidth=.3)+
  geom_point(data=d[Display_class=='Below display thresholds'],aes(colour=Display_class),alpha=.75,size=1.5)+
  geom_point(data=d[Display_class!='Below display thresholds'],aes(colour=Display_class),alpha=.9,size=2.05)+
  geom_text(data=labels,aes(x=label_x,y=label_y,label=label,colour=Display_class,hjust=hjust),
    family='Arial',fontface='italic',size=3.1,lineheight=1.0,show.legend=FALSE)+
  scale_colour_manual(values=COL,drop=FALSE,labels=c('R higher','NR higher','Below cutoffs'))+
  scale_x_continuous(breaks=-4:4,expand=expansion(mult=c(.10,.10)))+
  scale_y_continuous(breaks=0:5,limits=c(0,5.2),expand=expansion(mult=c(0,.015)))+
  coord_cartesian(xlim=c(-5.6,6.2),clip='off')+
  labs(title='NSCLC species differential abundance',
    subtitle='Three pooled cohorts · 432 R / 392 NR Run records',
    x=expression(log[2]~'fold change (NR/R)'),y=expression(-log[10]~'(BH q)'),
    caption='Pooled Wilcoxon test; BH across 541 species. Colored: q < 0.05 and |log₂FC| > 1.')+pub

# Each bar is one species' mean-abundance log2 ratio, not a patient mean with SE.
p<-ggplot(b,aes(log2FC_NR_vs_R,Species_display,fill=Display_class))+
  geom_vline(xintercept=0,colour='#535A5E',linewidth=.4)+
  geom_col(width=.76,colour=NA)+
  scale_fill_manual(values=COL[1:2],drop=FALSE)+
  scale_x_continuous(breaks=seq(-3,4,1),limits=c(-2.65,3.6),expand=expansion(mult=c(.01,.02)))+
  scale_y_discrete(expand=expansion(add=.55))+
  labs(title='NSCLC differential species',
    subtitle='Three pooled cohorts · 432 R / 392 NR Run records',
    x=expression(log[2]~'fold change (NR/R)'),y=NULL,
    caption='q < 0.05 and |log₂FC| > 1. Up to 15 per direction by q; 14 R-higher and 11 NR-higher qualify.')+
  pub+theme(axis.text.y=element_text(face='italic',size=11.5),
    axis.line.y=element_blank(),axis.ticks.y=element_blank(),
    legend.position='top',plot.margin=margin(12,15,12,12))
saveplot<-function(g,name,width,height){
  for(ext in c('png','pdf')) {
    tmp<-tempfile(fileext=paste0('.',ext));dest<-file.path(OUT,paste0(name,'.',ext))
    if(ext=='png')ggsave(tmp,g,width=width,height=height,units='in',dpi=600,bg='white',device=ragg::agg_png)
    else ggsave(tmp,g,width=width,height=height,units='in',bg='white',device=cairo_pdf)
    stopifnot(file.copy(tmp,dest,overwrite=TRUE),file.info(dest)$size>1000);unlink(tmp)
  }
}
saveplot(v,'NSCLC_species_DA_volcano',10.0,6.5)
saveplot(p,'NSCLC_species_DA_bar',9.0,8.1)
write_json(list(volcano=list(width_inches=10,height_inches=6.5),bar=list(width_inches=9,height_inches=8.1),
  dpi=600,colours=as.list(COL),background='white',labels_seed=42,
  display='same fixed CRC thresholds; no automatic biological selection changes'),
  file.path(OUT,'qa','render_spec.json'),auto_unbox=TRUE,pretty=TRUE)
fwrite(labels[,.(Species_full,log2FC_NR_vs_R,neg_log10_BH_q,label,label_x,label_y,hjust)],
  file.path(OUT,'qa','volcano_label_positions.csv'))
vb<-ggplot_build(v);bb<-ggplot_build(p)
stopifnot(nrow(vb$data[[4]])==516L,nrow(vb$data[[5]])==25L,
  nrow(bb$data[[2]])==nrow(b),
  max(abs(sort(bb$data[[2]]$x)-sort(b$log2FC_NR_vs_R)))<1e-12)
write_json(list(volcano_grey_points=516,colored_points=25,bar_rows=nrow(b),
  all_bar_coordinates_match_CSV=TRUE),file.path(OUT,'qa','plot_data_verification.json'),
  auto_unbox=TRUE,pretty=TRUE)
cat('Saved PNG and vector PDF for both plots.\n')
