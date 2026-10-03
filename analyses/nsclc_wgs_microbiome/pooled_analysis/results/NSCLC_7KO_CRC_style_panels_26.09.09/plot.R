#!/usr/bin/env Rscript
# Publication figures from analysis.R outputs only; no statistical reselection.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(ggvenn);library(ragg)})
script<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]));OUT<-dirname(script)
obj<-readRDS(file.path(OUT,'plot_data.rds'));list2env(obj,envir=environment())
FIG<-file.path(OUT,'figures');dir.create(FIG,showWarnings=FALSE)
TEAL<-'#1B9E8F';RED<-'#C43C3C';GREY<-'#9D9D9D';ROLE<-c(Degradation='#2879C9',Production='#C89E00',Both='#AD39C8')
theme_set(theme_classic(base_size=11,base_family='Arial'))
base_theme<-theme(plot.title=element_text(face='bold',size=14),plot.subtitle=element_text(size=10,color='#444444'),
 plot.caption=element_text(size=8,color='#444444',hjust=0),axis.text=element_text(color='black'),
 plot.margin=margin(12,16,12,12),legend.title=element_text(size=10))
qfmt<-function(q)ifelse(is.na(q),'NA',ifelse(q<.001,formatC(q,format='e',digits=1),formatC(q,format='f',digits=3)))
disp<-function(s)gsub('_',' ',s,fixed=TRUE)
savefig<-function(p,name,w,h){
 ggsave(file.path(FIG,paste0(name,'.png')),p,width=w,height=h,dpi=300,bg='white',device=ragg::agg_png,limitsize=FALSE)
 ggsave(file.path(FIG,paste0(name,'.pdf')),p,width=w,height=h,bg='white',device=grDevices::cairo_pdf,limitsize=FALSE)
 cat('Saved',name,'\n');flush.console()
}
forest_plot<-function(k,caption=TRUE){
 d<-copy(forest[KO==k]);d[,y:=match(Cohort,rev(CO))]
 yl<-paste0(d$Cohort,'\nR ',d$n_R,' / NR ',d$n_NR)
 p<-ggplot(d,aes(y=y))+
  geom_rect(data=d[y%%2==0],aes(xmin=-.65,xmax=.65,ymin=y-.45,ymax=y+.45),inherit.aes=FALSE,fill='#F2F2F2')+
  geom_vline(xintercept=0,linetype='dotted',linewidth=.6)+
  geom_segment(data=d[estimable==TRUE],aes(x=CI95_low,xend=CI95_high,yend=y),linewidth=.55,color='#333333')+
  geom_point(data=d[estimable==TRUE],aes(x=effect_NR_minus_R,color=color_hex),size=3.2)+
  geom_text(data=d[estimable==FALSE],aes(x=0,label='All zero — not estimable'),size=3,color='#777777')+
  annotate('text',x=-.34,y=3.7,label='R higher',color=TEAL,size=3.7)+
  annotate('text',x=.34,y=3.7,label='NR higher',color=RED,size=3.7)+
  scale_color_identity()+scale_y_continuous(breaks=d$y,labels=yl,limits=c(.5,4),expand=c(0,0))+
  scale_x_continuous(breaks=c(-.5,0,.5),limits=c(-.65,.65),expand=c(0,0))+
  labs(title=paste0(d$Gene[1],' (',k,')'),x=expression('Rank-biserial effect size ('*r[rb]*')'),y=NULL)+
  base_theme+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=9))
 if(caption)p<-p+labs(caption='Point: P(NR > R) − P(NR < R); bar: bootstrap 95% CI.\nTeal/red: BH q < 0.05; grey: q ≥ 0.05. BH within KO across estimable cohorts.\nPoint size is fixed; all estimable points are filled.')
 p
}
savefig(forest_plot('K00303'),'Main_soxB_3cohort_forest',6.7,3.9)
savefig(forest_plot('K08688'),'Main_Creatinase_3cohort_forest',6.7,3.9)
other<-KOS[!KOS%in%c('K00303','K08688')]
p_other<-wrap_plots(lapply(other,forest_plot,caption=FALSE),ncol=2)+
 plot_annotation(title='Other five sarcosine KOs | 3 NSCLC cohorts',
 caption='Point: P(NR > R) − P(NR < R); bar: bootstrap 95% CI. Teal/red: BH q < 0.05; grey: q ≥ 0.05.\nBH within each KO across its estimable cohort contrasts. All-zero contrasts have no estimate, CI or test.\nAll circles have the same size and are filled.',theme=base_theme)
savefig(p_other,'Supplement_other5KO_3cohort_forests',13,10.5)
d<-copy(pooled);d[,y:=8-seq_len(.N)];d[,label:=paste0(Gene,' (',KO,')\n',Role,' · prevalence ',sprintf('%.1f',prevalence_percent),'%')]
d[,qtext:=paste0('q = ',qfmt(BH_q))]
p_pool<-ggplot(d,aes(y=y))+
 geom_rect(data=d[y%%2==0],aes(xmin=-2.35,xmax=2,ymin=y-.43,ymax=y+.43),inherit.aes=FALSE,fill='#F2F2F2')+
 geom_vline(xintercept=0,linewidth=.45)+
 geom_segment(aes(x=0,xend=log2FC_NR_vs_R,yend=y,color=color_hex),linewidth=.7)+
 geom_point(aes(x=log2FC_NR_vs_R,fill=color_hex,size=prevalence_percent),shape=21,stroke=.5,color='#555555')+
 geom_text(aes(x=log2FC_NR_vs_R+ifelse(log2FC_NR_vs_R<0,-.14,.14),label=qtext,hjust=ifelse(log2FC_NR_vs_R<0,1,0)),size=3)+
 annotate('text',x=-1.2,y=7.75,label='R higher',color=TEAL,size=4)+annotate('text',x=1,y=7.75,label='NR higher',color=RED,size=4)+
 scale_color_identity()+scale_fill_identity()+scale_size_area(max_size=10,limits=c(0,100),breaks=c(1,10,50,100),name='Prevalence (%)')+
 scale_x_continuous(limits=c(-2.35,2),breaks=-2:2,expand=c(0,0))+
 scale_y_continuous(breaks=d$y,labels=d$label,limits=c(.5,8.15),expand=c(0,0))+
 labs(title='Seven sarcosine KOs | pooled NSCLC',subtitle='824 Runs: R 432 / NR 392',
 x=expression(log[2]*' fold change (NR / R)'),y=NULL,
 caption='Teal/red: BH q < 0.05; grey: q ≥ 0.05. BH across all seven KOs. Circle area = pooled prevalence.\nAll seven targets retained, including prevalence <10%; mean-abundance fold change uses pseudocount 10⁻⁸ (%).')+
 base_theme+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.position='right',axis.text.y=element_text(size=10))
savefig(p_pool,'Supplement_7KO_pooled_lollipop',11,6.5)
da_plot<-function(a,compact=FALSE){
 all<-copy(associated_DA[Analysis==a]);d<-all[Plot_bar==TRUE];setorder(d,-log2FC_NR_vs_R);d[,y:=rev(seq_len(.N))]
 ttl<-if(a=='Pooled')'Sarcosine-associated species | pooled NSCLC' else a
 sub<-paste0('R ',all$n_R[1],' / NR ',all$n_NR[1],' Runs; ',nrow(d),' of ',nrow(all),' candidates pass BH q < 0.05')
 if(nrow(d)==0)return(ggplot()+annotate('text',x=.5,y=.5,label='No species pass BH q < 0.05',size=4,color='#555555')+
   xlim(0,1)+ylim(0,1)+labs(title=ttl,subtitle=sub,caption=paste0('BH across ',all$BH_family_n[1],' estimable candidates; same pooled candidate set.'))+theme_void(base_family='Arial')+base_theme+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank()))
 # Separate role squares in a narrow aligned strip, avoiding axis/data space.
 span<-max(abs(associated_DA[Plot_bar==TRUE]$log2FC_NR_vs_R))
 lim<-max(1,ceiling(span*1.25*2)/2)
 p<-ggplot(d,aes(y=y))+
  geom_vline(xintercept=0,linewidth=.45)+geom_col(aes(x=log2FC_NR_vs_R,fill=bar_color_hex),width=.65,orientation='y',color='#333333',linewidth=.3)+
  geom_text(aes(x=log2FC_NR_vs_R+ifelse(log2FC_NR_vs_R<0,-.025,.025)*lim,label=paste0('q=',qfmt(BH_q_selected)),
   hjust=ifelse(log2FC_NR_vs_R<0,1,0)),size=if(compact)2.8 else 3.1)+
  annotate('text',x=-lim*.52,y=nrow(d)+.9,label='R-enriched',color=TEAL,size=3.6)+
  annotate('text',x=lim*.52,y=nrow(d)+.9,label='NR-enriched',color=RED,size=3.6)+
  scale_fill_identity()+scale_x_continuous(limits=c(-lim,lim),breaks=c(-lim,0,lim),expand=c(0,0))+
  scale_y_continuous(breaks=d$y,labels=disp(d$Species),limits=c(.4,nrow(d)+1.45),expand=c(0,0))+
  labs(title=ttl,subtitle=sub,x=expression(log[2]*' fold change (NR / R)'),y=NULL,
   caption=paste0('BH across ',all$BH_family_n[1],' candidates; no fold-change cutoff. Role: KO correlation, not demonstrated metabolism.'))+
  base_theme+theme(axis.text.y=element_text(face='italic',size=if(compact)9 else 11),axis.line.y=element_blank(),axis.ticks.y=element_blank())
 # All selected NSCLC candidates have one qualifying degradation KO in these data.
 # Add a role square beside each species through a separate annotation column.
 p<-p+geom_point(aes(x=-lim*.97,color=Role),shape=15,size=2.7)+
  scale_color_manual(values=ROLE,drop=FALSE,name='Associated KO pathway')+
  theme(legend.position='bottom',legend.text=element_text(size=9))
 p
}
savefig(da_plot('Pooled'),'Main_sarcosine_species_DA_pooled',11,6.4)
for(co in CO)savefig(da_plot(co),paste0('Supplement_species_DA_',co),10.5,4.2)
p_da<-wrap_plots(lapply(CO,da_plot,compact=TRUE),ncol=1)+plot_annotation(title='Sarcosine-associated species | individual NSCLC cohorts',
 caption='Species selected once from the pooled seven-KO correlation screen; the same candidate set is tested in each cohort.\nDA bars use BH within this candidate set; this is distinct from the full-species DA testing used for the enrichment Venn.',theme=base_theme)
savefig(p_da,'Supplement_species_DA_3cohorts',11,11.5)
cats<-unique(venn$Category)
venn_plot<-function(cat){
 d<-venn[Category==cat];sets<-setNames(lapply(CO,function(co)d[[co]]),CO)
 # Fixed circles with explicit exact-region counts, including zero-sized sets.
 # Area is schematic and never used to imply set size.
 labels<-c('A','B','C')
 t<-seq(0,2*pi,length.out=400);centers<-data.table(Set=labels,cx=c(-.55,.55,0),cy=c(.35,.35,-.55))
 circles<-rbindlist(lapply(1:3,function(i)data.table(Set=labels[i],x=centers$cx[i]+cos(t),y=centers$cy[i]+sin(t))))
 pos<-data.table(A=c(T,F,F,T,T,F,T),B=c(F,T,F,T,F,T,T),C=c(F,F,T,F,T,T,T),
  x=c(-1.03,1.03,0,0,-.55,.55,0),y=c(.55,.55,-1.15,.9,-.42,-.42,.03))
 pos[,N:=vapply(seq_len(.N),function(i)sum(d[[CO[1]]]==A[i]&d[[CO[2]]]==B[i]&d[[CO[3]]]==C[i]),integer(1))]
 labs<-data.table(x=c(-1,1,0),y=c(1.63,1.63,-1.86),label=paste0(CO,'\nn = ',vapply(CO,function(co)sum(d[[co]]),integer(1))))
 p<-ggplot()+geom_polygon(data=circles,aes(x,y,group=Set,fill=Set),alpha=.27,color='#777777',linewidth=.5)+
  geom_text(data=pos,aes(x,y,label=N),size=5)+geom_text(data=labs,aes(x,y,label=label),size=3.1,lineheight=1)+
  scale_fill_manual(values=c(A='#E69F00',B='#56B4E9',C='#009E73'))+coord_fixed(xlim=c(-1.9,1.9),ylim=c(-2.2,2),clip='off')+
  labs(title=cat)+theme_void(base_family='Arial')+base_theme+theme(legend.position='none',plot.title=element_text(size=13,hjust=.5),axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
 p
}
p_venn<-wrap_plots(lapply(cats,venn_plot),ncol=2)+plot_annotation(title='Species overlap across three NSCLC cohorts',
 caption='Production/degradation-associated: Spearman ρ > 0.3 and BH q < 0.05 with the corresponding 2-KO/5-KO score.\nNR/R-enriched: full-species DA BH q < 0.05 and log₂FC > 1 / < −1. Species prevalence ≥10% within cohort.\nBH separately within each cohort and test family. Circles are schematic; numbers are exact species counts. No three-cohort intersections.',theme=base_theme)
savefig(p_venn,'Supplement_species_4criteria_3cohort_Venn',11.5,11)
stopifnot(nrow(selected)>0)
sl<-rbindlist(lapply(seq_len(nrow(selected)),function(si)rbindlist(lapply(seq_along(KOS),function(j){
 data.table(KO=KOS[j],Species_full=selected$Species_full[si],Species=selected$Species[si],
  KO_abundance_percent=K[,j],Species_abundance_percent=S[,selected$Species_full[si]],Group=meta$Group,Run_ID=meta$Run_ID)
}))))
sl[,KOlabel:=factor(paste0(genes[match(KO,KOS)],' (',KO,')'),levels=paste0(genes,' (',KOS,')'))]
sl[,Specieslabel:=factor(vapply(disp(Species),function(x)paste(strwrap(x,width=25),collapse='\n'),''),
 levels=vapply(disp(selected$Species),function(x)paste(strwrap(x,width=25),collapse='\n'),''))]
stat<-copy(ko_cor[Analysis=='Pooled' & Species_full%in%selected$Species_full])
stat[,KOlabel:=factor(paste0(genes[match(KO,KOS)],' (',KO,')'),levels=levels(sl$KOlabel))]
stat[,Specieslabel:=factor(vapply(disp(Species),function(x)paste(strwrap(x,width=25),collapse='\n'),''),levels=levels(sl$Specieslabel))]
stat[,label:=paste0('ρ = ',sprintf('%.2f',rho),'\nq = ',qfmt(BH_q))]
# facet_grid shares each KO column x scale and each species row y scale if
# transposed; requested CRC layout is KO rows/species columns, so draw an
# explicit patchwork grid with consistent limits per KO and per species.
cells<-list()
for(j in seq_along(KOS))for(si in seq_len(nrow(selected))){
 d<-sl[KO==KOS[j]&Species_full==selected$Species_full[si]];st<-stat[KO==KOS[j]&Species_full==selected$Species_full[si]]
 p<-ggplot(d,aes(KO_abundance_percent,Species_abundance_percent,color=Group))+
  geom_point(size=.7,alpha=.55,stroke=0)+
  annotate('text',x=Inf,y=Inf,label=st$label,hjust=1.05,vjust=1.15,size=2.8,color='black',lineheight=1.05)+
  scale_color_manual(values=c(R=TEAL,NR=RED),breaks=c('R','NR'))+
  scale_x_continuous(limits=c(0,max(d$KO_abundance_percent)*1.05),breaks=scales::breaks_pretty(3),expand=expansion(mult=c(.02,.04)))+
  scale_y_continuous(limits=c(0,max(d$Species_abundance_percent)*1.15),breaks=scales::breaks_pretty(3),expand=expansion(mult=c(.02,.02)))+
  labs(title=if(j==1)as.character(d$Specieslabel[1]) else NULL,
   x=paste0(genes[j],' (',KOS[j],') (%)'),y=if(si==1)'Species abundance (%)' else NULL)+
  theme_classic(base_size=8,base_family='Arial')+theme(legend.position='none',axis.text=element_text(color='black',size=7),
   axis.title.x=element_text(size=8),axis.title.y=element_text(size=8),plot.title=element_text(face='italic',size=10,hjust=.5),
   plot.margin=margin(4,6,6,6))
 cells[[length(cells)+1]]<-p
}
p_scatter<-wrap_plots(cells,ncol=nrow(selected))+plot_annotation(title=paste0('Seven sarcosine KOs × ',nrow(selected),' NSCLC-selected species'),
 subtitle='Pooled 824 Runs (R 432 / NR 392) | teal: R; red: NR | raw relative abundance (%)',
 caption='Each panel shows all 824 Runs, including zero abundances. Spearman ρ; BH q from the full 7 KO × 541 species family (3,787 tests).\nSpecies selected by KO association (|ρ| > 0.3, BH q < 0.05) plus global species DA BH q < 0.05; all seven KO pairs are shown.\nExploratory pooled correlations; no cohort adjustment. WGS KO abundance is not gene expression. Runs are not verified unique patients.',theme=base_theme)
savefig(p_scatter,'Supplement_7KO_selected_species_scatter',16.5,19.5)
