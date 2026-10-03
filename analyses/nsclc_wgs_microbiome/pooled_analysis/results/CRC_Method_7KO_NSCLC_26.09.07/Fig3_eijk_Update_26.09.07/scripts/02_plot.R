#!/usr/bin/env Rscript
# Visualization only. All displayed numerical values read from 01_analysis.R.
suppressPackageStartupMessages({library(ggplot2);library(ggrepel);library(patchwork);library(data.table);library(ragg);library(svglite)})
set.seed(42); options(warn=1)
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE);out<-dirname(dirname(normalizePath(sub('^--file=','',arg[1]))))
d<-readRDS(file.path(out,'tables/plot_objects.rds'))
BLUE<-'#2D628D';ORANGE<-'#C67B37';TEAL<-'#199D90';INK<-'#26343D';GREY<-'#AEB6BB';PURPLE<-'#8848A4'
base<-theme_classic(base_size=10,base_family='Arial')+theme(
 text=element_text(colour=INK),axis.text=element_text(colour=INK,size=9),axis.title=element_text(size=10),
 axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.3),
 plot.title=element_text(size=14,face='bold',margin=margin(b=5)),
 plot.subtitle=element_text(size=9.5,colour='#64717B',margin=margin(b=9)),
 plot.caption=element_text(size=7.5,colour='#64717B',hjust=0,lineheight=1.15,margin=margin(t=8)),
 legend.title=element_blank(),legend.position='top',legend.text=element_text(size=9),
 plot.title.position='plot',plot.caption.position='plot',plot.margin=margin(9,11,9,10))
fmt<-function(p)ifelse(p<.001,formatC(p,format='e',digits=1),sprintf('%.3f',p))
save<-function(p,stem,w,h){
 # Scientific data plots generated from tables; no AI image synthesis.
 f<-file.path(out,'figures',stem)
 ragg::agg_png(paste0(f,'.png'),width=w,height=h,units='mm',res=600,background='white');print(p);dev.off()
 # macOS Quartz embeds native fonts without the optional XQuartz/Cairo libraries.
 if (capabilities('aqua')) grDevices::quartz(type='pdf',file=paste0(f,'.pdf'),width=w/25.4,height=h/25.4,family='Arial')
 else grDevices::cairo_pdf(paste0(f,'.pdf'),width=w/25.4,height=h/25.4,family='Arial')
 print(p);dev.off()
 svglite::svglite(paste0(f,'.svg'),width=w/25.4,height=h/25.4,bg='white');print(p);dev.off()
}
# e: All observations, including zero abundances. Same untransformed sums as f.
eplots<-lapply(seq_len(nrow(d$e_stats)),function(i){
 s<-d$e_stats$ID[i];z<-d$scores;z$Group<-factor(z$Group,levels=c('NR','R'));z$value<-z[[s]]
 title<-c('Degradation','Production','Production / degradation')[i]
 ylab<-if(i<3)'Summed KO relative abundance (%)' else expression(log[2]*' ratio')
 label<-paste0('P = ',fmt(d$e_stats$wilcox_p[i]),'\nBH q = ',fmt(d$e_stats$q[i]))
 yr<-range(z$value);span<-diff(yr);ann<-yr[2]+span*.16
 ggplot(z,aes(Group,value,colour=Group,fill=Group))+
  geom_boxplot(width=.48,outlier.shape=NA,alpha=.17,linewidth=.42,colour=INK)+
  geom_point(position=position_jitter(width=.18,height=0,seed=42),size=.43,alpha=.5,stroke=0)+
  annotate('segment',x=1,xend=2,y=yr[2]+span*.04,yend=yr[2]+span*.04,linewidth=.32)+
  annotate('text',x=1.5,y=ann,label=label,size=2.7,lineheight=1.1,colour=INK)+
  scale_colour_manual(values=c(NR=ORANGE,R=BLUE))+scale_fill_manual(values=c(NR=ORANGE,R=BLUE))+
  scale_y_continuous(expand=expansion(mult=c(.04,.04)),limits=c(if(i<3)0 else yr[1]-span*.04,yr[2]+span*.30))+
  labs(title=title,x=NULL,y=ylab)+base+theme(legend.position='none',plot.title=element_text(size=9.7),
  plot.title.position='panel',axis.title.y=element_text(size=8.4),axis.text=element_text(size=8.5),plot.margin=margin(6,7,5,3))
})
pe<-wrap_plots(eplots,nrow=1)+plot_annotation(
 title='Sarcosine pathway scores',subtitle='Three NSCLC cohorts pooled · 392 NR / 432 R Run records',
 caption='Five degradation KOs and two production KOs. Boxes: median and interquartile range; whiskers: 1.5 × IQR.\nAll observations shown. Two-sided Wilcoxon; BH across three scores. Ratio includes the same 1e-8 pseudocount as Figure 3f.',
 theme=base)
save(pe,'Fig3e_Pathway_Scores_7KO',174,81)
# i: Keep every selected taxon; labels use the frozen six-per-response rule.
zi<-as.data.table(d$i_data);labs_i<-zi[label==TRUE]
labs_i[,display:=gsub('_',' ',species)]
labs_i[,`:=`(anchor_y=NA_real_,anchor_x=NA_real_)]
# Anchors in clear upper/lower bands; leader lines never alter point coordinates.
for(g in c('NR','R')){
 ids<-which(labs_i$enriched==g);ids<-ids[order(-labs_i$Production_rho[ids])]
 labs_i$anchor_y[ids]<-seq(.42,-.36,length.out=length(ids))
 labs_i$anchor_x[ids]<-if(g=='NR')-.38 else .37
}
zi[,`:=`(display='',nudge_x=0,nudge_y=0)]
li<-match(labs_i$species,zi$species)
zi$display[li]<-labs_i$display
zi$nudge_x[li]<-labs_i$anchor_x-labs_i$Degradation_rho
zi$nudge_y[li]<-labs_i$anchor_y-labs_i$Production_rho
it<-d$i_tests
pi<-ggplot(zi,aes(Degradation_rho,Production_rho,colour=enriched))+
 geom_hline(yintercept=0,colour='#B7BEC2',linetype='dashed',linewidth=.38)+
 geom_vline(xintercept=0,colour='#B7BEC2',linetype='dashed',linewidth=.38)+
 geom_point(size=2,alpha=.85)+
 geom_text_repel(data=zi,aes(label=display),nudge_x=zi$nudge_x,
  nudge_y=zi$nudge_y,direction='y',size=2.65,fontface='italic',seed=42,
  box.padding=.45,point.padding=.5,min.segment.length=0,segment.colour='#9BA6AE',segment.size=.25,
  max.overlaps=Inf,max.iter=20000,max.time=10,show.legend=FALSE)+
 scale_colour_manual(values=c(R=BLUE,NR=ORANGE),breaks=c('R','NR'),labels=c('R-enriched','NR-enriched'))+
 coord_cartesian(xlim=c(-.55,.59),ylim=c(-.44,.47),clip='off')+
 labs(title='Species–function associations',subtitle=sprintf('%d response-associated species · five degradation / two production KOs',nrow(zi)),
  x=expression('Degradation association (meta '*rho*')'),y=expression('Production association (meta '*rho*')'),
  caption=paste0('Species-level R-enriched vs NR-enriched comparison: degradation BH q = ',fmt(it$q_two_comparisons[1]),
   '; production BH q = ',fmt(it$q_two_comparisons[2]),'.\nUnadjusted within-cohort Spearman correlations; Fisher-z DL meta-analysis.\nPoints are species, not patients. Associations do not establish species-resolved KO carriage.'))+base+
  theme(plot.caption=element_text(size=7.2),plot.margin=margin(9,12,9,12))
save(pi,'Fig3i_Species_Function_7KO',168,133)
# j: Recomputed full eligible-species set, preserving old fixed focal labels.
zj<-as.data.table(d$j_data);js<-d$j_stat
zj[,significance:=factor(significance,levels=c('Neither','CRC only','NSCLC only','Both'))]
focal<-c('Lachnospira_eligens','Roseburia_faecis','Clostridium_sp_AF36_4','Clostridiales_bacterium_KLE1615')
jlab<-zj[species %in% focal];jlab[,label:=c(Lachnospira_eligens='L. eligens',Roseburia_faecis='R. faecis',Clostridium_sp_AF36_4='Clostridium sp. AF36_4',Clostridiales_bacterium_KLE1615='KLE1615')[species]]
pj<-ggplot(zj,aes(CRC_rho,NSCLC_rho))+
 geom_hline(yintercept=0,colour='#BCC2C5',linewidth=.35)+geom_vline(xintercept=0,colour='#BCC2C5',linewidth=.35)+
 geom_abline(slope=1,intercept=0,linetype='dashed',colour='#ABB4BA',linewidth=.45)+
 geom_point(aes(fill=significance),shape=21,colour='white',stroke=.3,size=2.05,alpha=.9)+
 geom_text_repel(data=jlab,aes(label=label),size=3.15,fontface='italic',seed=42,box.padding=.65,
  point.padding=.4,min.segment.length=0,max.overlaps=Inf,segment.colour='#8B969E',max.iter=10000)+
 annotate('label',x=.10,y=-.43,hjust=0,label=sprintf('Shared species: %d\nSpearman ρ = %.3f\nSame direction: %.1f%%\nPermutation P = %.4f',nrow(zj),js$spearman_rho,100*js$sign_concordance,js$conditional_permutation_p),
  size=3.1,linewidth=0,fill='white',colour=INK,lineheight=1.2)+
 scale_fill_manual(values=c('Neither'='#C9CDD0','CRC only'=TEAL,'NSCLC only'='#4E92C1','Both'=PURPLE),drop=FALSE)+
 coord_equal(xlim=c(-.65,.67),ylim=c(-.63,.67),expand=FALSE)+
 labs(title='Cross-cancer species–degradation concordance',subtitle='The same five-KO degradation score in CRC and NSCLC',
  x=expression('CRC association (partial Spearman meta '*rho*')'),y=expression('NSCLC association (partial Spearman meta '*rho*')'),
  caption='Colours: BH q < 0.05 within disease. Phenotype-adjusted correlations; DL meta-analysis.\n1,000 within-cohort/phenotype permutations. NSCLC cross-project overlap retained; P is conditional.\nCRC uses the original paired-library BioSample aggregation. Association, not taxonomic KO attribution.')+
 base+theme(legend.position='bottom',plot.title=element_text(size=12.3),plot.caption=element_text(size=7.4))
save(pj,'Fig3j_Cross_Cancer_Concordance_7KO',148,153)
# k: Actual recomputed membership, including the uniform-FDR sensitivity.
zk<-as.data.table(d$k_data);ks<-d$k_stats[1,];kb<-d$k_stats[2,]
zk<-zk[order(-shared,-in_CRC,species)];zk[,row:=rev(seq_len(.N))]
short<-c(Roseburia_faecis='R. faecis',Lachnospira_eligens='L. eligens',Clostridiales_bacterium_KLE1615='KLE1615',
 Coprococcus_eutactus='C. eutactus',Lachnospiraceae_bacterium_AM48_27BH='AM48 27BH',Clostridium_sp_AF36_4='Clostridium sp. AF36_4')
zk[,display:=ifelse(species %in% names(short),short[species],gsub('_',' ',species))]
axis_labels<-vapply(zk$display,function(x)if(x %in% c('L. eligens','R. faecis','C. eutactus'))
 paste0("italic('",x,"')") else if(x=='Clostridium sp. AF36_4')"italic('Clostridium')~'sp. AF36_4'" else paste0("'",x,"'"),'')
members<-rbind(zk[in_CRC==TRUE,.(species,row,x=1,cohort='CRC')],zk[in_NSCLC==TRUE,.(species,row,x=2,cohort='NSCLC')])
pk<-ggplot()+geom_rect(data=zk[shared==TRUE],aes(xmin=.68,xmax=2.32,ymin=row-.46,ymax=row+.46),fill='#F1F4F6')+
 geom_segment(data=zk[shared==TRUE],aes(x=1,xend=2,y=row,yend=row),colour='#AAB4BA',linewidth=.8)+
 geom_point(data=members,aes(x,row,fill=cohort),shape=21,size=5,stroke=.5,colour=INK)+
 scale_fill_manual(values=c(CRC=TEAL,NSCLC=BLUE),guide='none')+
 scale_x_continuous(breaks=c(1,2),labels=c('CRC\n4 cohorts','NSCLC\n3 cohorts'),position='top',limits=c(.65,2.35),expand=c(0,0))+
 scale_y_continuous(breaks=zk$row,labels=parse(text=axis_labels),limits=c(.4,max(zk$row)+.6),expand=c(0,0))+
 labs(title='Shared degradation-associated species',subtitle=sprintf('%d of %d NSCLC recurrent species shared with CRC',ks$overlap,ks$n_NSCLC),
  x=NULL,y=NULL,caption=sprintf('All cohorts: ρ > 0.3 and prevalence ≥10%%. CRC: BH q < 0.05; NSCLC: nominal P < 0.05.\nUniform per-cohort BH criterion: %d NSCLC core species. Recurrence is exploratory.',kb$n_NSCLC))+
 base+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(size=11,face='bold'),
  axis.text.y=element_text(size=11),plot.title=element_text(size=13),plot.subtitle=element_text(size=10),plot.caption=element_text(size=7.6))
save(pk,'Fig3k_Shared_Species_7KO',152,112)
saveRDS(list(e=pe,i=pi,j=pj,k=pk),file.path(out,'qa/plot_builds.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'qa/plot_sessionInfo.txt'))
message('Four scientific panels exported as PNG/PDF/SVG.')
