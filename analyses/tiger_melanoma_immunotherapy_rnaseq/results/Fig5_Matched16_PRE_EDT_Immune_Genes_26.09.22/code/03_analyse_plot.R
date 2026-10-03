options(stringsAsFactors=FALSE);set.seed(20260922)
args<-commandArgs(trailingOnly=TRUE); o<-normalizePath(if(length(args)) args[1] else '.')
library(data.table);library(ggplot2);library(patchwork)
d<-fread(file.path(o,'tables/patient_data_long32.csv'))
metrics<-c('CD8_percent','CYT','SARDH','PIPOX','GNMT','DMGDH')
labels<-c(CD8_percent='CD8+ T-cell estimate',CYT='Cytolytic expression (CYT)',SARDH='SARDH',PIPOX='PIPOX',GNMT='GNMT',DMGDH='DMGDH')
units<-c(CD8_percent='quanTIseq CD8+ T cells (%)',CYT='CYT (log2 geometric mean TPM)',SARDH='Expression (log2[FPKM + 1])',PIPOX='Expression (log2[FPKM + 1])',GNMT='Expression (log2[FPKM + 1])',DMGDH='Expression (log2[FPKM + 1])')
pre<-d[timepoint=='PRE'];edt<-d[timepoint=='EDT'][match(pre$patient_id,patient_id)]
stopifnot(nrow(pre)==16,nrow(edt)==16,identical(pre$patient_id,edt$patient_id))
wide<-pre[,.(patient_id,response,therapy,age,sex,sample_PRE=sample_id)];wide[,sample_EDT:=edt$sample_id]
for(m in metrics){wide[[paste0('PRE_',m)]]<-pre[[m]];wide[[paste0('EDT_',m)]]<-edt[[m]];wide[[paste0('Delta_',m)]]<-edt[[m]]-pre[[m]]}
fwrite(wide,file.path(o,'tables/patient_data_paired16.csv'))
paired<-rbindlist(lapply(c('All','R','NR'),function(g)rbindlist(lapply(metrics,function(m){
 ix<-if(g=='All') rep(TRUE,16) else pre$response==g;x<-pre[[m]][ix];y<-edt[[m]][ix];dl<-y-x
 t<-wilcox.test(y,x,paired=TRUE,exact=FALSE,correct=TRUE,digits.rank=12)
 data.table(metric=m,group=g,n=length(x),PRE_median=median(x),EDT_median=median(y),delta_median=median(dl),delta_mean=mean(dl),delta_SD=sd(dl),n_increased=sum(dl>0),n_decreased=sum(dl<0),n_unchanged=sum(dl==0),V=as.numeric(t$statistic),P=t$p.value)
}))))
paired[,BH_family:=ifelse(group=='All','overall_paired_6','stratified_paired_12')];paired[,BH_q:=p.adjust(P,'BH'),by=BH_family]
delta<-rbindlist(lapply(metrics,function(m){x<-wide[[paste0('Delta_',m)]][wide$response=='R'];y<-wide[[paste0('Delta_',m)]][wide$response=='NR'];t<-wilcox.test(x,y,exact=FALSE,correct=TRUE)
 data.table(metric=m,n_R=length(x),n_NR=length(y),R_delta_median=median(x),NR_delta_median=median(y),R_delta_mean=mean(x),NR_delta_mean=mean(y),U_R=as.numeric(t$statistic),rank_biserial=2*as.numeric(t$statistic)/(length(x)*length(y))-1,P=t$p.value)
}));delta[,BH_q:=p.adjust(P,'BH')];delta[,BH_family:='R_vs_NR_delta_6']
cross<-rbindlist(lapply(c('PRE','EDT'),function(tp)rbindlist(lapply(metrics,function(m){x<-d[timepoint==tp & response=='R'][[m]];y<-d[timepoint==tp & response=='NR'][[m]];t<-wilcox.test(x,y,exact=FALSE,correct=TRUE)
 data.table(metric=m,timepoint=tp,n_R=length(x),n_NR=length(y),R_median=median(x),NR_median=median(y),U_R=as.numeric(t$statistic),P=t$p.value)
}))));cross[,BH_q:=p.adjust(P,'BH')];cross[,BH_family:='R_vs_NR_timepoint_12']
fwrite(paired,file.path(o,'tables/paired_statistics.csv'));fwrite(delta,file.path(o,'tables/delta_R_vs_NR_statistics.csv'));fwrite(cross,file.path(o,'tables/timepoint_R_vs_NR_statistics.csv'))
# ID-bearing long table and paired data for each endpoint, with both immune values retained in master tables.
for(m in metrics){z<-wide[,c('patient_id','sample_PRE','sample_EDT','response','therapy','age','sex',paste0(c('PRE_','EDT_','Delta_'),m)),with=FALSE];fwrite(z,file.path(o,'tables',paste0('PRISM_',m,'_paired_ID.csv')))}
fmt<-function(x)ifelse(x<.001,formatC(x,format='e',digits=2),formatC(x,format='f',digits=3))
cols<-c(PRE='#2455D6',EDT='#E41A1C',R='#238B9A',NR='#CC7655')
theme_set(theme_classic(base_size=15)+theme(text=element_text(colour='#203844'),axis.text=element_text(colour='#203844'),plot.title=element_text(face='bold',size=17),plot.subtitle=element_text(size=12,colour='#546B78'),plot.caption=element_text(size=10,colour='#546B78',hjust=0),strip.background=element_blank(),strip.text=element_text(face='bold'),legend.position='none',plot.margin=margin(12,14,12,12)))
allplots<-list();stratplots<-list();delplots<-list();crossplots<-list()
for(m in metrics){
 z<-copy(d);z[,value:=get(m)];z[,timepoint:=factor(timepoint,levels=c('PRE','EDT'))];z[,response:=factor(response,levels=c('NR','R'))]
 pa<-paired[metric==m & group=='All'];ps<-paired[metric==m & group!='All'];st<-delta[metric==m];cs<-cross[metric==m]
 p1<-ggplot(z,aes(timepoint,value,group=patient_id))+geom_line(colour='#BBC4CA',linewidth=.5)+geom_point(aes(colour=timepoint),size=2.8,alpha=.9)+scale_colour_manual(values=cols)+labs(title=labels[m],subtitle=paste0('All matched patients | n = 16\nPaired BH q = ',fmt(pa$BH_q)),x=NULL,y=units[m])+scale_y_continuous(expand=expansion(mult=c(.07,.13)))
 ps[,response:=factor(group,levels=c('NR','R'))];ps[,lab:=paste0('n = ',n,' | BH q = ',fmt(BH_q))]
 p2<-ggplot(z,aes(timepoint,value,group=patient_id))+geom_line(colour='#BBC4CA',linewidth=.5)+geom_point(aes(colour=timepoint),size=2.6,alpha=.9)+scale_colour_manual(values=cols)+facet_wrap(~response,nrow=1)+geom_text(data=ps,aes(x=1.5,y=Inf,label=lab),inherit.aes=FALSE,vjust=1.3,size=3.6,colour='#546B78')+scale_y_continuous(expand=expansion(mult=c(.07,.24)))+labs(title=labels[m],subtitle='Paired PRE–EDT within each response group',x=NULL,y=units[m])
 zz<-wide[,.(patient_id,response,value=get(paste0('Delta_',m)))];zz[,response:=factor(response,levels=c('NR','R'))]
 p3<-ggplot(zz,aes(response,value,colour=response,fill=response))+geom_hline(yintercept=0,linetype='dashed',colour='#AAB4BD',linewidth=.5)+geom_boxplot(width=.48,outlier.shape=NA,alpha=.12,linewidth=.8)+geom_point(position=position_jitter(width=.10,height=0,seed=20260922),size=3)+scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+scale_x_discrete(labels=c(NR='NR\nn = 7',R='R\nn = 9'))+labs(title=labels[m],subtitle=paste0('EDT − PRE | R vs NR\nBH q = ',fmt(st$BH_q)),x=NULL,y=if(m=='CD8_percent') 'Change in CD8 estimate (percentage points)' else if(m=='CYT') 'Change in CYT (log2 scale)' else 'Change in expression (log2[FPKM + 1])')+scale_y_continuous(expand=expansion(mult=c(.10,.14)))
 cs[,timepoint:=factor(timepoint,levels=c('PRE','EDT'))];cs[,lab:=paste0('BH q = ',fmt(BH_q))]
 p4<-ggplot(z,aes(response,value,colour=response,fill=response))+geom_boxplot(width=.5,outlier.shape=NA,alpha=.12)+geom_point(position=position_jitter(width=.09,height=0,seed=20260922),size=2.7)+facet_wrap(~timepoint,nrow=1)+scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+geom_text(data=cs,aes(x=1.5,y=Inf,label=lab),inherit.aes=FALSE,vjust=1.3,size=3.6)+scale_y_continuous(expand=expansion(mult=c(.08,.22)))+scale_x_discrete(labels=c(NR='NR (7)',R='R (9)'))+labs(title=labels[m],subtitle='R vs NR at each timepoint',x=NULL,y=units[m])
 allplots[[m]]<-p1;stratplots[[m]]<-p2;delplots[[m]]<-p3;crossplots[[m]]<-p4
 combined<-(p1|p2|p3)+plot_layout(widths=c(1,1.4,1))+plot_annotation(title=paste0(labels[m],' | Matched melanoma patients'),caption='PRE: blue; EDT: red. Overall paired tests: BH across 6 endpoints; response-stratified paired tests: BH across 12 tests.\nChange comparison: two-sided Mann–Whitney; BH across 6 endpoints. Same 16 patients; all therapies retained.')
 for(ext in c('png','pdf'))ggsave(file.path(o,'figures',paste0(m,'_paired_and_delta.',ext)),combined,width=17,height=6,dpi=350,bg='white',device=if(ext=='png')ragg::agg_png else cairo_pdf,limitsize=FALSE)
}
views<-list(All_paired=allplots,By_response_paired=stratplots,Delta_R_vs_NR=delplots,R_vs_NR_at_each_timepoint=crossplots)
captions<-c(All_paired='Two-sided paired Wilcoxon signed-rank; BH across 6 endpoints. Same 16 patients; all therapies retained.',By_response_paired='Two-sided paired Wilcoxon signed-rank; BH across 12 tests (6 endpoints × 2 response groups).',Delta_R_vs_NR='Two-sided Mann–Whitney comparing patient-level EDT − PRE; BH across 6 endpoints. Dashed line: no change.',R_vs_NR_at_each_timepoint='Two-sided Mann–Whitney; BH across 12 tests (6 endpoints × 2 timepoints). These are cross-sectional comparisons.')
for(n in names(views)){
 fig<-wrap_plots(views[[n]],ncol=2)+plot_annotation(title='Melanoma | 16 matched PRE–EDT pairs',subtitle=paste0('Responders: 9 | Nonresponders: 7 | ',gsub('_',' ',n)),caption=captions[n])
 for(ext in c('png','pdf'))ggsave(file.path(o,'figures',paste0(n,'.',ext)),fig,width=12,height=15,dpi=300,bg='white',device=if(ext=='png')ragg::agg_png else cairo_pdf,limitsize=FALSE)
}
print(paired);print(delta);print(cross)
