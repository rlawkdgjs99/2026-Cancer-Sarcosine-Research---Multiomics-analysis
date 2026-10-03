#!/usr/bin/env Rscript
# Targeted, patient-level relative-proportion comparisons; frozen plan is authoritative.
options(device=function(...)grDevices::cairo_pdf(file=tempfile(fileext='.pdf'),family='Arial',...),stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(gridExtra);library(grid);library(jsonlite)})
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
base <- normalizePath(file.path(dirname(script),'..')); tabs<-file.path(base,'results/tables'); figs<-file.path(base,'results/figures')
plan<-fromJSON(file.path(tabs,'Fig8_CellFrac_plan.json'))
md<-fread(plan$inputs$patients$path); comp<-fread(plan$inputs$composition$path)
post<-md[treatment=='Post'][order(Patient)]
stopifnot(nrow(post)==12L,!anyDuplicated(post$Patient),sum(post$group=='High')==7L,sum(post$group=='Low')==5L)
stopifnot(all(ifelse(md$degradation_mean_z>md$cutoff,'High','Low')==md$group))
post[,histology:=factor(histology)];post[,High:=as.integer(group=='High')]
stopifnot(identical(levels(post$histology),c('Adeno','Squamous')))
counts<-comp[Patient %in% post$Patient];stopifnot(!anyDuplicated(counts[,.(Patient,final_lineage)]))
ct<-counts[,.(total=sum(cells)),by=Patient]
stopifnot(all(ct$total[match(post$Patient,ct$Patient)]==post$all_cells))
# Enumerate unique assignments, preserving each histology's observed High total.
ids<-split(seq_len(nrow(post)),post$histology)
choices<-lapply(ids,function(ix)combn(ix,sum(post$High[ix]),simplify=FALSE))
perms<-lapply(choices[[1]],function(a)lapply(choices[[2]],function(z){h<-integer(nrow(post));h[c(a,z)]<-1L;h}))
perms<-do.call(rbind,unlist(perms,recursive=FALSE));stopifnot(nrow(perms)==300L,!anyDuplicated(as.data.frame(perms)),sum(apply(perms,1,function(x)all(x==post$High)))==1L)
model<-function(y,h){m<-lm(y~histology+h,data=data.frame(y=y,histology=post$histology,h=h));c(beta=coef(m)['h'],t=coef(summary(m))['h','t value'])}
rows<-list(); stats<-list(); perm_rows<-list();summ<-list();kk<-0L
for(scope in c('all_singlets','immune_singlets')){
 den<-if(scope=='all_singlets')post$all_cells else counts[!final_lineage %in% c('Epithelial','CAF'),.(den=sum(cells)),by=Patient][match(post$Patient,Patient),den]
 stopifnot(all(is.finite(den)),all(den>0))
 for(cell in c('CD8 T cell','Conventional DC')){
  kk<-kk+1L
  num<-counts[final_lineage==cell][match(post$Patient,Patient),cells]
  stopifnot(length(num)==12L,all(num>0),all(num<den))
  pct<-100*num/den;y<-log(num/(den-num));ob<-model(y,post$High)
  ts<-apply(perms,1,function(h)model(y,h)[2]);n_extreme<-sum(abs(ts)>=abs(ob[2])-1e-12);pv<-n_extreme/length(ts)
  d<-data.table(Patient=post$Patient,Sample=post$Sample,histology=as.character(post$histology),group=post$group,scope=scope,cell_type=cell,target_cells=num,denominator_cells=den,percent=pct,log_ratio=y)
  rows[[kk]]<-d
  m<-d[,.(n=.N,mean_percent=mean(percent),SD_percent=sd(percent),median_percent=median(percent)),by=.(scope,cell_type,group)]
  summ[[kk]]<-m
  delta<-m[group=='High',mean_percent]-m[group=='Low',mean_percent]
  stats[[kk]]<-data.table(scope=scope,cell_type=cell,n_Low=5L,n_High=7L,log_ratio_beta_High_minus_Low=unname(ob[1]),t_statistic=unname(ob[2]),raw_mean_difference_pp=delta,extreme_assignments=n_extreme,total_assignments=length(ts),raw_P=pv)
  perm_rows[[kk]]<-data.table(scope=scope,cell_type=cell,assignment=seq_along(ts),t_statistic=ts,is_observed=apply(perms,1,function(x)all(x==post$High)))
 }
}
plotdata<-rbindlist(rows);res<-rbindlist(stats);res[,BH_q:=p.adjust(raw_P,'BH'),by=scope]
summary_table<-rbindlist(summ)
fwrite(plotdata,file.path(tabs,'Fig8_CellFrac_patients.csv'));fwrite(res,file.path(tabs,'Fig8_CellFrac_statistics.csv'));fwrite(summary_table,file.path(tabs,'Fig8_CellFrac_summary.csv'));fwrite(rbindlist(perm_rows),file.path(tabs,'Fig8_CellFrac_permutations.csv'))
fwrite(data.table(Patient=post$Patient,histology=as.character(post$histology),group=post$group),file.path(tabs,'Fig8_CellFrac_roster.csv'))
# Minimal Prism plot input: one percentage column per group, patient IDs retained.
for(cell in c('CD8 T cell','Conventional DC')){
 d<-plotdata[scope=='all_singlets' & cell_type==cell]
 lo<-d[group=='Low']; hi<-d[group=='High'];n<-max(nrow(lo),nrow(hi))
 pad<-function(x)c(x,rep(NA,n-length(x)))
 z<-data.table(Low_Patient=pad(lo$Patient),Low_percent=pad(lo$percent),High_Patient=pad(hi$Patient),High_percent=pad(hi$percent))
 fwrite(z,file.path(tabs,paste0('Fig8_',if(cell=='CD8 T cell')'CD8' else 'cDC','_Prism.csv')),na='')
}
# Fixed x offsets prevent random jitter changes between exports.
cols<-c(Low='#2E8798',High='#C3466C');fills<-c(Low='#D7EAED',High='#F2D8E0')
make_plot<-function(cell){
 d<-copy(plotdata[scope=='all_singlets' & cell_type==cell]);m<-copy(summary_table[scope=='all_singlets' & cell_type==cell]);q<-res[scope=='all_singlets' & cell_type==cell,BH_q]
 d[,group:=factor(group,levels=c('Low','High'))];m[,group:=factor(group,levels=c('Low','High'))]
 d[,x:=as.integer(group)+seq(-.14,.14,length.out=.N),by=group];m[,x:=as.integer(group)]
 yy<-max(c(d$percent,m$mean_percent+m$SD_percent));top<-if(cell=='CD8 T cell')ceiling(yy*1.22/10)*10 else ceiling(yy*1.22*2)/2
 bracket<-top*.88
 g<-ggplot()+
 geom_col(data=m,aes(x=x,y=mean_percent,fill=group,colour=group),width=.58,linewidth=.7)+
 geom_errorbar(data=m,aes(x=x,ymin=mean_percent-SD_percent,ymax=mean_percent+SD_percent,colour=group),width=.16,linewidth=.65)+
 geom_point(data=d,aes(x=x,y=percent,fill=group),shape=21,size=2.55,stroke=.4,colour='white')+
 geom_segment(aes(x=1,xend=2,y=bracket,yend=bracket),linewidth=.4,colour='#526470')+
 geom_segment(aes(x=1,xend=1,y=bracket,yend=bracket-top*.022),linewidth=.4,colour='#526470')+
 geom_segment(aes(x=2,xend=2,y=bracket,yend=bracket-top*.022),linewidth=.4,colour='#526470')+
 annotate('text',x=1.5,y=bracket+top*.06,label=sprintf('BH q = %.3f',q),family='Arial',size=3.1,colour='#20343F')+
 scale_x_continuous(breaks=1:2,labels=c('Low\nn = 5','High\nn = 7'),limits=c(.48,2.52),expand=c(0,0))+
 scale_y_continuous(limits=c(0,top),expand=expansion(mult=c(0,0)),breaks=pretty(c(0,top),n=5))+
 scale_colour_manual(values=cols)+scale_fill_manual(values=fills)+
 labs(title=if(cell=='CD8 T cell')'CD8⁺ T cells' else 'Conventional DC',x='Degradation group',y='Fraction of all captured cells (%)')+
 theme_classic(base_size=11,base_family='Arial')+
 theme(legend.position='none',plot.title=element_text(size=13,face='bold',colour='#20343F',margin=margin(b=12)),axis.title=element_text(size=10,colour='#20343F'),axis.text=element_text(size=9.5,colour='#304550'),axis.line=element_line(colour='#526470',linewidth=.45),axis.ticks=element_line(colour='#526470',linewidth=.45),axis.ticks.length=unit(1.5,'mm'),plot.margin=margin(7,13,5,7))
 # Use darker fills for the patient point layer; bars retain pale fills.
 g$layers[[3]]$aes_params$fill<-unname(cols[as.character(d$group)])
 g
}
p1<-make_plot('CD8 T cell');p2<-make_plot('Conventional DC')
# A null graphics fallback prevents overwriting project-root Rplots.pdf.
header<-grobTree(textGrob('CD8⁺ T-cell and DC proportions',x=.025,y=.76,just='left',gp=gpar(fontfamily='Arial',fontsize=15,fontface='bold',col='#20343F')),textGrob('Post-treatment NSCLC · fixed whole-tumour degradation groups',x=.025,y=.26,just='left',gp=gpar(fontfamily='Arial',fontsize=9,col='#617582')))
foot<-textGrob('Bars: mean ± SD. Each point represents one patient.\nHistology-stratified permutation tests on log-ratios; BH across the two cell types.',x=.025,y=.68,just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=8,col='#617582',lineheight=1.25))
plots<-list(CD=list(width=6.8,height=4.15,content=arrangeGrob(header,arrangeGrob(p1,p2,ncol=2),foot,ncol=1,heights=c(.55,3.08,.52))),C8=list(width=3.45,height=3.5,content=ggplotGrob(p1)),DC=list(width=3.45,height=3.5,content=ggplotGrob(p2)))
for(stem in names(plots)){
 z<-plots[[stem]]
 grDevices::cairo_pdf(file.path(figs,paste0(stem,'.pdf')),width=z$width,height=z$height,family='Arial',bg='white');grid.newpage();grid.draw(z$content);dev.off()
 pngname<-if(stem=='CD')'Fig8_CD8_cDC_Proportions.png' else paste0(stem,'.png')
 ragg::agg_png(file.path(figs,pngname),width=z$width,height=z$height,units='in',res=600,background='white');grid.newpage();grid.draw(z$content);dev.off()
}
sink(file.path(tabs,'Fig8_CellFrac_sessionInfo.txt'));print(sessionInfo());sink()
print(res)
print(summary_table)
cat('Completed fixed patient-level tests and bar plots. No prior outputs overwritten.\n')

if(!is.null(warnings()))print(warnings())
