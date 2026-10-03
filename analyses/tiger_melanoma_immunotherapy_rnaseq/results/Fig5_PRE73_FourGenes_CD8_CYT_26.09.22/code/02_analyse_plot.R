options(stringsAsFactors=FALSE,digits=17)
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
args<-commandArgs(FALSE);script<-normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]));out<-dirname(dirname(script))
d<-read.csv(file.path(out,'tables/patient_data_PRE73.csv'));stopifnot(nrow(d)==73,!anyDuplicated(d$sample_id),!anyDuplicated(d$patient_id))
metrics<-c('SARDH','PIPOX','GNMT','DMGDH');ends<-c('CD8_percent','CYT')
labels<-c(SARDH='SARDH expression',PIPOX='PIPOX expression',GNMT='GNMT expression',DMGDH='DMGDH expression')
prefix<-c(SARDH='GS',PIPOX='GP',GNMT='GG',DMGDH='GD')
d$age_rank<-rank(d$age,ties.method='average');d$sex<-factor(d$sex,levels=c('Female','Male'));d$therapy<-factor(d$therapy,levels=c('antiPD1','combo'))
X<-model.matrix(~age_rank+sex+therapy,d);stopifnot(qr(X)$rank==4)
write.csv(X,file.path(out,'qa/design_matrix.csv'),row.names=FALSE)
cors<-list();groups<-list();resids<-list()
for(m in metrics)for(v in ends){
 stopifnot(all(is.finite(d[[m]])),all(is.finite(d[[v]])))
 for(model in c('Adjusted','Unadjusted')){
  rx<-rank(d[[m]],ties.method='average');ry<-rank(d[[v]],ties.method='average')
  Z<-if(model=='Adjusted')X else matrix(1,nrow(d),1)
  ex<-qr.resid(qr(Z),rx);ey<-qr.resid(qr(Z),ry);rho<-cor(ex,ey);k<-ncol(Z)-1;df<-73-k-2
  t<-rho*sqrt(df/(1-rho^2));P<-2*pt(-abs(t),df);ci<-tanh(atanh(rho)+c(-1,1)*qnorm(.975)/sqrt(73-k-3))
  olsP<-coef(summary(lm(ry~0+Z+rx)))['rx','Pr(>|t|)'];stopifnot(abs(log(P)-log(olsP))<1e-10)
  cors[[length(cors)+1]]<-data.frame(metric=m,endpoint=v,model=model,n=73,rho=rho,P=P,covariate_df=k,df=df,t=t,Fisher_lower95=ci[1],Fisher_upper95=ci[2])
  resids[[length(resids)+1]]<-data.frame(sample_id=d$sample_id,metric=m,endpoint=v,model=model,x_residual=ex,y_residual=ey)
 }
 gr<-d[[paste0(m,'_group')]];lo<-d[gr=='Low',v];hi<-d[gr=='High',v]
 stopifnot(all(gr==ifelse(d[[m]]>median(d[[m]]),'High','Low')))
 tt<-wilcox.test(hi,lo,alternative='two.sided',exact=FALSE,correct=TRUE)
 groups[[length(groups)+1]]<-data.frame(metric=m,endpoint=v,n_Low=length(lo),n_High=length(hi),Low_Q1=quantile(lo,.25,names=FALSE),Low_median=median(lo),Low_Q3=quantile(lo,.75,names=FALSE),High_Q1=quantile(hi,.25,names=FALSE),High_median=median(hi),High_Q3=quantile(hi,.75,names=FALSE),U_High=unname(tt$statistic),rank_biserial=2*unname(tt$statistic)/(length(lo)*length(hi))-1,P=tt$p.value)
}
cr<-do.call(rbind,cors);cr$BH_q<-NA_real_
for(model in unique(cr$model)){ii<-cr$model==model;cr$BH_q[ii]<-p.adjust(cr$P[ii],'BH')}
gr<-do.call(rbind,groups);gr$BH_q<-p.adjust(gr$P,'BH')
write.csv(cr,file.path(out,'tables/correlation_statistics.csv'),row.names=FALSE)
write.csv(gr,file.path(out,'tables/group_statistics.csv'),row.names=FALSE)
write.csv(do.call(rbind,resids),file.path(out,'tables/rank_residuals.csv'),row.names=FALSE)
fmt<-function(x)if(x<.001)formatC(x,format='e',digits=2)else formatC(x,format='g',digits=3)
INK<-'#223843';MUTED<-'#647783';COL<-c(CD8_percent='#258E94',CYT='#C25458');GC<-c(Low='#2455D6',High='#E41A1C')
base_theme<-theme_classic(base_size=15,base_family='Arial')+theme(text=element_text(colour=INK),legend.position='none',plot.title=element_text(size=21,face='bold',margin=margin(b=7)),plot.subtitle=element_text(size=12.5,lineheight=1.5,colour=MUTED,margin=margin(b=15)),axis.title=element_text(size=15),axis.title.x=element_text(margin=margin(t=10)),axis.title.y=element_text(margin=margin(r=10)),axis.text=element_text(size=13,colour=INK),axis.line=element_line(colour=MUTED,linewidth=.5),axis.ticks=element_line(colour=MUTED),plot.margin=margin(12,15,10,12))
save_plot<-function(p,stem,w,h){
 ggsave(file.path(out,'figures',paste0(stem,'.png')),p,width=w,height=h,dpi=450,bg='white',device=ragg::agg_png)
 ggsave(file.path(out,'figures',paste0(stem,'.pdf')),p,width=w,height=h,bg='white',device=cairo_pdf)
}
all_sc<-list();all_gp<-list();exports<-list()
for(m in metrics){
 sc<-list();gp<-list()
 for(v in ends){
  title<-if(v=='CD8_percent')'CD8+ T-cell estimate' else 'Cytolytic expression (CYT)'
  yl<-if(v=='CD8_percent')'quanTIseq CD8+ T cells (%)' else 'CYT score (log2 geometric mean)'
  rr<-cr[cr$metric==m&cr$endpoint==v&cr$model=='Adjusted',];gg<-gr[gr$metric==m&gr$endpoint==v,]
  s<-ggplot(d,aes(x=.data[[m]],y=.data[[v]]))+
   geom_smooth(method='lm',formula=y~x,se=FALSE,colour=COL[v],linetype='dashed',linewidth=.75)+
   geom_point(shape=21,fill=COL[v],colour='white',stroke=.35,size=3.2,alpha=.88)+
   scale_x_continuous(labels=function(x) formatC(x,format='f',digits=1),expand=expansion(mult=c(.06,.06)))+scale_y_continuous(expand=expansion(mult=c(.06,.1)))+
   labs(title=title,subtitle=paste0('Melanoma | Pretreatment | n = 73\n',sprintf('Adjusted rho = %+.2f   |   BH q = %s',rr$rho,fmt(rr$BH_q))),x=paste0(labels[m],'\nlog2(FPKM + 1)'),y=yl)+base_theme
  pd<-ggplot_build(s)$data[[2]];stopifnot(nrow(pd)==73,max(abs(pd$x-d[[m]]))<1e-12,max(abs(pd$y-d[[v]]))<1e-12)
  write.csv(data.frame(sample_id=d$sample_id,x=pd$x,y=pd$y),file.path(out,'qa',paste0(m,'_',v,'_scatter.csv')),row.names=FALSE)
  write.csv(d[,c('sample_id',m,v)],file.path(out,'tables',paste0('PRISM_',m,'_',v,'_XY.csv')),row.names=FALSE)
  sc[[v]]<-s;all_sc[[paste(m,v)]]<-s
  sc_cap<-'Partial Spearman: adjusted for age, sex and therapy; BH across eight comparisons.\nPoints: observed values. Dashed line: unadjusted linear trend.'
  stem<-paste0(m,'_',v,'_scatter');save_plot(s+labs(caption=sc_cap)+theme(plot.caption=element_text(size=9,hjust=0,colour=MUTED,margin=margin(t=14))),stem,6.5,5.9)
  short<-paste0(prefix[m],'S',if(v=='CD8_percent')'C' else 'Y','.png');exports[[length(exports)+1]]<-data.frame(stem=stem,short=short)
  dd<-data.frame(sample_id=d$sample_id,group=factor(d[[paste0(m,'_group')]],levels=c('Low','High')),value=d[[v]])
  set.seed(20260922);dd$x<-as.integer(dd$group)+runif(73,-.15,.15)
  write.csv(dd,file.path(out,'qa',paste0(m,'_',v,'_group.csv')),row.names=FALSE)
  p<-ggplot(dd,aes(x=as.integer(group),y=value))+
   geom_boxplot(aes(group=group,fill=group,colour=group),width=.48,linewidth=.75,outlier.shape=NA,alpha=.12)+
   geom_point(aes(x=x,colour=group),size=2.8,alpha=.8)+scale_color_manual(values=GC)+scale_fill_manual(values=GC)+
   scale_x_continuous(breaks=1:2,labels=c(paste0('Low\nn = ',gg$n_Low),paste0('High\nn = ',gg$n_High)),limits=c(.55,2.45))+
   scale_y_continuous(expand=expansion(mult=c(.045,.09)))+
   labs(title=title,subtitle=paste0('Melanoma | Pretreatment | n = 73\nBH q = ',fmt(gg$BH_q)),x=labels[m],y=yl)+base_theme
  gp[[v]]<-p;all_gp[[paste(m,v)]]<-p
  gp_cap<-'Points: patients. Box: median/IQR; whiskers: 1.5 IQR.\nTwo-sided Mann–Whitney; BH across eight comparisons.'
  stem<-paste0(m,'_',v,'_groups');save_plot(p+labs(caption=gp_cap)+theme(plot.caption=element_text(size=9,hjust=0,colour=MUTED,margin=margin(t=14))),stem,6.5,6.3)
  short<-paste0(prefix[m],'G',if(v=='CD8_percent')'C' else 'Y','.png');exports[[length(exports)+1]]<-data.frame(stem=stem,short=short)
  lo<-dd$value[dd$group=='Low'];hi<-dd$value[dd$group=='High'];nn<-max(length(lo),length(hi))
  write.csv(data.frame(Low=c(lo,rep(NA,nn-length(lo))),High=c(hi,rep(NA,nn-length(hi)))),file.path(out,'tables',paste0('PRISM_',m,'_',v,'_groups.csv')),row.names=FALSE,na='')
 }
 annotation_theme<-theme(plot.caption=element_text(size=10.5,hjust=0,colour=MUTED),plot.margin=margin(5,10,8,10))
 for(kind in c('scatter','groups')){
  pp<-if(kind=='scatter')sc else gp
  cap<-if(kind=='scatter')sc_cap else 'Low ≤ pooled median; High > pooled median. All 73 patients retained. Points: patients; boxes: median/IQR; whiskers: 1.5 IQR.\nTwo-sided Mann–Whitney (no clinical covariate adjustment); BH correction across eight comparisons.'
  stem<-paste0(m,'_',kind,'_combined');save_plot(wrap_plots(pp,ncol=2)+plot_annotation(caption=cap,theme=annotation_theme),stem,13,6.3)
  exports[[length(exports)+1]]<-data.frame(stem=stem,short=paste0(prefix[m],if(kind=='scatter')'S' else 'G','.png'))
 }
}
for(kind in c('scatter','groups')){
 pp<-if(kind=='scatter')all_sc else all_gp
 stem<-paste0('All_metrics_',kind)
 cap<-if(kind=='scatter')'Adjusted partial Spearman (age, sex, therapy). BH across eight comparisons. Observed points; dashed lines: unadjusted linear trends.' else 'Pooled median splits; Low ≤ median, High > median. Two-sided Mann–Whitney; BH across eight comparisons. Boxes: median/IQR.'
 save_plot(wrap_plots(pp,ncol=2)+plot_annotation(caption=cap,theme=theme(plot.caption=element_text(size=11,hjust=0,colour=MUTED))),stem,13,22.1)
 exports[[length(exports)+1]]<-data.frame(stem=stem,short=if(kind=='scatter')'G4S.png' else 'G4G.png')
}
write.csv(do.call(rbind,exports),file.path(out,'qa/figure_exports.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'qa/sessionInfo.txt'))
print(cr[cr$model=='Adjusted',]);print(gr)
