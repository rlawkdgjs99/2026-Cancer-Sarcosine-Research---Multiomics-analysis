args <- commandArgs(trailingOnly=TRUE)
out <- normalizePath(args[1]); base <- normalizePath(file.path(out, '../..'))
library(ggplot2); library(patchwork)
options(digits=17)
src <- file.path(base,'results/Fig5_PRE73_Balance_CD8_CYT_Correlations_26.09.22/tables/patient_data_PRE73.csv')
osfile <- file.path(base,'results/Fig5_PRE73_Balance_HighLow_OS_26.09.22/tables/patient_data_PRE73.csv')
d <- read.csv(src,check.names=FALSE); os <- read.csv(osfile)
stopifnot(nrow(d)==73, !anyDuplicated(d$sample_id), !anyDuplicated(d$patient_id))
ix <- match(d$sample_id,os$sample_id)
stopifnot(!anyNA(ix),all(d$patient_id==os$patient_id[ix]),all(d$Balance==os$Balance[ix]))
cut <- median(d$Balance)
d$Balance_group <- factor(ifelse(d$Balance>cut,'High','Low'),levels=c('Low','High'))
stopifnot(all(as.character(d$Balance_group)==os$group[ix]),sum(d$Balance_group=='Low')==37,
          sum(d$Balance_group=='High')==36,all(is.finite(d$CD8_percent)),all(is.finite(d$CYT)),
          max(abs(d$CD8_percent-100*d$CD8_fraction))<1e-10,
          max(abs(d$CYT-log2(sqrt(d$GZMA_TPM*d$PRF1_TPM))))<1e-10)
write.csv(d,file.path(out,'tables/patient_data_PRE73.csv'),row.names=FALSE)
write.csv(data.frame(median=cut,Low_rule='<= median',High_rule='> median'),file.path(out,'tables/threshold.csv'),row.names=FALSE)
outcomes <- c('CD8_percent','CYT')
res <- do.call(rbind,lapply(outcomes,function(v){
 lo <- d[d$Balance_group=='Low',v]; hi <- d[d$Balance_group=='High',v]
 tt <- wilcox.test(hi,lo,alternative='two.sided',exact=FALSE,correct=TRUE)
 data.frame(outcome=v,n_Low=length(lo),n_High=length(hi),Low_Q1=quantile(lo,.25,names=FALSE),
 Low_median=median(lo),Low_Q3=quantile(lo,.75,names=FALSE),High_Q1=quantile(hi,.25,names=FALSE),
 High_median=median(hi),High_Q3=quantile(hi,.75,names=FALSE),U_High=unname(tt$statistic),
 rank_biserial=2*unname(tt$statistic)/(length(lo)*length(hi))-1,p=tt$p.value,method=tt$method)
}))
res$BH_q <- p.adjust(res$p,'BH');write.csv(res,file.path(out,'tables/comparison_statistics.csv'),row.names=FALSE)
plots <- list(); colors <- c(Low='#2455D6',High='#E41A1C')
for(v in outcomes){
 rr <- res[res$outcome==v,]; dd <- data.frame(sample_id=d$sample_id,group=d$Balance_group,value=d[[v]])
 set.seed(20260922);dd$x <- as.integer(dd$group)+runif(nrow(dd),-.15,.15)
 write.csv(dd,file.path(out,paste0('qa/plotted_',v,'.csv')),row.names=FALSE)
 lo <- dd$value[dd$group=='Low'];hi <- dd$value[dd$group=='High']
 write.csv(data.frame(Low=c(lo,rep(NA,37-length(lo))),High=c(hi,rep(NA,37-length(hi)))),
 file.path(out,paste0('tables/PRISM_',v,'.csv')),row.names=FALSE,na='')
 title <- if(v=='CD8_percent') 'CD8+ T-cell estimate' else 'Cytolytic expression (CYT)'
 yl <- if(v=='CD8_percent') 'quanTIseq CD8+ T cells (%)' else 'CYT score (log2 geometric mean)'
 subtitle <- sprintf('Melanoma | Pretreatment | n = 73\nBH q = %.3g',rr$BH_q)
 p <- ggplot(dd,aes(x=as.integer(group),y=value))+geom_boxplot(aes(group=group,fill=group,color=group),
 width=.48,linewidth=.75,outlier.shape=NA,alpha=.12)+
 geom_point(aes(x=x,color=group),size=2.8,alpha=.8)+
 scale_color_manual(values=colors)+scale_fill_manual(values=colors)+
 scale_x_continuous(breaks=1:2,labels=c('Balance Low\nn = 37','Balance High\nn = 36'),limits=c(.55,2.45))+
 scale_y_continuous(expand=expansion(mult=c(.045,.09)))+
 labs(title=title,subtitle=subtitle,x=NULL,y=yl)+
 theme_classic(base_size=16,base_family='Arial')+
 theme(legend.position='none',plot.title=element_text(size=22,face='bold',color='#213844',margin=margin(b=10)),
 plot.subtitle=element_text(size=13.5,lineheight=1.4,color='#4C6370',margin=margin(b=18)),
 axis.text=element_text(size=15,color='#213844'),axis.title.y=element_text(size=16,margin=margin(r=12)),
 axis.line=element_line(color='#617581',linewidth=.6),plot.margin=margin(18,18,12,15))
 plots[[v]] <- p
 single <- p+labs(caption='Points: individual patients. Box: median and IQR.\nTwo-sided Mann–Whitney; BH correction across two outcomes.')+
 theme(plot.caption=element_text(size=10.5,hjust=0,color='#617581',margin=margin(t=15)))
 for(ext in c('png','pdf')) ggsave(file.path(out,paste0('figures/Balance_HighLow_',v,'.',ext)),single,
 width=6.4,height=6.5,units='in',dpi=450,bg='white',device=if(ext=='png') ragg::agg_png else cairo_pdf)
}
combined <- (plots[[1]]+plots[[2]])+plot_annotation(caption=paste0(
 'Balance Low ≤ pooled median; High > pooled median. All 73 pretreatment patients retained.\n',
 'Points: individual patients; boxes: median and IQR; whiskers: 1.5 IQR. Two-sided Mann–Whitney; BH correction across two outcomes.'),
 theme=theme(plot.caption=element_text(size=11,hjust=0,color='#617581',margin=margin(t=12)),plot.margin=margin(5,12,8,12)))
for(ext in c('png','pdf')) ggsave(file.path(out,paste0('figures/Balance_HighLow_Combined.',ext)),combined,
 width=12.8,height=6.8,units='in',dpi=450,bg='white',device=if(ext=='png') ragg::agg_png else cairo_pdf)
writeLines(capture.output(sessionInfo()),file.path(out,'qa/sessionInfo.txt'))
print(res)
