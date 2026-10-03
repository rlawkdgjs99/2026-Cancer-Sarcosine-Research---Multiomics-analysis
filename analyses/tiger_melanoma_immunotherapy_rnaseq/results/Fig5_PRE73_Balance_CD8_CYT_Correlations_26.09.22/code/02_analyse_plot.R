options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
args<-commandArgs(FALSE);script<-normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]));out<-dirname(dirname(script))
d<-read.csv(file.path(out,'tables/patient_data_PRE73.csv'))
stopifnot(nrow(d)==73,!anyDuplicated(d$sample_id),all(complete.cases(d)))
d$age_rank<-rank(d$age,ties.method='average');d$sex<-factor(d$sex,levels=c('Female','Male'));d$therapy<-factor(d$therapy,levels=c('antiPD1','combo'))
X<-model.matrix(~age_rank+sex+therapy,d);stopifnot(ncol(X)==4,qr(X)$rank==4)
write.csv(X,file.path(out,'qa/design_matrix.csv'),row.names=FALSE)
results<-list();residuals<-list()
for(v in c('CD8_percent','CYT'))for(model in c('Adjusted','Unadjusted')){
 rx<-rank(d$Balance,ties.method='average');ry<-rank(d[[v]],ties.method='average')
 Z<-if(model=='Adjusted')X else matrix(1,nrow(d),1)
 ex<-qr.resid(qr(Z),rx);ey<-qr.resid(qr(Z),ry);rho<-cor(ex,ey);k<-ncol(Z)-1;df<-nrow(d)-k-2
 t<-rho*sqrt(df/(1-rho^2));P<-2*pt(-abs(t),df);ci<-tanh(atanh(rho)+c(-1,1)*qnorm(.975)/sqrt(nrow(d)-k-3))
 olsP<-coef(summary(lm(ry~0+Z+rx)))['rx','Pr(>|t|)'];stopifnot(abs(log(P)-log(olsP))<1e-10)
 if(model=='Adjusted'){
  prec<-solve(cor(cbind(rx,ry,X[,-1])));stopifnot(abs(rho+prec[1,2]/sqrt(prec[1,1]*prec[2,2]))<1e-12)
  deltares<-qr.resid(qr(Z),rank(d$Delta_sampleSD));stopifnot(abs(rho-cor(deltares,ey))<1e-12)
 }
 results[[length(results)+1]]<-data.frame(endpoint=v,model=model,n=nrow(d),rho=rho,P=P,covariate_df=k,df=df,t=t,Fisher_lower95=ci[1],Fisher_upper95=ci[2])
 residuals[[length(residuals)+1]]<-data.frame(sample_id=d$sample_id,endpoint=v,model=model,Balance_rank=rx,outcome_rank=ry,Balance_residual=ex,outcome_residual=ey)
}
r<-do.call(rbind,results);r$BH_q<-NA_real_
for(m in unique(r$model)){ii<-r$model==m;r$BH_q[ii]<-p.adjust(r$P[ii],'BH')}
write.csv(r,file.path(out,'tables/correlation_statistics.csv'),row.names=FALSE)
write.csv(do.call(rbind,residuals),file.path(out,'tables/rank_residuals.csv'),row.names=FALSE)
fmt<-function(x)if(x<.001)formatC(x,format='e',digits=2)else sprintf('%.3f',x)
INK<-'#223843';MUTED<-'#647783';COL<-c(CD8_percent='#258E94',CYT='#C25458')
plots<-list()
for(v in c('CD8_percent','CYT')){
 a<-r[r$endpoint==v&r$model=='Adjusted',];title<-if(v=='CD8_percent')'CD8+ T-cell estimate' else 'Cytolytic expression (CYT)'
 yl<-if(v=='CD8_percent')'quanTIseq CD8+ T cells (%)' else 'CYT score (log2 geometric mean)'
 lab<-sprintf('Adjusted rho = %+.2f   |   BH q = %s',a$rho,fmt(a$BH_q))
 p<-ggplot(d,aes(x=Balance,y=.data[[v]]))+
 geom_smooth(method='lm',formula=y~x,se=FALSE,colour=COL[v],linetype='dashed',linewidth=.75)+
 geom_point(shape=21,fill=COL[v],colour='white',stroke=.35,size=3.2,alpha=.88)+
 scale_x_continuous(breaks=seq(-3,3,1),expand=expansion(mult=c(.05,.06)))+
 scale_y_continuous(expand=expansion(mult=c(.06,.1)))+
 labs(title=title,subtitle=paste0('Melanoma | Pretreatment | n = 73\n',lab),x='Sarcosine metabolic balance\nDegradation - Production',y=yl)+
 theme_classic(base_size=15,base_family='Arial')+
 theme(text=element_text(colour=INK),plot.title=element_text(size=21,face='bold',margin=margin(b=7)),plot.subtitle=element_text(size=12.5,lineheight=1.5,colour=MUTED,margin=margin(b=15)),axis.title=element_text(size=15),axis.title.x=element_text(margin=margin(t=10),lineheight=1.2),axis.title.y=element_text(margin=margin(r=10)),axis.text=element_text(size=13,colour=INK),axis.line=element_line(colour=MUTED,linewidth=.5),axis.ticks=element_line(colour=MUTED),plot.margin=margin(12,15,10,12))
 if(v=='CD8_percent')p<-p+expand_limits(y=0)
 plots[[v]]<-p
 stem<-if(v=='CD8_percent')'Balance_CD8' else 'Balance_CYT'
 cap<-'Partial Spearman: adjusted for age, sex and therapy; BH across two outcomes.\nPoints show observed values; dashed line is an unadjusted linear trend.'
 p<-p+labs(caption=cap)+theme(plot.caption=element_text(size=9,hjust=0,lineheight=1.2,colour=MUTED,margin=margin(t=14)))
 ggsave(file.path(out,'figures',paste0(stem,'.png')),p,width=6.5,height=5.8,dpi=450,bg='white',device=ragg::agg_png)
 ggsave(file.path(out,'figures',paste0(stem,'.pdf')),p,width=6.5,height=5.8,bg='white',device=cairo_pdf)
 # Assert all observed points survived the plotting layer without transformation or filtering.
 pd<-ggplot_build(p)$data[[2]];stopifnot(nrow(pd)==73,max(abs(pd$x-d$Balance))<1e-12,max(abs(pd$y-d[[v]]))<1e-12)
 write.csv(data.frame(sample_id=d$sample_id,x=pd$x,y=pd$y),file.path(out,'qa',paste0(stem,'_plotted.csv')),row.names=FALSE)
}
combo<-wrap_plots(plots,ncol=2)+plot_annotation(caption='Partial Spearman correlations adjusted for age, sex and therapy; BH correction across two outcomes.\nAll 73 patients retained. Points show observed values; dashed lines are unadjusted linear trends.',theme=theme(plot.caption=element_text(size=10.5,hjust=0,lineheight=1.25,colour=MUTED),plot.margin=margin(5,10,8,10)))
ggsave(file.path(out,'figures/Balance_Combined.png'),combo,width=13,height=5.8,dpi=450,bg='white',device=ragg::agg_png)
ggsave(file.path(out,'figures/Balance_Combined.pdf'),combo,width=13,height=5.8,bg='white',device=cairo_pdf)
capture.output(sessionInfo(),file=file.path(out,'qa/sessionInfo.txt'));print(r)
