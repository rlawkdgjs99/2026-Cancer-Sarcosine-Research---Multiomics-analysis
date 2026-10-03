#!/usr/bin/env Rscript
# AUTHOR requested the actual CRC seven-KO workflow, not the old NSCLC
# prevalence-gated MaAsLin2 meta-estimate display. Separate exploratory branch.
# Reference: CRC scripts 121 (seven-KO inclusion/BH) and 126 (pooled r_rb/CI).
# Unit: original selected Run record. No cohort/covariate adjustment or dedup.
# Known cross-project Sample-name overlap is retained, NOT claimed independent.
suppressPackageStartupMessages({library(ggplot2);library(ggtext);library(ggrepel);
  library(patchwork);library(ragg);library(digest)})
options(stringsAsFactors=FALSE)
arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
script <- normalizePath(sub('^--file=','',arg))
pool <- dirname(dirname(script)); nsroot <- dirname(pool); ar <- dirname(nsroot)
crcroot <- file.path(ar,'HGMT_CRC_WGS-Healthy_vs_Cancer')
crci <- file.path(crcroot,'Healthy_vs_Cancer_4_CRC_cohorts_integrated')
out <- file.path(pool,'results','CRC_Method_7KO_NSCLC_26.09.07')
dir.create(out,showWarnings=FALSE)
for(d in c('figures','tables','qa')) dir.create(file.path(out,d),showWarnings=FALSE)
target <- c('K00301','K00302','K00303','K00305','K00306','K00315','K08688')
gene <- c('Sarcosine oxidase','soxA','soxB','soxG','PIPOX','DMGDH','Creatinase')
role <- c(rep('Degradation',5),rep('Production',2))
cohorts <- c('PRJNA751792','PRJNA1023797','PRJEB22863')
used <- file.path(nsroot,'사용데이터_모음')
# macOS NFD pathname, resolve by an exact normalized match if needed.
if(!dir.exists(used)) used <- list.dirs(nsroot,recursive=FALSE,full.names=TRUE)[
  stringi::stri_trans_nfc(basename(list.dirs(nsroot,recursive=FALSE,full.names=TRUE)))=='사용데이터_모음']
crcused <- file.path(crcroot,'사용데이터_모음','CRC_WGS_4cohort_pooled_26.09.01')
if(!dir.exists(crcused)) {
  parent <- list.dirs(crcroot,recursive=FALSE,full.names=TRUE)
  crcused <- file.path(parent[stringi::stri_trans_nfc(basename(parent))=='사용데이터_모음'],'CRC_WGS_4cohort_pooled_26.09.01')
}
nsfiles <- file.path(used,cohorts,paste0(cohorts,'_KEGG_KO_relative_abundance_matrix.csv'))
mdfiles <- file.path(used,cohorts,paste0(cohorts,'_patient_metadata_WGS_NSCLC_R_NR.csv'))
crcfile <- file.path(crcused,'CRC_WGS_4cohort_pooled_KEGG_KO_relative_abundance_matrix.csv')
crctarget <- file.path(crci,'results_integrated','sarcosine','sarcosine_KO_per_sample_pooled.csv')
crcfrozen <- file.path(crci,'results_integrated','sarcosine','sarcosine_KO_comparison_filtered_pooled.csv')
inputs <- c(nsfiles,mdfiles,crcfile,crctarget,crcfrozen,
  file.path(crci,'scripts',c('121_CRC_WGS_sarcosine_7KO_lollipop_26.09.01.R','126_remaining_5KO_pooled_forests_26.09.01.R','integrated_analysis_pooled.R')),
  file.path(used,'SOURCE_MANIFEST.csv'),file.path(used,'POST_EXPORT_CELL_LEVEL_VALIDATION.csv'))
stopifnot(all(file.exists(inputs)))
hash <- function(p) digest(file=p,algo='sha256',serialize=FALSE)
hashes <- vapply(inputs,hash,'')
write.csv(data.frame(path=inputs,SHA256=hashes),file.path(out,'qa','input_sha256.csv'),row.names=FALSE)
read <- function(p) read.csv(p,check.names=FALSE,fileEncoding='UTF-8-BOM')
write <- function(d,n) write.csv(d,file.path(out,'tables',n),row.names=FALSE,na='')
message('Reading the three NSCLC analysis-used KO matrices...')
xs <- lapply(nsfiles,read); mds <- lapply(mdfiles,read)
for(i in seq_along(xs)) {
  x <- xs[[i]]; m <- mds[[i]]
  stopifnot(!anyDuplicated(x$Run_ID),!anyDuplicated(m[['Run ID']]),
    setequal(x$Run_ID,m[['Run ID']]))
  m <- m[match(x$Run_ID,m[['Run ID']]),,drop=FALSE]
  stopifnot(all(m$Cohort==cohorts[i]),all(m[['Analysis Group']] %in% c('R','NR')))
  mds[[i]] <- m
  values <- as.matrix(x[-1]);storage.mode(values)<-'double'
  stopifnot(all(is.finite(values)),all(values>=0),max(abs(rowSums(values)-100))<.001)
}
allko <- sort(unique(unlist(lapply(xs,function(x) names(x)[-1]))))
ns <- matrix(0,sum(vapply(xs,nrow,integer(1))),length(allko),dimnames=list(NULL,allko))
offset <- 0L
for(x in xs) {rows <- seq_len(nrow(x))+offset; ns[rows,match(names(x)[-1],allko)]<-as.matrix(x[-1]);offset<-max(rows)}
meta <- do.call(rbind,lapply(seq_along(xs),function(i) data.frame(Run_ID=xs[[i]]$Run_ID,Cohort=cohorts[i],
  Group=mds[[i]][['Analysis Group']],Sample_name=mds[[i]][['Sample name']])))
stopifnot(nrow(ns)==824,!anyDuplicated(meta$Run_ID),sum(meta$Group=='R')==432,sum(meta$Group=='NR')==392,
  identical(as.integer(table(factor(meta$Cohort,levels=cohorts))),c(338L,421L,65L)),all(target %in% colnames(ns)))
rm(xs);gc()
nr <- meta$Group=='R'; nnr <- !nr
write(cbind(meta,as.data.frame(ns[,target,drop=FALSE],check.names=FALSE)),'NSCLC_seven_KO_per_Run.csv')
det <- do.call(rbind,lapply(cohorts,function(c) do.call(rbind,lapply(target,function(k)
  data.frame(Cohort=c,KO=k,n=sum(meta$Cohort==c),positive=sum(ns[meta$Cohort==c,k]>0),
             positive_R=sum(ns[meta$Cohort==c & nr,k]>0),positive_NR=sum(ns[meta$Cohort==c & nnr,k]>0))))))
write(det,'cohort_detection_counts.csv')
# Exactly the CRC five-KO degradation/two-KO production definitions and ratio.
scores <- cbind(meta,Degradation=rowSums(ns[,target[1:5],drop=FALSE]),Production=rowSums(ns[,target[6:7],drop=FALSE]))
scores$Production_Degradation <- log2((scores$Production+1e-8)/(scores$Degradation+1e-8))
stopifnot(all(is.finite(as.matrix(scores[c('Degradation','Production','Production_Degradation')]))))
write(scores,'NSCLC_CRC_defined_pathway_scores.csv')

rbc <- function(a,b) {
 n<-length(a);m<-length(b);r<-rank(c(a,b),ties.method='average')
 2*(sum(r[seq_len(n)])-n*(n+1)/2)/(n*m)-1
}
bootstrap <- function(a,b) {
 set.seed(42,kind='Mersenne-Twister',normal.kind='Inversion',sample.kind='Rejection')
 replicate(5000,rbc(sample(a,length(a),replace=TRUE),sample(b,length(b),replace=TRUE)))
}
one <- function(v,id) {
 a<-v[nr];b<-v[nnr];eff<-rbc(a,b)
 wt<-suppressWarnings(wilcox.test(a,b,exact=FALSE,correct=TRUE))
 stopifnot(abs(eff-mean(outer(a,b,function(x,y)sign(x-y))))<1e-12,
  abs(eff-(2*unname(wt$statistic)/(length(a)*length(b))-1))<1e-12)
 boot<-bootstrap(a,b);ci<-unname(quantile(boot,c(.025,.975),type=7))
 write(data.frame(replicate=seq_along(boot),effect=boot),paste0('bootstrap_',id,'.csv'))
 data.frame(ID=id,n_R=length(a),n_NR=length(b),effect=eff,ci_low=ci[1],ci_high=ci[2],
  wilcox_p=wt$p.value,mean_R=mean(a),mean_NR=mean(b),prevalence=mean(v>0),
  bootstrap_reps=5000L,seed=42L)
}
message('Computing seven pooled KO effects and group-stratified bootstrap intervals...')
ko <- do.call(rbind,lapply(target,function(k)one(ns[,k],k)))
ko$gene<-gene;ko$role<-role;ko$q<-p.adjust(ko$wilcox_p,'BH')
ko$log2FC_NR_R<-log2((ko$mean_NR+1e-8)/(ko$mean_R+1e-8))
ko$below_10pct<-ko$prevalence<.1
stopifnot(all(is.finite(as.matrix(ko[c('effect','ci_low','ci_high','wilcox_p','q')]))))
write(ko,'NSCLC_seven_KO_pooled_statistics.csv')
pathnames <- c('Degradation','Production','Production_Degradation')
pa<-do.call(rbind,lapply(pathnames,function(s)one(scores[[s]],s)));pa$q<-p.adjust(pa$wilcox_p,'BH')
write(pa,'NSCLC_pathway_pooled_statistics.csv')

message('Computing pooled CRC and NSCLC genome-wide rank-biserial effects...')
crc <- read(crcfile); stopifnot(nrow(crc)==1647,!anyDuplicated(crc$Run_ID),sum(crc$Analysis_Group=='Healthy')==745)
crcidx<-crc$Analysis_Group=='Healthy'; cv<-as.matrix(crc[grep('^K[0-9]{5}$',names(crc))])
stopifnot(all(is.finite(cv)),all(cv>=0))
ct<-read(crctarget); ix<-match(crc$Run_ID,ct$Run.ID)
stopifnot(!anyNA(ix),max(abs(cv[,target]-as.matrix(ct[ix,target])))<1e-12)
frozen<-read(crcfrozen);frozen<-frozen[match(target,frozen$KO),]
crc7<-data.frame(KO=target,effect=sapply(target,function(k)rbc(cv[crcidx,k],cv[!crcidx,k])),
  p=sapply(target,function(k)suppressWarnings(wilcox.test(cv[crcidx,k],cv[!crcidx,k],exact=FALSE)$p.value)))
stopifnot(max(abs(crc7$p-frozen$p_value))<1e-12,
 max(abs(p.adjust(crc7$p,'BH')-frozen$p_adj))<1e-12,
 max(abs(colMeans(cv[crcidx,target])-frozen$Mean_Healthy))<1e-12)
write(crc7,'CRC_target_reconstruction.csv')
ncprev<-colMeans(ns>0);ccprev<-colMeans(cv>0)
common<-intersect(colnames(ns)[ncprev>=.1],colnames(cv)[ccprev>=.1])
plotkos<-sort(union(common,target))
cross<-data.frame(KO=plotkos,CRC_effect=sapply(plotkos,function(k)rbc(cv[crcidx,k],cv[!crcidx,k])),
  NSCLC_effect=sapply(plotkos,function(k)rbc(ns[nr,k],ns[nnr,k])),
  CRC_prevalence=ccprev[plotkos],NSCLC_prevalence=ncprev[plotkos],genomewide_background=plotkos %in% common,
  targeted=plotkos %in% target)
cross$gene<-gene[match(cross$KO,target)];cross$role<-role[match(cross$KO,target)]
stopifnot(sum(cross$targeted)==7,all(is.finite(cross$CRC_effect)),all(is.finite(cross$NSCLC_effect)))
write(cross,'CRC_NSCLC_pooled_KO_concordance.csv')
bg<-cross[cross$genomewide_background,]
# Descriptive concordance. No copied P value from the old meta-analysis.
concordance<-data.frame(n_background=nrow(bg),n_display=nrow(cross),n_targets=7,
  spearman_rho=cor(bg$CRC_effect,bg$NSCLC_effect,method='spearman'),
  same_direction=sum(sign(bg$CRC_effect)==sign(bg$NSCLC_effect)),
  direction_agreement_percent=100*mean(sign(bg$CRC_effect)==sign(bg$NSCLC_effect)))
write(concordance,'concordance_summary.csv')

message('Rendering three PNG figures...')
cols<-c('R higher'='#2E5F8A','NR higher'='#C47B3B')
theme_plot<-function()theme_classic(base_family='Arial',base_size=9)+theme(
 axis.text=element_text(colour='#293944'),plot.title=element_text(size=12,face='bold',colour='#263640'),
 plot.subtitle=element_text(size=8.5,colour='#61707A'),plot.caption=element_text(size=7,hjust=0,colour='#61707A',lineheight=1.1),
 plot.margin=margin(8,10,8,8),legend.position='none')
fmt<-function(p) ifelse(p<.001,sprintf('%.2e',p),sprintf('%.3f',p))
ko$y<-7:1;ko$direction<-ifelse(ko$effect>=0,'R higher','NR higher')
ko$fill<-ifelse(ko$below_10pct,'white',cols[ko$direction])
ko$label<-paste0(c('Sarcosine oxidase','<i>soxA</i>','<i>soxB</i>','<i>soxG</i>','PIPOX','DMGDH','Creatinase'),
  " <span style='color:#7A858D'>(",target,')</span>')
xf<-max(abs(c(ko$ci_low,ko$ci_high)))*1.2
forest<-ggplot(ko,aes(x=effect,y=y))+
 geom_hline(yintercept=2.5,colour='#DFE3E6',linewidth=.35)+
 geom_vline(xintercept=0,colour='#AAB3B9',linetype='22',linewidth=.4)+
 geom_segment(aes(x=ci_low,xend=ci_high,yend=y,colour=direction),linewidth=.85)+
 geom_point(aes(fill=fill,colour=direction),shape=21,size=3,stroke=.8)+
 scale_fill_identity()+scale_colour_manual(values=cols)+
 scale_x_continuous(limits=c(-xf,xf),breaks=scales::breaks_pretty(n=5))+
 scale_y_continuous(limits=c(.5,7.95),breaks=ko$y,labels=ko$label,expand=expansion(mult=0))+
 labs(x='Pooled rank-biserial effect (R − NR)',y=NULL)+theme_plot()+
 theme(axis.text.y=ggtext::element_markdown(size=8.5),axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(3,2,3,0))
right<-ggplot(ko,aes(y=y))+
 geom_hline(yintercept=2.5,colour='#DFE3E6',linewidth=.35)+
 geom_text(aes(x=.25,label=sprintf('%.1f',100*prevalence)),family='Arial',size=2.8)+
 geom_text(aes(x=1.15,label=fmt(q)),family='Arial',size=2.8)+
 annotate('text',x=c(.25,1.15),y=7.65,label=c('Prev. (%)','BH q'),family='Arial',fontface='bold',size=2.8)+
 scale_y_continuous(limits=c(.5,7.95),expand=expansion(mult=0))+scale_x_continuous(limits=c(-.2,1.7))+theme_void()+theme(plot.margin=margin(3,0,3,0))
fig1<-forest+right+plot_layout(widths=c(4,1.25))+plot_annotation(
 title='Individual microbial sarcosine genes',subtitle='Three NSCLC cohorts pooled · 432 R / 392 NR Run records',
 caption='95% CI: 5,000 bootstrap resamples. BH correction across 7 KOs.\nOpen circles: pooled prevalence <10%. Unadjusted Run-level comparison.',theme=theme_plot())

pa$y<-3:1;pa$direction<-ifelse(pa$effect>=0,'R higher','NR higher')
pa$label<-c('Degradation\n5-KO sum','Production\n2-KO sum','Production / Degradation\nlog₂ ratio')
xx<-max(abs(c(pa$ci_low,pa$ci_high)))*1.35
pp<-ggplot(pa,aes(x=effect,y=y))+
 geom_vline(xintercept=0,colour='#AAB3B9',linetype='22',linewidth=.4)+
 geom_segment(aes(x=ci_low,xend=ci_high,yend=y,colour=direction),linewidth=.85)+
 geom_point(aes(colour=direction),size=3)+scale_colour_manual(values=cols)+
 scale_x_continuous(limits=c(-xx,xx),breaks=scales::breaks_pretty(n=4))+
 scale_y_continuous(limits=c(.5,3.8),breaks=pa$y,labels=pa$label,expand=expansion(mult=0))+
 labs(x='Pooled rank-biserial effect (R − NR)',y=NULL)+theme_plot()+
 theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=8.5),plot.margin=margin(3,2,3,0))
pr<-ggplot(pa,aes(y=y))+
 geom_text(aes(x=.4,label=fmt(wilcox_p)),family='Arial',size=2.8)+
 geom_text(aes(x=1.5,label=fmt(q)),family='Arial',size=2.8)+
 annotate('text',x=c(.4,1.5),y=3.65,label=c('P','BH q'),family='Arial',size=2.8,fontface='bold')+
 scale_y_continuous(limits=c(.5,3.8),expand=expansion(mult=0))+scale_x_continuous(limits=c(0,2))+theme_void()+theme(plot.margin=margin(3,0,3,0))
fig2<-pp+pr+plot_layout(widths=c(4,1.5))+plot_annotation(title='Sarcosine pathway balance',
 subtitle='CRC-defined seven-KO set · Three NSCLC cohorts pooled',
 caption='95% bootstrap CI. Two-sided Wilcoxon P. BH correction across 3 scores.\nRatio = log₂[(production + 10⁻⁸)/(degradation + 10⁻⁸)]. Unadjusted Run-level comparison.',theme=theme_plot())

hi<-cross[cross$targeted,];hi<-hi[match(target,hi$KO),]
hi$pointfill<-ifelse(pmin(hi$CRC_prevalence,hi$NSCLC_prevalence)<.1,'white',ifelse(hi$role=='Degradation','#2E5F8A','#C47B3B'))
li<-ceiling(max(abs(c(cross$CRC_effect,cross$NSCLC_effect)))*10)/10
cross$label<-ifelse(cross$targeted,paste0(cross$gene,' (',cross$KO,')'),'')
fig3<-ggplot(cross,aes(CRC_effect,NSCLC_effect))+
 geom_hline(yintercept=0,colour='#BAC2C7',linewidth=.35,linetype='22')+
 geom_vline(xintercept=0,colour='#BAC2C7',linewidth=.35,linetype='22')+
 geom_abline(slope=1,colour='#D8DDE0',linewidth=.3,linetype='dotted')+
 geom_point(data=cross[!cross$targeted,],colour='#91A0AA',alpha=.25,size=.6,stroke=0)+
 geom_point(data=hi,aes(fill=pointfill,colour=role),shape=21,size=3,stroke=.8)+
 scale_fill_identity()+scale_colour_manual(values=c(Degradation='#2E5F8A',Production='#C47B3B'))+
 geom_text_repel(aes(label=label),seed=42,family='Arial',size=2.9,colour='#2D3D47',
  point.size=.6,point.padding=.2,box.padding=.3,min.segment.length=0,max.overlaps=Inf,
  max.iter=10000,max.time=10,segment.colour='#71808A',segment.size=.3)+
 annotate('text',x=-.96*li,y=.96*li,hjust=0,vjust=1,family='Arial',size=2.8,lineheight=1.15,
  label=sprintf('%s shared prevalent KOs\nSpearman ρ = %.3f\nSame direction: %.1f%%',format(nrow(bg),big.mark=','),concordance$spearman_rho,concordance$direction_agreement_percent))+
 coord_equal(xlim=c(-li,li),ylim=c(-li,li))+
 labs(title='Cross-disease concordance of microbial KO effects',subtitle='All seven CRC-defined sarcosine KOs highlighted',
 x='CRC pooled rank-biserial effect (Healthy − Cancer)',y='NSCLC pooled rank-biserial effect (R − NR)',
 caption='Blue: degradation. Orange: production. Open circles: prevalence <10% in either disease.\nBackground and concordance: KOs with pooled prevalence ≥10% in both diseases.\nUnadjusted Run-level comparisons. Low-prevalence targets retained separately.')+theme_plot()+theme(axis.title=element_text(size=8.5))

save_png<-function(p,name,w,h) {
 seed<-.Random.seed; temp<-tempfile(fileext='.png')
 ggsave(temp,p,width=w,height=h,dpi=600,units='in',device=ragg::agg_png,bg='white')
 .Random.seed<<-seed
 dst<-file.path(out,'figures',name);stopifnot(file.copy(temp,dst,overwrite=TRUE),file.exists(dst),file.info(dst)$size>10000);unlink(temp)
}
save_png(fig1,'NSCLC_CRC7_Individual_KOs.png',6.6,3.9)
save_png(fig2,'NSCLC_CRC7_Pathway_Balance.png',6.6,2.85)
save_png(fig3,'CRC_NSCLC_CRC7_Concordance.png',5.8,5.95)
stopifnot(identical(hashes,vapply(inputs,hash,'')))
writeLines(capture.output(sessionInfo()),file.path(out,'qa','sessionInfo.txt'))
writeLines(c('# CRC-method seven-KO NSCLC branch','',
 'Author explicitly requested the CRC pooled analysis strategy on 2026-09-07, correcting the rejected prevalence-gated meta redraw.',
 'Three discovery NSCLC cohorts only. 824 unique Run IDs: 432 R, 392 NR. These are not asserted to be 824 independent patients.',
 'Original exported Sample names overlap across two projects (283 repeated names). No silent deduplication was performed.',
 'This is an exploratory Run-level pooled comparison analogous to the CRC code, without cohort/covariate adjustment. It is NOT the old MaAsLin2/DL meta-analysis.',
 'Seven target KOs retained regardless of prevalence: K00301 K00302 K00303 K00305 K00306 K00315 K08688.',
 'Source zeros are retained. Missing cohort KO columns are zero-filled only within documented complete KO profiles. No patient, feature or value was fabricated.',
 'PIPOX uses all selected Run records, including source zero observations in two cohorts. Its pooled contrast is not independent cross-cohort replication.',
 'Effect: rank-biserial P(R>NR)-P(R<NR); 95% percentile group-stratified bootstrap, 5000 replicates, per-feature seed 42, quantile type 7.',
 'P: two-sided Wilcoxon rank-sum, normal approximation with ties and continuity correction. BH across seven target KO tests; separately across three score tests.',
 'Degradation = first five target KO sums; Production = last two. Ratio copies CRC: log2((Production+1e-8)/(Degradation+1e-8)). All 824 ratios are retained.',
 'CRC concordance reference uses its original 1647 Run records (745 Healthy/902 Cancer), matching frozen seven-KO values.',
 'Concordance uses newly computed pooled rank-biserial effects on both axes. Background/rho use >=10% pooled prevalence in both diseases; the seven targets are an explicit display override.',
 'No old 7400-KO meta rho/P or I2 was carried into these pooled graphics. No permutation P is reported for this descriptive concordance.',
 'Original data, figures, manuscript and PPT remain unchanged. These PNGs are new unadopted candidates. See numerical and raster QA for completed checks.'),file.path(out,'ANALYSIS_NOTES.md'))
cat('COMPLETE R ANALYSIS:',out,'\n')
