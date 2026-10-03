#!/usr/bin/env Rscript
# Author-requested individual-cohort extension of the verified CRC-method pooled DA.
# Frozen design: analyze the same three discovery Run rosters separately.
# Within each cohort: raw MetaPhlAn4 percentage abundance; prevalence >=10%;
# x=log2((arithmetic mean_NR+1e-6)/(arithmetic mean_R+1e-6)).
# Two-sided unpaired Wilcoxon, exact=FALSE, correct=TRUE, including tie correction.
# BH across all prevalence-eligible species WITHIN each cohort, not the pooled set.
# Display q<.05 AND |log2FC|>1. Negative/teal=R, positive/red=NR.
# Bars: up to 15 per direction by q; volcano labels: up to 15 by q.
# Exact q ties ordered by full taxonomy ID. Never pad an empty/short result.
# No covariate adjustment, relabeling, deduplication or outlier removal.
# Unique Run and Sample-name strings within each cohort do not establish that
# different cohorts are independent replications: cross-project overlap remains.
suppressPackageStartupMessages({
  library(data.table);library(digest);library(jsonlite);library(stringi)
})
setDTthreads(2);options(stringsAsFactors=FALSE)
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
OUT<-dirname(dirname(normalizePath(sub('^--file=','',arg[1]))))
POOLED<-dirname(OUT);BASE<-dirname(dirname(dirname(POOLED)))
PROJECT<-dirname(dirname(BASE));QA<-file.path(OUT,'qa')
dir.create(QA,recursive=TRUE,showWarnings=FALSE)
USED<-list.dirs(BASE,recursive=FALSE,full.names=TRUE)
USED<-USED[stri_trans_nfc(basename(USED))=='사용데이터_모음'];stopifnot(length(USED)==1L)
COHORTS<-c('PRJNA751792','PRJNA1023797','PRJEB22863')
NRUN<-c(338L,421L,65L);NR<-c(174L,225L,33L);NNR<-c(164L,196L,32L)
hash<-function(p)digest(file=p,algo='sha256',serialize=FALSE)
manifest<-fread(file.path(POOLED,'qa/input_sha256.csv'))
stopifnot(all(vapply(file.path(PROJECT,manifest$path),hash,'')==manifest$sha256))
# Protect all completed pooled outputs/scripts/QA, as well as the source inputs.
protected<-list.files(POOLED,recursive=TRUE,full.names=TRUE)
protected<-protected[!startsWith(protected,paste0(OUT,.Platform$file.sep))]
extra<-data.table(path=substring(protected,nchar(PROJECT)+2L),
  sha256=vapply(protected,hash,''),bytes=file.info(protected)$size)
manifest<-unique(rbind(manifest,extra),by='path')
fwrite(manifest,file.path(QA,'input_sha256.csv'))
write_json(list(
  design='Same CRC-method contrast and thresholds, applied independently within each cohort',
  cohorts=COHORTS,unit='Run record',raw_abundance_units='percent',
  prevalence_min=.10,prevalence_scope='Within each cohort',
  contrast='log2((Mean_NR+1e-6)/(Mean_R+1e-6)); R negative, NR positive',
  test='Two-sided unpaired Wilcoxon rank sum; exact=FALSE, correct=TRUE, tie correction',
  BH_family='All prevalence-eligible species within the corresponding cohort',
  display_q_less_than=.05,display_abs_log2FC_greater_than=1,
  bars='Up to 15 per direction by q, ordered by log2FC',
  labels='Up to 15 meeting both display thresholds, by q',
  covariates='None',source_audit='Previously audited original/cached inputs; all source hashes rechecked',
  limitation='Cross-project Sample-name overlap and conflicting labels remain unresolved; do not claim three independent patient cohorts.'
),file.path(QA,'analysis_design.json'),pretty=TRUE,auto_unbox=TRUE)
roster<-fread(file.path(POOLED,'qa/sample_metadata.csv'))
summaries<-list();eligibility<-list()
for(i in seq_along(COHORTS)) {
  co<-COHORTS[i];dest<-file.path(OUT,co);dir.create(dest,showWarnings=FALSE)
  dat<-fread(file.path(USED,co,paste0(co,'_species_relative_abundance_matrix.csv')))
  md<-fread(file.path(USED,co,paste0(co,'_patient_metadata_WGS_NSCLC_R_NR.csv')))
  stopifnot(nrow(dat)==NRUN[i],nrow(md)==NRUN[i],!anyDuplicated(dat$Run_ID),
    !anyDuplicated(md[['Run ID']]),!anyDuplicated(md[['Sample name']]),
    setequal(dat$Run_ID,md[['Run ID']]),all(md$Cohort==co),
    all(md[['Assay type']]=='WGS'),all(md[['Phenotype name']]=='Carcinoma, Non-Small-Cell Lung'))
  md<-md[match(dat$Run_ID,md[['Run ID']])]
  old<-roster[Cohort==co][match(dat$Run_ID,Run_ID)]
  stopifnot(identical(old$Run_ID,dat$Run_ID),
    identical(old$Group,md[['Analysis Group']]),
    identical(old$Sample_name,md[['Sample name']]))
  m<-as.matrix(dat[,-1]);storage.mode(m)<-'double'
  stopifnot(all(is.finite(m)),all(m>=0),max(abs(rowSums(m)-100))<.001,
    !anyDuplicated(colnames(m)),all(grepl('|s__',colnames(m),fixed=TRUE)),
    !any(grepl('|t__',colnames(m),fixed=TRUE)))
  ri<-which(md[['Analysis Group']]=='R');ni<-which(md[['Analysis Group']]=='NR')
  stopifnot(length(ri)==NR[i],length(ni)==NNR[i],length(ri)+length(ni)==nrow(m))
  prev<-colMeans(m>0);keep<-prev>=.10
  eligibility[[co]]<-data.table(Cohort=co,Species_full=colnames(m),
    nonzero_n=colSums(m>0),n=nrow(m),prevalence=prev,eligible=keep)
  d<-m[,keep,drop=FALSE]
  # No undefined constant-feature test is silently assigned a fabricated P.
  stopifnot(all(apply(d,2,function(x)length(unique(x))>1L)))
  tt<-lapply(seq_len(ncol(d)),function(j)wilcox.test(d[ni,j],d[ri,j],
    alternative='two.sided',paired=FALSE,exact=FALSE,correct=TRUE))
  ans<-data.table(Cohort=co,Species_full=colnames(d),Species=sub('.*\\|s__','',colnames(d)),
    n_R=length(ri),n_NR=length(ni),prevalence=prev[keep],
    prevalence_R=colMeans(d[ri,,drop=FALSE]>0),prevalence_NR=colMeans(d[ni,,drop=FALSE]>0),
    Mean_R=colMeans(d[ri,,drop=FALSE]),Mean_NR=colMeans(d[ni,,drop=FALSE]),
    Wilcoxon_U_NR=vapply(tt,function(x)as.numeric(x$statistic),numeric(1)),
    p_value=vapply(tt,function(x)x$p.value,numeric(1)))
  stopifnot(!anyDuplicated(ans$Species),all(is.finite(ans$p_value)),
    all(ans$p_value>0 & ans$p_value<=1))
  ans[,log2FC_NR_vs_R:=log2((Mean_NR+1e-6)/(Mean_R+1e-6))]
  ans[,BH_q:=p.adjust(p_value,method='BH')]
  ans[,neg_log10_BH_q:=-log10(BH_q)]
  ans[,Direction:=fifelse(log2FC_NR_vs_R<0,'R higher',fifelse(log2FC_NR_vs_R>0,'NR higher','Equal'))]
  ans[,Display_class:=fifelse(BH_q<.05 & log2FC_NR_vs_R< -1,'R higher',
    fifelse(BH_q<.05 & log2FC_NR_vs_R>1,'NR higher','Below display thresholds'))]
  setorder(ans,BH_q,Species_full)
  bar<-rbindlist(lapply(c('R higher','NR higher'),function(z)head(ans[Display_class==z],15L)))
  lab<-head(ans[Display_class!='Below display thresholds'],15L)
  setorder(bar,log2FC_NR_vs_R,Species_full)
  ans[,volcano_label:=Species_full %in% lab$Species_full]
  ans[,bar_selected:=Species_full %in% bar$Species_full]
  bar[,Species_display:=gsub('_',' ',Species)]
  bar[,bar_order_bottom_to_top:=seq_len(.N)]
  fwrite(ans,file.path(dest,paste0(co,'_species_DA_results.csv')))
  fwrite(bar,file.path(dest,paste0(co,'_species_DA_bar_data.csv')))
  summaries[[co]]<-data.table(Cohort=co,n_runs=nrow(d),n_R=length(ri),n_NR=length(ni),
    input_species=ncol(m),tested_species=ncol(d),q_below_05=sum(ans$BH_q<.05),
    R_display_hits=sum(ans$Display_class=='R higher'),NR_display_hits=sum(ans$Display_class=='NR higher'),
    bar_R=sum(bar$Display_class=='R higher'),bar_NR=sum(bar$Display_class=='NR higher'),
    volcano_labels=nrow(lab),min_q=min(ans$BH_q),
    min_log2FC=min(ans$log2FC_NR_vs_R),max_log2FC=max(ans$log2FC_NR_vs_R))
}
fwrite(rbindlist(eligibility),file.path(QA,'species_prevalence.csv'))
fwrite(rbindlist(summaries),file.path(OUT,'cohort_summary.csv'))
stopifnot(all(vapply(file.path(PROJECT,manifest$path),hash,'')==manifest$sha256))
writeLines(capture.output(sessionInfo()),file.path(QA,'sessionInfo.txt'))
print(rbindlist(summaries))

