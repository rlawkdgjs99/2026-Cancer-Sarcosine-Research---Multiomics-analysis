#!/usr/bin/env Rscript
# Author-requested CRC-method species DA, three NSCLC discovery cohorts.
# Frozen BEFORE testing (2026-09-09):
# - Unit: existing 824 Run records; R/NR labels preserved from source metadata.
# - Exclude Korean validation and RCC. No biological-subject deduplication.
# - Raw MetaPhlAn4 species relative abundance (%) including zero observations.
# - Union of full species taxonomy IDs, absent entries zero; pooled prevalence >=10%.
# - x = log2((arithmetic mean_NR + 1e-6)/(arithmetic mean_R + 1e-6)).
# - Two-sided unpaired Wilcoxon rank-sum, asymptotic, ties and continuity correction.
# - BH across ALL prevalence-eligible species. No cohort or covariate adjustment.
# - Display class: q<.05 AND |log2FC|>1; R negative/teal, NR positive/red.
# - Bar: up to 15 species per direction by smallest q; then order by log2FC.
# - Volcano labels: up to 15 display-eligible species by smallest q.
# - No outcome-dependent filters, effect cutoffs, pseudocounts or reranking rules.
# Known limitation: overlapping Sample names and discordant cross-project response
# labels remain in the legacy Run set. Pooled P/q do not establish independent-
# patient inference. This is a requested unadjusted pooled comparison.
# CRC sources: integrated_analysis_pooled.R lines348-405; regenerate_plots.R
# lines105-107,184-203; 99_publication_style_plots_26.07.23.R lines598-605.

suppressPackageStartupMessages({
  library(data.table); library(digest); library(jsonlite)
})
setDTthreads(2); options(stringsAsFactors=FALSE); set.seed(42)
arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
SCRIPT <- normalizePath(sub('^--file=', '', arg[1]))
OUT <- dirname(dirname(SCRIPT)); QA <- file.path(OUT,'qa')
POOL <- dirname(dirname(OUT)); BASE <- dirname(POOL)
COLLECTION <- dirname(BASE); PROJECT <- dirname(COLLECTION)
CRC <- file.path(COLLECTION,'HGMT_CRC_WGS-Healthy_vs_Cancer','Healthy_vs_Cancer_4_CRC_cohorts_integrated')
USED <- list.dirs(BASE,recursive=FALSE,full.names=TRUE)
USED <- USED[stringi::stri_trans_nfc(basename(USED))=='사용데이터_모음']
stopifnot(length(USED)==1L)
DIRS <- c(PRJNA751792='NSCLC_PRJNA751792',PRJNA1023797='NSCLC_PRJNA1023797',PRJEB22863='NSCLC_RCC_PRJEB22863')
N_EXPECTED <- c(PRJNA751792=338L,PRJNA1023797=421L,PRJEB22863=65L)
PREVALENCE <- .10; PC <- 1e-6; Q_CUTOFF <- .05; LFC_CUTOFF <- 1
BAR_PER_DIRECTION <- 15L; VOLCANO_LABELS <- 15L
rel <- function(p) substring(normalizePath(p),nchar(PROJECT)+2L)
hash <- function(p) digest(file=p,algo='sha256',serialize=FALSE)
one <- function(d,pattern) {p<-list.files(d,pattern=pattern,full.names=TRUE);stopifnot(length(p)==1L);p}
read_meta <- function(p) {
  ln<-sub('\r$','',readLines(p,warn=FALSE));hdr<-strsplit(ln[2],'\t',fixed=TRUE)[[1]]
  rows<-ln[-c(1,2)];rows<-rows[nzchar(rows)]
  ans<-as.data.table(do.call(rbind,lapply(strsplit(rows,'\t',fixed=TRUE),function(x)x[seq_along(hdr)])))
  setnames(ans,hdr);ans
}
response <- function(s) vapply(s,function(x){
  if(is.na(x))return(NA_character_)
  v<-regmatches(x,regexpr('response_group:[^;]*',x))
  if(!length(v)||!nzchar(v))v<-regmatches(x,regexpr('response:[^;]*',x))
  if(!length(v)||!nzchar(v))return(NA_character_)
  v<-trimws(sub('^response(_group)?:[[:space:]]*','',v))
  if(v %in% c('R','NR'))v else NA_character_
},character(1))
inputs <- c(file.path(CRC,'scripts',c('integrated_analysis_pooled.R','regenerate_plots.R','99_publication_style_plots_26.07.23.R')),
            file.path(CRC,'results_integrated/bacteria/diff_abundance_species_pooled.csv'),
            file.path(POOL,'results/pooled_species_meta.csv'),
            file.path(POOL,'results/CRC_Method_7KO_NSCLC_26.09.07/tables/NSCLC_CRC_defined_pathway_scores.csv'))
for(co in names(DIRS)) inputs<-c(inputs,
  file.path(USED,co,paste0(co,c('_species_relative_abundance_matrix.csv','_patient_metadata_WGS_NSCLC_R_NR.csv'))),
  one(file.path(BASE,DIRS[co]),'^Bacteria_.*\\.txt$'),
  one(file.path(BASE,DIRS[co]),'^selected_project_.*\\.txt$'))
stopifnot(all(file.exists(inputs)));before<-vapply(inputs,hash,'')
manifest<-data.table(path=vapply(inputs,rel,''),sha256=before,bytes=file.info(inputs)$size)
matrices<-list();metas<-list();preflight<-list()
for(co in names(DIRS)) {
  cat('Checking original vs exported inputs:',co,'\n');flush.console()
  cached<-fread(file.path(USED,co,paste0(co,'_species_relative_abundance_matrix.csv')))
  md<-fread(file.path(USED,co,paste0(co,'_patient_metadata_WGS_NSCLC_R_NR.csv')),na.strings=c('','NA'))
  original<-read_meta(one(file.path(BASE,DIRS[co]),'^selected_project_.*\\.txt$'))
  original<-original[`Assay type`=='WGS' & `Phenotype name`=='Carcinoma, Non-Small-Cell Lung']
  original[,parsed_group:=response(`Sample description`)]
  original<-original[!is.na(parsed_group)]
  stopifnot(nrow(cached)==N_EXPECTED[co],nrow(md)==N_EXPECTED[co],nrow(original)==N_EXPECTED[co],
    !anyDuplicated(cached$Run_ID),!anyDuplicated(md[['Run ID']]),!anyDuplicated(original[['Run ID']]),
    setequal(cached$Run_ID,md[['Run ID']]),setequal(cached$Run_ID,original[['Run ID']]))
  md<-md[match(cached$Run_ID,md[['Run ID']])]
  original<-original[match(cached$Run_ID,original[['Run ID']])]
  stopifnot(identical(md[['Analysis Group']],original$parsed_group),
    identical(md[['Sample name']],original[['Sample name']]),all(md$Cohort==co),
    all(md[['Assay type']]=='WGS'),all(md[['Phenotype name']]=='Carcinoma, Non-Small-Cell Lung'))
  m<-as.matrix(cached[,-1]);storage.mode(m)<-'double';rownames(m)<-cached$Run_ID
  stopifnot(all(is.finite(m)),all(m>=0),max(abs(rowSums(m)-100))<.001,
    all(grepl('|s__',colnames(m),fixed=TRUE)),!any(grepl('|t__',colnames(m),fixed=TRUE)),
    !anyDuplicated(colnames(m)))
  raw<-fread(one(file.path(BASE,DIRS[co]),'^Bacteria_.*\\.txt$'),sep='\t',quote='',showProgress=FALSE)
  stopifnot(ncol(raw)==3L);setnames(raw,c('Taxa','Run_ID','Abundance'))
  raw<-raw[Run_ID %in% cached$Run_ID & grepl('|s__',Taxa,fixed=TRUE) & !grepl('|t__',Taxa,fixed=TRUE)]
  dup_n<-sum(duplicated(raw,by=c('Run_ID','Taxa')));stopifnot(dup_n==0L)
  rw<-dcast(raw,Run_ID~Taxa,value.var='Abundance',fill=0)
  stopifnot(setequal(rw$Run_ID,cached$Run_ID),all(names(rw)[-1] %in% colnames(m)))
  # Cached matrices retain columns detected only in excluded source records.
  # Such columns must be all-zero in this Run set; verify rather than dropping them.
  cache_only<-setdiff(colnames(m),names(rw)[-1])
  stopifnot(!length(cache_only)||all(m[,cache_only,drop=FALSE]==0))
  rm<-matrix(0,nrow(m),ncol(m),dimnames=dimnames(m))
  rm[,names(rw)[-1]]<-as.matrix(rw[match(cached$Run_ID,rw$Run_ID),-1])
  raw_difference<-max(abs(m-rm));stopifnot(raw_difference<1e-12)
  matrices[[co]]<-m
  metas[[co]]<-data.table(Cohort=co,Run_ID=cached$Run_ID,Sample_name=md[['Sample name']],Group=md[['Analysis Group']])
  preflight[[co]]<-data.table(Cohort=co,n=nrow(m),n_R=sum(md[['Analysis Group']]=='R'),n_NR=sum(md[['Analysis Group']]=='NR'),
    species=ncol(m),missing=sum(is.na(m)),zero_values=sum(m==0),row_sum_min=min(rowSums(m)),row_sum_max=max(rowSums(m)),
    raw_export_max_abs_difference=raw_difference,duplicate_Run_taxa=dup_n,cache_only_zero_columns=length(cache_only))
}
metadata<-rbindlist(metas);stopifnot(nrow(metadata)==824L,!anyDuplicated(metadata$Run_ID),
  sum(metadata$Group=='R')==432L,sum(metadata$Group=='NR')==392L)
roster<-fread(file.path(POOL,'results/CRC_Method_7KO_NSCLC_26.09.07/tables/NSCLC_CRC_defined_pathway_scores.csv'))
stopifnot(setequal(roster$Run_ID,metadata$Run_ID))
ri<-match(metadata$Run_ID,roster$Run_ID)
stopifnot(identical(metadata$Group,roster$Group[ri]),identical(metadata$Cohort,roster$Cohort[ri]))
ids<-unique(unlist(lapply(matrices,colnames),use.names=FALSE))
short<-sub('.*\\|s__','',ids);stopifnot(!anyDuplicated(short))
m<-matrix(0,nrow(metadata),length(ids),dimnames=list(metadata$Run_ID,ids))
for(co in names(DIRS))m[rownames(matrices[[co]]),colnames(matrices[[co]])]<-matrices[[co]]
stopifnot(identical(rownames(m),metadata$Run_ID),all(is.finite(m)),max(abs(rowSums(m)-100))<.001)
overlap<-metadata[,.(n_runs=.N,n_cohorts=uniqueN(Cohort),n_groups=uniqueN(Group)),by=Sample_name]
metadata<-merge(metadata,overlap,by='Sample_name',all.x=TRUE,sort=FALSE)
metadata<-metadata[match(rownames(m),Run_ID)]
metadata[,cross_project_name_overlap:=n_cohorts>1L]
metadata[,cross_project_group_conflict:=n_cohorts>1L & n_groups>1L]
stopifnot(nrow(overlap[n_cohorts>1L])==283L,nrow(overlap[n_cohorts>1L & n_groups>1L])==45L)
prev<-colMeans(m>0);keep<-prev>=PREVALENCE
cat('Pooled:',nrow(m),'Run records;',ncol(m),'union species;',sum(keep),'pass prevalence\n')
plan<-list(frozen_design_date='2026-09-09',cohorts=names(DIRS),unit='Run record',n_R=432,n_NR=392,
  abundance_units='MetaPhlAn4 relative abundance percent; no renormalization or log transformation for Wilcoxon',
  pooled_prevalence_min=PREVALENCE,pseudocount_percent_units=PC,
  contrast='log2((mean_NR+1e-6)/(mean_R+1e-6)); negative R higher, positive NR higher',
  test='Unpaired two-sided Wilcoxon rank-sum; exact=FALSE, correct=TRUE; tie correction',
  covariates='None; cohort unadjusted',BH_family='All pooled prevalence-eligible species',
  display_q_less_than=Q_CUTOFF,display_abs_log2FC_greater_than=LFC_CUTOFF,
  bar_selection='Up to 15 per direction, lowest BH q; order by log2FC',
  volcano_labels='Up to 15 display-eligible species, lowest BH q',
  limitation='Legacy Run records retain 283 cross-project repeated Sample names, including 45 R/NR conflicts. Independent-patient inference is not established.')
write_json(plan,file.path(QA,'analysis_design.json'),pretty=TRUE,auto_unbox=TRUE)
fwrite(manifest,file.path(QA,'input_sha256.csv'))
fwrite(rbindlist(preflight),file.path(QA,'input_verification.csv'))
fwrite(metadata,file.path(QA,'sample_metadata.csv'))
fwrite(data.table(Species_full=ids,Species=short,prevalence=prev,eligible=keep),file.path(QA,'species_prevalence.csv'))
if('--preflight-only' %in% commandArgs(TRUE))quit(status=0)

# Analysis begins only after all scientific input checks and design recording.
d<-m[,keep,drop=FALSE];ri<-which(metadata$Group=='R');ni<-which(metadata$Group=='NR')
tests<-lapply(seq_len(ncol(d)),function(j)stats::wilcox.test(d[ni,j],d[ri,j],
  alternative='two.sided',paired=FALSE,exact=FALSE,correct=TRUE))
ans<-data.table(Species_full=colnames(d),Species=sub('.*\\|s__','',colnames(d)),
  n_R=length(ri),n_NR=length(ni),prevalence=prev[keep],
  prevalence_R=colMeans(d[ri,,drop=FALSE]>0),prevalence_NR=colMeans(d[ni,,drop=FALSE]>0),
  Mean_R=colMeans(d[ri,,drop=FALSE]),Mean_NR=colMeans(d[ni,,drop=FALSE]),
  Wilcoxon_U_NR=vapply(tests,function(x)as.numeric(x$statistic),numeric(1)),
  p_value=vapply(tests,function(x)x$p.value,numeric(1)))
stopifnot(all(is.finite(ans$p_value)),all(ans$p_value>0 & ans$p_value<=1))
ans[,log2FC_NR_vs_R:=log2((Mean_NR+PC)/(Mean_R+PC))]
ans[,BH_q:=p.adjust(p_value,method='BH')]
ans[,neg_log10_BH_q:=-log10(BH_q)]
ans[,Direction:=fifelse(log2FC_NR_vs_R<0,'R higher',fifelse(log2FC_NR_vs_R>0,'NR higher','Equal'))]
ans[,Display_class:=fifelse(BH_q<Q_CUTOFF & log2FC_NR_vs_R< -LFC_CUTOFF,'R higher',
  fifelse(BH_q<Q_CUTOFF & log2FC_NR_vs_R>LFC_CUTOFF,'NR higher','Below display thresholds'))]
# Resolve exact q ties deterministically, without a secondary biological selection rule.
setorder(ans,BH_q,Species_full)
bar<-rbindlist(lapply(c('NR higher','R higher'),function(z)head(ans[Display_class==z],BAR_PER_DIRECTION)))
lab<-head(ans[Display_class!='Below display thresholds'],VOLCANO_LABELS)
setorder(bar,log2FC_NR_vs_R,Species_full)
ans[,volcano_label:=Species_full %in% lab$Species_full]
ans[,bar_selected:=Species_full %in% bar$Species_full]
bar[,Species_display:=gsub('_',' ',Species)]
bar[,bar_order_bottom_to_top:=seq_len(.N)]
fwrite(ans,file.path(OUT,'NSCLC_species_DA_results.csv'))
fwrite(bar,file.path(OUT,'NSCLC_species_DA_bar_data.csv'))
saveRDS(list(metadata=metadata,matrix=d,statistics=ans),'/private/tmp/nsclc_species_da_260909/analysis_cache.rds')
summary<-list(n_runs=nrow(d),n_R=length(ri),n_NR=length(ni),union_species=ncol(m),tested_species=ncol(d),
  BH_q_below_05=sum(ans$BH_q<Q_CUTOFF),
  R_display_hits=sum(ans$Display_class=='R higher'),NR_display_hits=sum(ans$Display_class=='NR higher'),
  bar_R=sum(bar$Display_class=='R higher'),bar_NR=sum(bar$Display_class=='NR higher'),
  volcano_labels=nrow(lab),minimum_q=min(ans$BH_q),log2FC_range=range(ans$log2FC_NR_vs_R),
  source_files_unchanged=identical(before,vapply(inputs,hash,'')))
stopifnot(summary$source_files_unchanged)
write_json(summary,file.path(QA,'analysis_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
writeLines(capture.output(sessionInfo()),file.path(QA,'sessionInfo.txt'))
print(summary)
