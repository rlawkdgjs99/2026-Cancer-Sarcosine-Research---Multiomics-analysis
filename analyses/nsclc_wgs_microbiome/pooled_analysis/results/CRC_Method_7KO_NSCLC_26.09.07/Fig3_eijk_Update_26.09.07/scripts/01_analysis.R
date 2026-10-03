#!/usr/bin/env Rscript
# Frozen statistical design: see ../ANALYSIS_PLAN_FROZEN.md.
# Uses only relative paths derived from this script. Seed42; no new cohorts.
suppressPackageStartupMessages({library(data.table);library(metafor);library(digest);library(jsonlite)})
set.seed(42); options(stringsAsFactors=FALSE); setDTthreads(2)
arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
script <- normalizePath(sub('^--file=','',arg[1]))
OUT_DIR <- dirname(dirname(script)); parent <- dirname(OUT_DIR)
pool <- dirname(dirname(parent)); NSCLC_ROOT <- dirname(pool)
COLLECTION_ROOT <- dirname(NSCLC_ROOT); PROJECT_ROOT <- dirname(COLLECTION_ROOT)
CRC_ROOT <- file.path(COLLECTION_ROOT,'HGMT_CRC_WGS-Healthy_vs_Cancer','Healthy_vs_Cancer_4_CRC_cohorts_integrated')
ENA_DIR <- file.path(pool,'resources','ENA_CRC_run_sample_mapping_26.08.23')
PREV_MIN <- .10; N_PERM <- 1000L
DEG_KOS <- c('K00301','K00302','K00303','K00305','K00306'); PROD_KOS <- c('K00315','K08688')
NSCLC_DISC_DIRS <- c(PRJNA751792='NSCLC_PRJNA751792',PRJNA1023797='NSCLC_PRJNA1023797',PRJEB22863='NSCLC_RCC_PRJEB22863')
CRC_DIRS <- c(PRJEB10878='PRJEB10878_CRC',PRJEB27928='PRJEB27928_CRC',PRJEB6070='PRJEB6070_CRC_AdenomatousPolyps',PRJNA429097='PRJNA429097_CRC')
source(file.path(dirname(script),'00_verified_helpers.R'))
say <- function(...) {cat(sprintf(...),'\n');flush.console()}
rd <- function(p) as.data.frame(fread(p,check.names=FALSE))
wr <- function(x,n) fwrite(as.data.table(x),file.path(OUT_DIR,'tables',paste0(n,'.csv')),na='')
used <- list.dirs(NSCLC_ROOT,recursive=FALSE,full.names=TRUE)
used <- used[stringi::stri_trans_nfc(basename(used))=='사용데이터_모음'];stopifnot(length(used)==1L)
scores_path <- file.path(parent,'tables','NSCLC_CRC_defined_pathway_scores.csv')
ko_path <- file.path(parent,'tables','NSCLC_seven_KO_per_Run.csv')
score_stats_path <- file.path(parent,'tables','NSCLC_pathway_pooled_statistics.csv')
da_path <- file.path(pool,'results','pooled_species_meta.csv')
KO_CACHE_FILE <- file.path(CRC_ROOT,'results_integrated','sarcosine','sarcosine_KO_per_sample_pooled.csv')
source_files <- c(scores_path,ko_path,score_stats_path,da_path,KO_CACHE_FILE)
hash <- function(p) digest(file=p,algo='sha256',serialize=FALSE)
rel <- function(p) substring(normalizePath(p),nchar(PROJECT_ROOT)+2L)
protected <- list.files(parent,recursive=TRUE,full.names=TRUE)
protected <- protected[!startsWith(protected,paste0(OUT_DIR,'/')) & !dir.exists(protected)]
md <- list.dirs(PROJECT_ROOT,recursive=FALSE,full.names=TRUE)
md <- md[startsWith(stringi::stri_trans_nfc(basename(md)),'Manuscript')]
protected <- c(protected,list.files(file.path(md,'FigDesign&Manuscript'),pattern='26.09.07\\.(pptx|docx)$',full.names=TRUE))
protected <- protected[!startsWith(basename(protected),'~$')]
protected_sha <- vapply(protected,hash,'')
say('[1] Verify the exact seven-KO exposure and sample identities')
sc <- rd(scores_path); ko <- rd(ko_path)
stopifnot(nrow(sc)==824L,!anyDuplicated(sc$Run_ID),identical(sc$Run_ID,ko$Run_ID),identical(sc$Group,ko$Group),
          sum(sc$Group=='R')==432L,sum(sc$Group=='NR')==392L,all(is.finite(as.matrix(ko[c(DEG_KOS,PROD_KOS)]))))
reconstructed <- cbind(Degradation=rowSums(ko[DEG_KOS]),Production=rowSums(ko[PROD_KOS]))
reconstructed <- cbind(reconstructed,Production_Degradation=log2((reconstructed[,2]+1e-8)/(reconstructed[,1]+1e-8)))
scorecols <- colnames(reconstructed)
stopifnot(max(abs(reconstructed-as.matrix(sc[scorecols])))<1e-12)
e_stats <- rbindlist(lapply(scorecols,function(s){
  x<-sc[[s]][sc$Group=='R'];y<-sc[[s]][sc$Group=='NR'];w<-wilcox.test(x,y,exact=FALSE,correct=TRUE)
  data.table(ID=s,n_R=length(x),n_NR=length(y),wilcox_p=w$p.value,effect=2*as.numeric(w$statistic)/(length(x)*length(y))-1)
}))
e_stats[,q:=p.adjust(wilcox_p,'BH')];oldstat<-rd(score_stats_path)
stopifnot(max(abs(e_stats$wilcox_p-oldstat$wilcox_p))<1e-12,max(abs(e_stats$q-oldstat$q))<1e-12,max(abs(e_stats$effect-oldstat$effect))<1e-12)
wr(sc,'Fig3e_per_Run_scores');wr(e_stats,'Fig3e_statistics')
ns <- list(); inspection <- list()
for(co in names(NSCLC_DISC_DIRS)){
  spf<-file.path(used,co,paste0(co,'_species_relative_abundance_matrix.csv'))
  metaf<-file.path(used,co,paste0(co,'_patient_metadata_WGS_NSCLC_R_NR.csv'))
  source_files<-c(source_files,spf,metaf)
  sp<-rd(spf);meta<-rd(metaf);d<-sc[sc$Cohort==co,]
  stopifnot(!anyDuplicated(sp$Run_ID),setequal(sp$Run_ID,d$Run_ID),!anyDuplicated(meta[['Run ID']]))
  ix<-match(d$Run_ID,sp$Run_ID);mi<-match(d$Run_ID,meta[['Run ID']]);stopifnot(!anyNA(ix),!anyNA(mi),all(d$Group==meta[['Analysis Group']][mi]))
  m<-as.matrix(sp[ix,-1,drop=FALSE]); storage.mode(m)<-'double'
  stopifnot(all(is.finite(m)),all(m>=0),max(abs(rowSums(m)-100))<.001)
  colnames(m)<-sub('.*\\|s__','',colnames(m));stopifnot(!anyDuplicated(colnames(m)))
  rownames(m)<-d$Run_ID
  ns[[co]]<-list(cohort=co,species=m,score=d$Degradation,production=d$Production,group=d$Group,ids=d$Sample_name,run_ids=d$Run_ID)
  inspection[[co]]<-data.table(disease='NSCLC',cohort=co,n=nrow(m),features=ncol(m),missing=sum(is.na(m)),min=min(m),max=max(m),min_sum=min(rowSums(m)),max_sum=max(rowSums(m)))
}
wr(rbindlist(inspection),'input_species_inspection')
overlap<-as.data.table(sc)[,.(n_runs=.N,n_cohorts=uniqueN(Cohort),n_groups=uniqueN(Group)),by=Sample_name][n_cohorts>1]
wr(overlap,'NSCLC_cross_project_repeated_names');say('Retained overlap: %d Sample names; %d group conflicts',nrow(overlap),sum(overlap$n_groups>1))
stopifnot(nrow(overlap)==283L,sum(overlap$n_groups>1)==45L)
say('[2] Recompute Figure i species associations with both new scores')
ns_unadj<-lapply(ns,association_cohort,adjusted=FALSE)
ns_prod<-lapply(ns,function(d){d$score<-d$production;association_cohort(d,adjusted=FALSE)})
im_deg<-meta_disease(ns_unadj,'NSCLC');im_prod<-meta_disease(ns_prod,'NSCLC')
for(s in c('Degradation','Production')){
  a<-if(s=='Degradation')ns_unadj else ns_prod;m<-if(s=='Degradation')im_deg else im_prod
  wr(rbindlist(lapply(a,`[[`,'table')),paste0('Fig3i_',s,'_percohort'));wr(m,paste0('Fig3i_',s,'_meta'))
}
da<-as.data.table(rd(da_path));stopifnot(!anyDuplicated(da$species))
i_data<-merge(im_deg[,.(species,Degradation_rho=pooled_rho,Degradation_q=q_value_BH)],im_prod[,.(species,Production_rho=pooled_rho,Production_q=q_value_BH)],by='species')
i_data<-merge(i_data,da[,.(species,DA_direction=dir,DA_q=qval)],by='species')
i_data<-i_data[DA_q<.05 & DA_direction %in% c('higher in R','higher in NR')]
i_data[,enriched:=ifelse(DA_direction=='higher in R','R','NR')]
stopifnot(nrow(i_data)>0,all(is.finite(i_data$Degradation_rho)),all(is.finite(i_data$Production_rho)))
i_tests<-rbindlist(lapply(c('Degradation','Production'),function(s){
  x<-i_data[[paste0(s,'_rho')]];w<-wilcox.test(x[i_data$enriched=='R'],x[i_data$enriched=='NR'],exact=FALSE,correct=TRUE)
  data.table(score=s,n_R_species=sum(i_data$enriched=='R'),n_NR_species=sum(i_data$enriched=='NR'),p=w$p.value)
}));i_tests[,q_two_comparisons:=p.adjust(p,'BH')]
label_species<-c(head(i_data[enriched=='R'][order(-Degradation_rho)]$species,6),head(i_data[enriched=='NR'][order(Degradation_rho)]$species,6))
i_data[,label:=species %in% label_species]
wr(i_data,'Fig3i_plot_data');wr(i_tests,'Fig3i_descriptive_species_comparisons')
say('[3] Load original CRC BioSample preprocessing for Figure j')
ko_cache<-rd(KO_CACHE_FILE);stopifnot(!anyDuplicated(ko_cache$Run.ID),nrow(ko_cache)==1647L)
crc_raw<-Map(load_crc_cohort,names(CRC_DIRS),unname(CRC_DIRS),MoreArgs=list(ko_cache=ko_cache));names(crc_raw)<-names(CRC_DIRS)
source_files<-unique(c(source_files,unlist(lapply(crc_raw,`[[`,'input_files'))))
crc_primary<-lapply(crc_raw,aggregate_crc_cohort,layout='paired',fun='mean')
wr(rbindlist(lapply(crc_primary,`[[`,'unit_qc')),'Fig3j_CRC_BioSample_roster')
crc_assoc<-lapply(crc_primary,association_cohort,adjusted=TRUE)
ns_adj<-lapply(ns,association_cohort,adjusted=TRUE)
crc_meta<-meta_disease(crc_assoc,'CRC');ns_meta<-meta_disease(ns_adj,'NSCLC')
wr(rbindlist(lapply(crc_assoc,`[[`,'table')),'Fig3j_CRC_percohort');wr(rbindlist(lapply(ns_adj,`[[`,'table')),'Fig3j_NSCLC_percohort')
wr(crc_meta,'Fig3j_CRC_meta');wr(ns_meta,'Fig3j_NSCLC_meta')
j_data<-merge(crc_meta[,.(species,CRC_rho=pooled_rho,CRC_q=q_value_BH)],ns_meta[,.(species,NSCLC_rho=pooled_rho,NSCLC_q=q_value_BH)],by='species')
j_data[,significance:=ifelse(CRC_q<.05 & NSCLC_q<.05,'Both',ifelse(CRC_q<.05,'CRC only',ifelse(NSCLC_q<.05,'NSCLC only','Neither')))]
j_stat<-concordance_summary(crc_meta,ns_meta,'Primary: original CRC BioSample handling; new NSCLC five-KO score')
set.seed(42);per<-permute_concordance(crc_assoc,ns_adj,j_data$species,j_stat$spearman_rho,reps=N_PERM)
j_stat$permutation_reps<-N_PERM;j_stat$conditional_permutation_p<-per$p
wr(j_data,'Fig3j_plot_data');wr(j_stat,'Fig3j_statistics');wr(data.table(replicate=seq_along(per$null),rho=per$null),'Fig3j_permutation_null')
say('[4] Prespecified score/membership and cohort sensitivities')
sens<-list(j_stat[,1:7])
scenario<-function(cd,nd,adjusted,label){ca<-lapply(cd,association_cohort,adjusted=adjusted);na<-lapply(nd,association_cohort,adjusted=adjusted);concordance_summary(meta_disease(ca,'CRC'),meta_disease(na,'NSCLC'),label)}
sens[[2]]<-scenario(crc_primary,ns,FALSE,'Unadjusted Spearman')
sens[[3]]<-scenario(lapply(crc_raw,aggregate_crc_cohort,layout='all',fun='mean'),ns,TRUE,'CRC all-layout mean')
sens[[4]]<-scenario(lapply(crc_raw,aggregate_crc_cohort,layout='paired',fun='median'),ns,TRUE,'CRC paired-layout median')
for(co in names(ns_adj))sens[[paste0('dropNSCLC',co)]]<-concordance_summary(crc_meta,meta_disease(ns_adj[names(ns_adj)!=co],'NSCLC'),paste('Omit NSCLC',co))
for(co in names(crc_assoc))sens[[paste0('dropCRC',co)]]<-concordance_summary(meta_disease(crc_assoc[names(crc_assoc)!=co],'CRC'),ns_meta,paste('Omit CRC',co))
sens[['allcohorts']]<-concordance_summary(meta_disease(crc_assoc,'CRC',4L),meta_disease(ns_adj,'NSCLC',3L),'Species eligible in all4 CRC and all3 NSCLC datasets')
wr(rbindlist(sens,fill=TRUE),'Fig3j_sensitivity')
# i leave-one-project-out estimates: do not silently deduplicate unknown identities.
i_loo<-rbindlist(lapply(names(ns),function(co){
  de<-meta_disease(ns_unadj[names(ns)!=co],'NSCLC');pr<-meta_disease(ns_prod[names(ns)!=co],'NSCLC')
  m<-merge(de[,.(species,Degradation_rho=pooled_rho,Degradation_q=q_value_BH)],pr[,.(species,Production_rho=pooled_rho,Production_q=q_value_BH)],by='species')
  m[,omitted_cohort:=co];m
}));wr(i_loo,'Fig3i_leave_one_cohort_out')
# k intentionally preserves original Run-level CRC core scope, distinct from j BioSamples.
crc_run<-lapply(crc_raw,function(d){ids<-intersect(rownames(d$bact),names(d$score));ix<-match(ids,d$meta$Run.ID)
  list(cohort=d$cohort,species=d$bact[ids,,drop=FALSE],score=as.numeric(d$score[ids]),group=d$meta$Group[ix])})
k_crc<-lapply(crc_run,association_cohort,adjusted=FALSE)
kc<-rbindlist(lapply(k_crc,`[[`,'table'));kn<-rbindlist(lapply(ns_unadj,`[[`,'table'))
# cor.test exact=FALSE gives the same t approximation used by association_cohort.
kc[,selected:=rho>.3 & p_adj_BH<.05];kn[,selected:=rho>.3 & p_value<.05];kn[,selected_BH:=rho>.3 & p_adj_BH<.05]
wr(kc,'Fig3k_CRC_correlations');wr(kn,'Fig3k_NSCLC_correlations')
cs<-kc[selected==TRUE,.(n_selected=uniqueN(cohort)),by=species][n_selected==4]$species
ns_core<-kn[selected==TRUE,.(n_selected=uniqueN(cohort)),by=species][n_selected==3]$species
ns_bh<-kn[selected_BH==TRUE,.(n_selected=uniqueN(cohort)),by=species][n_selected==3]$species
univ<-intersect(unique(kc$species),union(im_deg$species,im_prod$species))
# Exact historical overlap-background construction: CRC per-cohort tested union
# intersect NSCLC meta-eligible species (union of the two score test families).
# This is a descriptive set statistic, not a test of independent microbial features.
stopifnot(all(cs %in% univ),all(ns_core %in% univ))
sets<-function(nset,scenario){
  species<-sort(union(cs,nset));membership<-data.table(species=species,in_CRC=species %in% cs,in_NSCLC=species %in% nset)
  membership[,shared:=in_CRC & in_NSCLC]
  a<-length(intersect(cs,univ));b<-length(intersect(nset,univ));ov<-length(intersect(cs,nset));N<-length(univ)
  data<-data.table(scenario=scenario,n_CRC=length(cs),n_NSCLC=length(nset),overlap=ov,universe=N,
    CRC_in_universe=a,NSCLC_in_universe=b,expected=a*b/N,
    fold_enrichment=if(a*b>0)ov/(a*b/N) else NA_real_,hypergeom_p=if(a*b>0)phyper(ov-1,a,N-a,b,lower.tail=FALSE) else 1,
    Leligens_CRC='Lachnospira_eligens' %in% cs,Leligens_NSCLC='Lachnospira_eligens' %in% nset)
  list(membership=membership,stats=data)
}
k<-sets(ns_core,'Original selection: CRC BH; NSCLC nominal');kb<-sets(ns_bh,'Uniform per-cohort BH sensitivity')
wr(k$membership,'Fig3k_plot_data');wr(kb$membership,'Fig3k_uniform_BH_membership');wr(rbind(k$stats,kb$stats),'Fig3k_statistics');wr(data.table(species=univ),'Fig3k_overlap_universe')
say('Core: CRC %d, NSCLC %d, shared %d; uniform-BH NSCLC %d',length(cs),length(ns_core),length(intersect(cs,ns_core)),length(ns_bh))
say('[5] Validate all DL fits against metafor and retain input hashes')
checks<-rbindlist(list(validate_meta_api(im_deg,ns_unadj,'i_Degradation'),validate_meta_api(im_prod,ns_prod,'i_Production'),validate_meta_api(crc_meta,crc_assoc,'j_CRC'),validate_meta_api(ns_meta,ns_adj,'j_NSCLC')))
wr(checks,'DL_vs_metafor_checks')
source_files<-unique(c(source_files,script,file.path(dirname(script),'00_verified_helpers.R'),file.path(OUT_DIR,'ANALYSIS_PLAN_FROZEN.md')))
source_hash<-data.table(path=vapply(source_files,rel,''),sha256=vapply(source_files,hash,''),bytes=file.info(source_files)$size)
fwrite(source_hash,file.path(OUT_DIR,'qa/input_sha256.csv'))
stopifnot(identical(protected_sha,vapply(protected,hash,'')))
fwrite(data.table(path=vapply(protected,rel,''),sha256=protected_sha,unchanged=TRUE),file.path(OUT_DIR,'qa/protected_files.csv'))
writeLines(capture.output(sessionInfo()),file.path(OUT_DIR,'qa/sessionInfo.txt'))
saveRDS(list(scores=sc,e_stats=e_stats,i_data=i_data,i_tests=i_tests,j_data=j_data,j_stat=j_stat,k_data=k$membership,k_stats=rbind(k$stats,kb$stats)),file.path(OUT_DIR,'tables/plot_objects.rds'))
saveRDS(list(ns=ns,crc_primary=crc_primary,crc_run=crc_run),file.path(OUT_DIR,'qa/verification_inputs.rds'))
say('DONE: e same as f; i %d taxa; j rho %.6f; k overlap %d',nrow(i_data),j_stat$spearman_rho,k$stats$overlap)
