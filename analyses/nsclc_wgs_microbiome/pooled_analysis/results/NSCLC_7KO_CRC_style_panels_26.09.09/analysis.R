#!/usr/bin/env Rscript
# NSCLC analogue of CRC Fig2h/i/j and Sup11--14. New outputs only.
# Frozen design (before inspecting new associations):
# 824 Run-level observations, 3 discovery cohorts; R and NR, NSCLC only.
# Percent relative abundances, zeros retained, no deduplication/renormalization.
# Seven KOs: degradation K00301/2/3/5/6; production K00315/K08688.
# Forest: P(NR>R)-P(NR<R), 5000 group-stratified bootstrap draws, seed42
# reset per KO/cohort, type7 percentile 95% CI; Wilcoxon two-sided,
# exact=FALSE, correct=TRUE; BH across estimable cohorts WITHIN each KO.
# All-zero contrasts: NA effect/CI/p/q, explicitly non-estimable.
# Current CRC PPT display: filled fixed-size forest points; color encodes
# BH significance and direction; pooled lollipop area encodes prevalence.
# Pooled KO BH is a different family: 7 KOs; log2FC PC1e-8 percent.
# Species: prevalence>=10% within each analysis scope; Spearman exact=FALSE.
# KO-species BH across all estimable 7KO x eligible-species pairs per scope.
# Pooled candidates: any pair |rho|>.3 & BHq<.05, top20 by max |rho|;
# ties species identifier alphabetical; pathway role from qualifying KO pairs.
# Associated-species DA: same pooled top20 tested in each scope, Wilcoxon,
# BH within estimable selected species per scope, log2FC PC1e-6 percent.
# DA bars show q<.05 selected species, NO additional fold-change cutoff.
# Also retain global DA BH across all prevalence-eligible species per scope.
# Venn pathway membership: score=sum(current 5deg or 2prod KOs), Spearman
# rho>.3 & BHq<.05, BH across eligible species separately per score/cohort.
# Venn NR/R membership: global DA BHq<.05 & log2FC >1 / <-1 respectively.
# Scatter selection: up to 6 pooled top20 candidates with global DA BHq<.05,
# ranked by max |rho| (no R/NR or production/degradation quotas); show ALL
# 7 KO pairs, including nonsignificant ones; q retained from full pair family.
# Selection and inference share these data: exploratory, not validation.
# Pooled results unadjusted for cohort/covariates; repeated cross-project
# Sample_name values retained, so Runs must not be called unique patients.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest);library(stringi)})
options(stringsAsFactors=FALSE); setDTthreads(2)
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
OUT <- dirname(script); NS <- dirname(dirname(dirname(OUT)))
used <- list.dirs(NS,recursive=FALSE,full.names=TRUE)
used <- used[stri_trans_nfc(basename(used))=='사용데이터_모음'];stopifnot(length(used)==1)
IN <- file.path(used,'NSCLC_Species_DA_PRISM_26.09.09')
CSV <- file.path(OUT,'csv');dir.create(CSV,showWarnings=FALSE)
savecsv <- function(d,n) fwrite(d,file.path(CSV,paste0(n,'.csv')),na='',bom=TRUE)
# Resolve the already verified exports by content: the user renamed files.
expected_sha <- c('18d011277f00c4d3d5f7a00599189dac440b84b35add4820795b5bfc71135c98',
 '7fec84c000d39477ac6a5d4f1decc8f24bbbfe042e1a0e7719aa7f54b45249f3',
 'a0c46ff4b35ac44710f9c3ab3797c500d2ecc0a298d9ed215ac4b3fbf1bec86c',
 '1e2f8808ceef82f4323e22fa6e8a9a14d2939e7972abcf2af05169990dd01fb2')
available <- list.files(IN,pattern='[.]csv$',full.names=TRUE)
available_sha <- vapply(available,function(p)digest(file=p,algo='sha256',serialize=FALSE),'')
files <- vapply(expected_sha,function(h){p<-available[available_sha==h];stopifnot(length(p)==1);p},'')
before <- vapply(files,function(p)digest(file=p,algo='sha256',serialize=FALSE),'')
meta <- fread(files[1]);sp <- fread(files[2]);ko <- fread(files[3]);sc <- fread(files[4])
CO <- c('PRJNA751792','PRJNA1023797','PRJEB22863')
KOS <- c('K00301','K00302','K00303','K00305','K00306','K00315','K08688')
genes <- c('Sarcosine oxidase','soxA','soxB','soxG','PIPOX','DMGDH','Creatinase')
roles <- c(rep('Degradation',5),rep('Production',2))
stopifnot(nrow(meta)==824,!anyDuplicated(meta$Run_ID),sum(meta$Group=='R')==432,sum(meta$Group=='NR')==392,
 identical(as.integer(table(factor(meta$Cohort,levels=CO))),c(338L,421L,65L)))
for(d in list(sp,ko,sc))stopifnot(identical(d$Run_ID,meta$Run_ID),identical(d$Group,meta$Group),identical(d$Cohort,meta$Cohort))
S <- as.matrix(sp[,-(1:4)]);K <- as.matrix(ko[,..KOS]);stopifnot(all(is.finite(S)),all(is.finite(K)),all(S>=0),all(K>=0),max(abs(rowSums(S)-100))<.001)
taxa <- colnames(S);short <- sub('.*\\|s__','',taxa); stopifnot(!anyDuplicated(short))
names(short) <- taxa
scores <- cbind(Degradation=rowSums(K[,1:5]),Production=rowSums(K[,6:7]))
stopifnot(max(abs(scores[,1]-sc$Degradation_5KO_percent))<1e-12,max(abs(scores[,2]-sc$Production_2KO_percent))<1e-12)
overlap <- meta[,.(n_cohorts=uniqueN(Cohort),n_groups=uniqueN(Group)),by=Sample_name][n_cohorts>1]
stopifnot(nrow(overlap)==283,sum(overlap$n_groups>1)==45)
scopes <- c(list(Pooled=seq_len(nrow(meta))),setNames(lapply(CO,function(co)which(meta$Cohort==co)),CO))
bh <- function(p){out<-rep(NA_real_,length(p));good<-is.finite(p);out[good]<-p.adjust(p[good],'BH');out}
rbc <- function(x,y){nx<-length(x);ny<-length(y);2*(sum(rank(c(x,y),ties.method='average')[seq_len(nx)])-nx*(nx+1)/2)/(nx*ny)-1}
color <- function(q,e)ifelse(is.na(q),'#9D9D9D',ifelse(q<.05,ifelse(e<0,'#1B9E8F','#C43C3C'),'#9D9D9D'))
cat('Computing 7 KO x 3 cohort contrasts ...\n');flush.console()
forest <- rbindlist(lapply(seq_along(KOS),function(j) rbindlist(lapply(CO,function(co){
 ix<-scopes[[co]];r<-K[ix[meta$Group[ix]=='R'],j];nr<-K[ix[meta$Group[ix]=='NR'],j];z<-all(c(r,nr)==0);estimable<-length(unique(c(r,nr)))>1
 ci<-c(NA_real_,NA_real_);p<-e<-NA_real_
 if(estimable){e<-rbc(nr,r);p<-wilcox.test(nr,r,exact=FALSE,correct=TRUE)$p.value;set.seed(42)
  bs<-replicate(5000,rbc(nr[sample.int(length(nr),length(nr),replace=TRUE)],r[sample.int(length(r),length(r),replace=TRUE)]))
  ci<-quantile(bs,c(.025,.975),type=7,names=FALSE)}
 data.table(KO=KOS[j],Gene=genes[j],Role=roles[j],Cohort=co,n_R=length(r),n_NR=length(nr),
  nonzero_R=sum(r>0),nonzero_NR=sum(nr>0),prevalence_cohort=mean(c(r,nr)>0),prevalence_pooled=mean(K[,j]>0),
  mean_R_percent=mean(r),mean_NR_percent=mean(nr),effect_NR_minus_R=e,CI95_low=ci[1],CI95_high=ci[2],
  p_value=p,estimable=estimable,status=if(z)'All zero' else if(!estimable)'Constant' else 'Estimable')
}))))
forest[,BH_q:=bh(p_value),by=KO];forest[,BH_family_n:=sum(estimable),by=KO]
forest[,`:=`(CI_error_minus=effect_NR_minus_R-CI95_low,CI_error_plus=CI95_high-effect_NR_minus_R,
 direction=ifelse(!estimable,'Not estimable',ifelse(effect_NR_minus_R<0,'R higher',ifelse(effect_NR_minus_R>0,'NR higher','Equal'))),
 point_size=3.2,fill_status=ifelse(!estimable,'None','Filled'),color_hex=color(BH_q,effect_NR_minus_R))]
savecsv(forest,'KO_7genes_3cohorts_forest')
cat('Computing species correlations and differential abundance ...\n');flush.console()
KOcor<-list();Scorecor<-list();DA<-list()
for(a in names(scopes)){
 ix<-scopes[[a]];sm<-S[ix,,drop=FALSE];km<-K[ix,,drop=FALSE];prev<-colMeans(sm>0);elig<-which(prev>=.10)
 cor_one<-function(x,y){if(length(unique(x))<2||length(unique(y))<2)return(c(NA_real_,NA_real_));w<-suppressWarnings(cor.test(x,y,method='spearman',exact=FALSE));c(unname(w$estimate),w$p.value)}
 kd<-rbindlist(lapply(seq_along(KOS),function(j)rbindlist(lapply(elig,function(k){v<-cor_one(km[,j],sm[,k]);data.table(Analysis=a,KO=KOS[j],Gene=genes[j],KO_role=roles[j],Species_full=taxa[k],Species=short[k],n=length(ix),prevalence=prev[k],rho=v[1],p_value=v[2])}))))
 kd[,BH_q:=bh(p_value)];kd[,BH_family_n:=sum(is.finite(p_value))];kd[,associated:=is.finite(BH_q)&BH_q<.05&abs(rho)>.3];KOcor[[a]]<-kd
 Scorecor[[a]]<-rbindlist(lapply(colnames(scores),function(s){
  d<-rbindlist(lapply(elig,function(k){v<-cor_one(scores[ix,s],sm[,k]);data.table(Analysis=a,Score=s,Species_full=taxa[k],Species=short[k],n=length(ix),prevalence=prev[k],rho=v[1],p_value=v[2])}))
  d[,BH_q:=bh(p_value)];d[,BH_family_n:=sum(is.finite(p_value))];d[,associated_positive:=is.finite(BH_q)&BH_q<.05&rho>.3];d
 }))
 DA[[a]]<-rbindlist(lapply(elig,function(k){r<-sm[meta$Group[ix]=='R',k];nr<-sm[meta$Group[ix]=='NR',k]
  data.table(Analysis=a,Species_full=taxa[k],Species=short[k],n_R=length(r),n_NR=length(nr),prevalence=prev[k],
   mean_R_percent=mean(r),mean_NR_percent=mean(nr),log2FC_NR_vs_R=log2((mean(nr)+1e-6)/(mean(r)+1e-6)),p_value=wilcox.test(nr,r,exact=FALSE,correct=TRUE)$p.value)
 }))
 DA[[a]][,BH_q_global:=bh(p_value)];DA[[a]][,BH_family_n:=sum(is.finite(p_value))]
 cat(a,':',length(elig),'eligible species,',sum(kd$associated),'KO correlations passing q and rho\n');flush.console()
}
ko_cor<-rbindlist(KOcor);score_cor<-rbindlist(Scorecor);da<-rbindlist(DA)
savecsv(ko_cor,'Species_7KO_Spearman_all_tests');savecsv(score_cor,'Species_pathway_score_Spearman_all_tests');savecsv(da,'Species_DA_all_eligible')
candidates<-KOcor$Pooled[associated==TRUE,.(max_abs_rho=max(abs(rho)),associated_KOs=paste(KO,collapse=';'),
 Role=if(uniqueN(KO_role)==2)'Both' else unique(KO_role)),by=.(Species_full,Species)]
setorder(candidates,-max_abs_rho,Species_full);candidates[,Correlation_rank:=seq_len(.N)];candidates[,In_top20:=Correlation_rank<=20]
top<-candidates[In_top20==TRUE];stopifnot(nrow(top)>0)
associated_DA<-rbindlist(lapply(names(scopes),function(a){
 ix<-scopes[[a]];ridx<-ix[meta$Group[ix]=='R'];nridx<-ix[meta$Group[ix]=='NR']
 d<-rbindlist(lapply(seq_len(nrow(top)),function(i){r<-S[ridx,top$Species_full[i]];nr<-S[nridx,top$Species_full[i]];est<-length(unique(c(r,nr)))>1
  data.table(Analysis=a,Species_full=top$Species_full[i],Species=top$Species[i],Role=top$Role[i],
   max_abs_rho=top$max_abs_rho[i],Correlation_rank=top$Correlation_rank[i],n_R=length(r),n_NR=length(nr),
   mean_R_percent=mean(r),mean_NR_percent=mean(nr),log2FC_NR_vs_R=log2((mean(nr)+1e-6)/(mean(r)+1e-6)),
   p_value=if(est)wilcox.test(nr,r,exact=FALSE,correct=TRUE)$p.value else NA_real_,estimable=est)
 }))
 d[,BH_q_selected:=bh(p_value)];d[,BH_family_n:=sum(estimable)];d[,Plot_bar:=is.finite(BH_q_selected)&BH_q_selected<.05]
 d[,`:=`(Direction=ifelse(log2FC_NR_vs_R<0,'R-enriched','NR-enriched'),bar_color_hex=ifelse(log2FC_NR_vs_R<0,'#1B9E8F','#C43C3C'))]
 d
}))
savecsv(associated_DA,'Sarcosine_associated_species_DA_pooled_and_3cohorts')
candidates<-merge(candidates,DA$Pooled[,.(Species_full,DA_BH_q_global=BH_q_global,DA_log2FC_NR_vs_R=log2FC_NR_vs_R)],by='Species_full',all.x=TRUE,sort=FALSE)
setorder(candidates,Correlation_rank)
selected<-head(candidates[In_top20==TRUE & DA_BH_q_global<.05],6)
candidates[,Selected_for_scatter:=Species_full %in% selected$Species_full]
savecsv(candidates,'Species_selection_candidates');savecsv(selected,'Scatter_selected_species')
scatter_stats<-ko_cor[Species_full %in% selected$Species_full]
scatter_stats[,Scatter_column:=match(Species_full,selected$Species_full)];savecsv(scatter_stats,'Scatter_7KO_statistics_pooled_and_3cohorts')
raw<-cbind(meta[,.(Cohort,Run_ID,Group,Sample_name)],as.data.table(K),
 data.table(Degradation_5KO_percent=scores[,1],Production_2KO_percent=scores[,2],log2_Production_Degradation_pc1e_8=log2((scores[,2]+1e-8)/(scores[,1]+1e-8))),
 as.data.table(S[,selected$Species_full,drop=FALSE]))
savecsv(raw,'Scatter_7KO_and_selected_species_per_Run')
# All tested species retained in a wide Boolean membership matrix; absence
# from one eligible universe is FALSE membership, NOT an estimated null effect.
categories<-c('Production-associated','NR-enriched','Degradation-associated','R-enriched')
venn<-rbindlist(lapply(categories,function(cat){
 universe<-sort(unique(unlist(lapply(CO,function(co)DA[[co]]$Species_full))))
 d<-data.table(Category=cat,Species_full=universe,Species=short[universe])
 for(co in CO){
  ids<-if(cat %in% c('Production-associated','Degradation-associated')) Scorecor[[co]][Score==sub('-associated','',cat)&associated_positive==TRUE]$Species_full else
   DA[[co]][BH_q_global<.05 & if(cat=='NR-enriched')log2FC_NR_vs_R>1 else log2FC_NR_vs_R< -1]$Species_full
  d[,(co):=Species_full %in% ids]
 }
 d[,n_cohorts:=rowSums(.SD),.SDcols=CO];d
}))
savecsv(venn,'Venn_3cohorts_membership')
venn_counts<-rbindlist(lapply(categories,function(cat){
 d<-venn[Category==cat];r<-CJ(A=c(FALSE,TRUE),B=c(FALSE,TRUE),C=c(FALSE,TRUE))[A|B|C]
 r[,N:=vapply(seq_len(.N),function(i)sum(d[[CO[1]]]==A[i]&d[[CO[2]]]==B[i]&d[[CO[3]]]==C[i]),integer(1))]
 setnames(r,c('A','B','C'),CO);r[,Category:=cat];setcolorder(r,c('Category',CO,'N'));r
}));savecsv(venn_counts,'Venn_3cohorts_region_counts')
# Recalculate pooled point estimates and verify frozen seven-KO statistics.
frozen_path<-file.path(dirname(OUT),'CRC_Method_7KO_NSCLC_26.09.07','tables','NSCLC_seven_KO_pooled_statistics.csv')
frozen<-fread(frozen_path);pooled<-rbindlist(lapply(seq_along(KOS),function(j){r<-K[meta$Group=='R',j];nr<-K[meta$Group=='NR',j]
 data.table(KO=KOS[j],Gene=genes[j],Role=roles[j],n_R=length(r),n_NR=length(nr),mean_R_percent=mean(r),mean_NR_percent=mean(nr),
  prevalence=mean(K[,j]>0),log2FC_NR_vs_R=log2((mean(nr)+1e-8)/(mean(r)+1e-8)),effect_NR_minus_R=rbc(nr,r),p_value=wilcox.test(nr,r,exact=FALSE,correct=TRUE)$p.value)
}))
pooled[,BH_q:=bh(p_value)];frozen<-frozen[match(pooled$KO,ID)]
stopifnot(max(abs(pooled$p_value-frozen$wilcox_p))<1e-12,max(abs(pooled$BH_q-frozen$q))<1e-12,
 max(abs(pooled$log2FC_NR_vs_R-frozen$log2FC_NR_R))<1e-12,max(abs(pooled$effect_NR_minus_R+frozen$effect))<1e-12)
pooled[,`:=`(prevalence_percent=prevalence*100,BH_family_n=7L,point_area_prevalence=prevalence,
 fill_status='Filled',color_hex=color(BH_q,log2FC_NR_vs_R))]
savecsv(pooled,'KO_7genes_pooled_lollipop')
summary<-data.table(Analysis=names(scopes),n=vapply(scopes,length,integer(1)),
 n_R=vapply(scopes,function(ix)sum(meta$Group[ix]=='R'),integer(1)),n_NR=vapply(scopes,function(ix)sum(meta$Group[ix]=='NR'),integer(1)),
 n_species_tested=vapply(DA,nrow,integer(1)),n_KO_species_tests=vapply(KOcor,function(d)sum(is.finite(d$p_value)),integer(1)),
 n_significant_associated_DA_bars=vapply(names(scopes),function(a)sum(associated_DA[Analysis==a]$Plot_bar),integer(1)))
savecsv(summary,'Analysis_sample_and_test_counts')
savecsv(data.table(Source=files,SHA256=before),'Input_source_SHA256')
stopifnot(identical(before,vapply(files,function(p)digest(file=p,algo='sha256',serialize=FALSE),'')))
saveRDS(list(meta=meta,S=S,K=K,short=short,scopes=scopes,CO=CO,KOS=KOS,genes=genes,roles=roles,
 forest=forest,pooled=pooled,associated_DA=associated_DA,selected=selected,ko_cor=ko_cor,venn=venn,summary=summary),file.path(OUT,'plot_data.rds'))
print(summary);cat('\nScatter species:\n');print(selected[,.(Species,Role,max_abs_rho,DA_BH_q_global)])
cat('\nVenn membership:\n');print(venn[,lapply(.SD,sum),by=Category,.SDcols=CO]);print(venn[n_cohorts==3,.N,by=Category])
cat('\nAll-zero contrasts:\n');print(forest[estimable==FALSE,.(KO,Cohort,status)])
