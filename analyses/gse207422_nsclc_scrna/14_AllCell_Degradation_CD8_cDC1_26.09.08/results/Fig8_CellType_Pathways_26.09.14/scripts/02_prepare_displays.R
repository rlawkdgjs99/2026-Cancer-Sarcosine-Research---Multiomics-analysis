#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])),"common.R"))
stopifnot(all(fread(file.path(out,'results/tables/06_verification_checks.csv'))$passed))
checks<-list();ck<-function(k,ok,d=''){checks[[length(checks)+1L]]<<-data.table(check=k,pass=isTRUE(ok),detail=as.character(d));if(!isTRUE(ok))stop(k,': ',d)}
plan<-file.path(out,'ANALYSIS_DISPLAY_PLAN.md');ck('display plan unchanged',sha(plan)==readLines(file.path(out,'logs/display_plan_sha256.txt')))
a<-readRDS(file.path(origin,'resources/allcell_prepared.rds'));prep<-readRDS(file.path(origin,'resources/outcome_prepared.rds'))
g15<-fread(file.path(src,'15_AllCell_Degradation_Pathway_Overview_26.09.14/results/tables/01_frozen_patient_groups.csv'))
ck('same groups and D as all-cell15',identical(a$metadata$Patient,g15$Patient)&&all(a$metadata$group==g15$group)&&max(abs(a$metadata$degradation_mean_z-g15$degradation_mean_z))<1e-14)
oldcd8<-readRDS(file.path(src,'02_lineage_reannotation/intermediate/10_lineage_pseudobulk_counts.rds'))[['CD8 T cell']]
oldc<-readRDS(file.path(src,'12_cDC1_CD8_Defined_Groups_Exploratory_26.09.08/resources/patient_pseudobulk.rds'))$counts
ck('CD8 independent earlier raw aggregate',all(a$counts$CD8==as.matrix(oldcd8[rownames(a$counts$CD8),a$metadata$Sample])))
ck('cDC1 independent earlier raw aggregate',all(a$counts$cDC1[,match(colnames(oldc),a$metadata$Patient),drop=FALSE]==oldc[rownames(a$counts$cDC1),,drop=FALSE]))
for(ct in names(prep)){
 pp<-prep[[ct]];yy<-calcNormFactors(pp$y,method='TMM');ck(paste(ct,'TMM reproduced'),max(abs(yy$samples$norm.factors-pp$y$samples$norm.factors))<1e-12)
}
m<-as.data.table(readRDS(member_path));sets<-lapply(split(m$gene_symbol,m$pathway),unique)
exact<-fread(file.path(src,'04_Fig7_Immune_Programs_26.09.06/resources/Fig6_exact_MSigDB_2026.1.Hs_membership.csv'))
ck('four pathway identities',all(vapply(focus$pathway[1:4],function(id)setequal(sets[[id]],exact[gs_name==id,gene_symbol]),logical(1))))
g<-fread(file.path(origin,'results/tables/02_GSEA_all_scopes.csv'))
ca<-fread(file.path(origin,'results/tables/03_CAMERA_primary.csv'));ss<-fread(file.path(origin,'results/tables/04_exact_scores_HC3.csv'))
v<-fread(file.path(origin,'results/tables/04_patient_score_values.csv'))
ck('full GSEA result count',nrow(g)==23836&&g[scope=='post_group',.N]==5959)
elig<-fread(file.path(origin,'results/tables/01_model_eligibility.csv'));wt(elig,'model_eligibility.csv')
selected<-list();audit<-list();counts<-list();es<-list()
for(ct in c('CD8','cDC1')){
 b<-readRDS(file.path(origin,'resources/models',paste0(ct,'_post_group.rds')));rank<-b$rank
 dt<-copy(b$patients);wt(dt[,.(Patient,Sample,treatment,histology,group,degradation_mean_z,CD8_cells,cDC1_cells)],paste0('patients_',ct,'.csv'))
 for(sc in unique(ss$scope)){
  bb<-readRDS(file.path(origin,'resources/models',paste0(ct,'_',sc,'.rds')))
  wt(data.table(Patient=rownames(bb$design),as.data.table(bb$design)),paste0('design_',ct,'_',sc,'.csv'))
 }
 measured<-lapply(sets,intersect,y=names(rank))
 for(dr in c('High','Low')){
  d<-copy(g[cell_type==ct&scope=='post_group'&collection=='GO:BP'&q_global<.05&is.finite(q_global)])
  d<-if(dr=='High')d[NES>0] else d[NES<0]
  d[,absNES:=abs(NES)];setorderv(d,c('q_global','absNES','pathway'),c(1,-1,1));chosen<-character();orders<-integer();reasons<-character();maxjac<-numeric()
  for(i in seq_len(nrow(d))){
   id<-d$pathway[i];jacs<-if(length(chosen))vapply(chosen,function(x)length(intersect(measured[[id]],measured[[x]]))/length(union(measured[[id]],measured[[x]])),numeric(1)) else 0
   jj<-max(jacs);decision<-if(length(chosen)>=8)'CAP_REACHED' else if(jj>=.5)'GENE_OVERLAP' else 'SELECTED'
   if(decision=='SELECTED')chosen<-c(chosen,id)
   orders<-c(orders,if(decision=='SELECTED')length(chosen) else NA_integer_);reasons<-c(reasons,decision);maxjac<-c(maxjac,jj)
  }
  d[,':='(direction=dr,display_order=orders,display_decision=reasons,max_Jaccard_previous=maxjac)]
  audit[[paste(ct,dr)]]<-d[,.(cell_type,direction,pathway,NES,q_global,display_order,display_decision,max_Jaccard_previous)]
  dd<-d[display_decision=='SELECTED'];dd[,measured_genes:=vapply(pathway,function(id)paste(sort(measured[[id]]),collapse=';'),character(1))]
  selected[[paste(ct,dr)]]<-dd
  counts[[paste(ct,dr)]]<-data.table(cell_type=ct,direction=dr,significant_GO=nrow(d),displayed_GO=nrow(dd))
  ck(paste(ct,dr,'display cap and q'),nrow(dd)<=8&&all(dd$q_global<.05))
  if(nrow(dd)>1){comb<-combn(dd$pathway,2);js<-apply(comb,2,function(x)length(intersect(measured[[x[1]]],measured[[x[2]]]))/length(union(measured[[x[1]]],measured[[x[2]]])));ck(paste(ct,dr,'pairwise display overlap'),all(js<.5))}
 }
 ids<-unique(c(unlist(lapply(selected,function(d)if(nrow(d)&&d$cell_type[1]==ct)d$pathway else character())),focus$pathway[1:4]))
 for(id in ids){
  hit<-names(rank)%in%sets[[id]];coef<-2^30/sum(abs(rank));if(coef>=1)coef<-floor(coef);w<-round(abs(rank)*coef)
  walk<-cumsum(ifelse(hit,w/sum(w[hit]),-1/sum(!hit)));ees<-if(abs(max(walk))>abs(min(walk)))max(walk) else min(walk);peak<-if(ees>0)which.max(walk) else which.min(walk)
  le<-if(ees>0)names(rank)[hit&seq_along(rank)<=peak]else names(rank)[hit&seq_along(rank)>peak]
  row<-g[cell_type==ct&scope=='post_group'&pathway==id]
  ck(paste(ct,id,'running ES and leading edge'),nrow(row)==1&&abs(ees-row$ES)<1e-9&&setequal(le,strsplit(row$leading_edge,';',fixed=TRUE)[[1]]),abs(ees-row$ES))
  es[[paste(ct,id)]]<-data.table(cell_type=ct,pathway=id,ES_source=row$ES,ES_verified=ees,error=abs(ees-row$ES),leading_edge_equal=TRUE)
 }
}
d<-rbindlist(selected);d<-d[,.(cell_type,direction,pathway,NES,ES,q_global,pval,size,leading_edge_n,leading_edge,measured_genes,display_order)]
wt(d,'plotdata_GO_overviews.csv');wt(rbindlist(audit),'GO_display_selection_audit.csv');wt(rbindlist(counts),'GO_display_counts.csv');wt(rbindlist(es),'selected_ES_checks.csv')
e<-merge(g[scope=='post_group'&pathway%in%focus$pathway[1:4]],ca[,.(cell_type,pathway,CAMERA_p=PValue,CAMERA_q=q_global)],by=c('cell_type','pathway'),all.x=TRUE)
e<-merge(e,ss[scope=='post_group',.(cell_type,pathway,score_beta=beta,score_ci_low=ci_low,score_ci_high=ci_high,score_p=p,score_q=q_eight_scores)],by=c('cell_type','pathway'),all.x=TRUE)
e[,display_order:=match(pathway,focus$pathway[1:4])];setorder(e,cell_type,display_order)
old<-fread(file.path(origin,'results/tables/07_primary_cross_method_summary.csv'));ix<-match(paste(e$cell_type,e$pathway),paste(old$cell_type,old$pathway))
for(k in c('NES','q_global','CAMERA_q','score_beta','score_q'))ck(paste('exact summary source',k),max(abs(e[[k]]-old[[k]][ix]))<1e-12)
wt(e[,.(cell_type,pathway,display_order,NES,q_global,pval,size,leading_edge_n,n,High,Low,CAMERA_p,CAMERA_q,score_beta,score_ci_low,score_ci_high,score_p,score_q)],'plotdata_Exact_Fig7_Pathways.csv')
wt(v[scope=='post_group'],'plotdata_Patient_Scores.csv')
# Compact sensitivity audit; cDC1 all-patient score NE is intentionally retained.
sens<-merge(g[pathway%in%focus$pathway[1:4],.(cell_type,scope,pathway,NES,q_global)],ss[,.(cell_type,scope,pathway,score_status=status,score_q=q_eight_scores)],by=c('cell_type','scope','pathway'),all.x=TRUE)
wt(sens,'exact_pathway_sensitivity.csv');wt(rbindlist(checks),'02_selection_checks.csv')
cat('SELECTION / ADDITIONAL VERIFICATION',length(checks),'PASS\n');print(rbindlist(counts));print(d[,.(cell_type,direction,pathway,NES,q_global)])
