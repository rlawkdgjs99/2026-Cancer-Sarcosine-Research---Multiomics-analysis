#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])),"common.R"))
cks<-list();ck<-function(name,pass,delta=NA_real_){cks[[name]]<<-data.table(check=name,pass=as.logical(pass),max_difference=delta);cat(name,pass,delta,'\n');flush.console();stopifnot(pass)}
bh<-function(p){ans<-rep(NA_real_,length(p));ii<-which(is.finite(p));oo<-ii[order(p[ii])];ans[oo]<-pmin(1,rev(cummin(rev(p[oo]*length(p)/seq_along(oo)))));ans}
g<-fread(file.path(out,'results/tables/02_GSEA_all_scopes.csv'));ca<-fread(file.path(out,'results/tables/03_CAMERA_primary.csv'));ss<-fread(file.path(out,'results/tables/04_exact_scores_HC3.csv'));sv<-fread(file.path(out,'results/tables/04_patient_score_values.csv'))
md<-fread(file.path(out,'results/tables/01_frozen_patient_groups.csv'));m<-as.data.table(readRDS(member_path));sets<-lapply(split(m$gene_symbol,m$pathway),unique)
ck('GSEA_four_scopes',setequal(unique(g$scope),c('post_group','post_continuous','all_group','all_continuous')));ck('no_duplicate_hypotheses',!anyDuplicated(g[,.(scope,pathway)]));ck('finite_p_NES',all(is.finite(g$pval))&&all(g$pval>0&g$pval<=1)&&all(is.finite(g$NES)))
for(sc in unique(g$scope)){dt<-g[scope==sc];de<-max(abs(bh(dt$pval)-dt$q_global));ck(paste0(sc,'_BH_global'),de<1e-12,de);for(co in unique(dt$collection)){v<-dt[collection==co];de<-max(abs(bh(v$pval)-v$q_collection_cell));ck(paste(sc,co,'BH_collection',sep='_'),de<1e-12,de)}}
ck('all_scope_BH',max(abs(bh(g$pval)-g$q_all_scopes))<1e-12,max(abs(bh(g$pval)-g$q_all_scopes)));ck('CAMERA_BH',max(abs(bh(ca$PValue)-ca$q_global))<1e-12,max(abs(bh(ca$PValue)-ca$q_global)))
for(sc in unique(ss$scope)){v<-ss[scope==sc];ck(paste(sc,'score_BH'),max(abs(bh(v$p)-v$q_four_scores))<1e-12,max(abs(bh(v$p)-v$q_four_scores)))}
# Freeze display selection using the prespecified rule. No biological-theme preferences.
mod<-readRDS(file.path(out,'resources/models/All_cells_post_group.rds'));universe<-names(mod$rank);measured<-lapply(sets,intersect,y=universe)
primary<-g[scope=='post_group'];candidates<-primary[collection=='GO:BP'&q_global<.05];candidates[,absNES:=abs(NES)];setorderv(candidates,c('q_global','absNES','pathway'),c(1L,-1L,1L))
sel<-list();audit<-list()
for(direction in c('High','Low')){v<-candidates[if(direction=='High') NES>0 else NES<0];chosen<-character();for(i in seq_len(nrow(v))){id<-v$pathway[i];ov<-if(length(chosen))vapply(chosen,function(j)length(intersect(measured[[id]],measured[[j]]))/length(union(measured[[id]],measured[[j]])),numeric(1)) else numeric();mx<-if(length(ov))max(ov) else 0;reason<-if(length(chosen)>=8)'DISPLAY_LIMIT' else if(mx>=.5)'GENE_OVERLAP_GE_0.50' else 'SELECTED';audit[[paste(direction,id)]]<-data.table(direction=direction,pathway=id,max_Jaccard_to_selected=mx,overlap_with=if(length(ov))names(which.max(ov)) else '',status=reason);if(reason=='SELECTED')chosen<-c(chosen,id)}
 zt<-v[match(chosen,pathway)];zt[,':='(direction=direction,display_order=seq_len(.N))];sel[[direction]]<-zt}
plotdt<-rbindlist(sel);wt(rbindlist(audit),'03_GO_display_selection_audit.csv');wt(plotdt,'plotdata_Fig8e_GO_overview.csv');wt(data.table(direction=c('High','Low'),significant_GOBP=vapply(c('High','Low'),function(d)nrow(candidates[if(d=='High') NES>0 else NES<0]),integer(1)),displayed=vapply(sel,nrow,integer(1))),'03_GO_display_counts.csv')
exact<-merge(primary[pathway%in%focus$pathway[1:4]],focus[,.(pathway,label,role)],by='pathway');exact<-merge(exact,ca[,.(pathway,CAMERA_p=PValue,CAMERA_q=q_global)],by='pathway');exact<-merge(exact,ss[scope=='post_group',.(pathway,score_beta=beta,score_ci_low=ci_low,score_ci_high=ci_high,score_p=p,score_q=q_four_scores)],by='pathway');exact<-exact[match(focus$pathway[1:4],pathway)];wt(exact,'plotdata_Fig8e_Exact_Fig7_Pathways.csv')
walks<-list()
for(sc in unique(g$scope)){
 mo<-readRDS(file.path(out,'resources/models',paste0('All_cells_',sc,'.rds')));X<-mo$design;E<-mo$voom$E;W<-mo$voom$weights;j<-mo$coefficient;df<-nrow(X)-ncol(X)
 ck(paste(sc,'fixed_labels'),all(as.character(mo$patients$group)==md$group[match(mo$patients$Patient,md$Patient)]));ck(paste(sc,'no_exposure_rank_genes'),!any(exposure_genes%in%names(mo$rank)))
 ind<-t(vapply(seq_len(nrow(E)),function(i){inv<-solve(crossprod(X,X*W[i,]));bb<-inv%*%crossprod(X,W[i,]*E[i,]);er<-E[i,]-as.numeric(X%*%bb);c(beta=bb[j],unscaled=sqrt(inv[j,j]),s2=sum(W[i,]*er^2)/df)},numeric(3)))
 d0<-mo$fit$df.prior;s0<-mo$fit$s2.prior;post<-(df*ind[,'s2']+d0*s0)/(df+d0);if(length(d0)==1&&is.infinite(d0))post<-rep_len(s0,nrow(E));tt<-ind[,'beta']/sqrt(post)/ind[,'unscaled']
 for(nm in c('beta','unscaled','s2','posterior','t')){ref<-switch(nm,beta=mo$fit$coefficients[,j],unscaled=mo$fit$stdev.unscaled[,j],s2=mo$fit$sigma^2,posterior=mo$fit$s2.post,t=mo$fit$t[,j]);val<-switch(nm,beta=ind[,'beta'],unscaled=ind[,'unscaled'],s2=ind[,'s2'],posterior=post,t=tt);delta<-max(abs(val-ref));ck(paste(sc,'independent_WLS',nm),delta<1e-8,delta)}
 wt(data.table(Patient=rownames(X),X),paste0('verification_design_',sc,'.csv'))
 paths<-if(sc=='post_group')unique(c(plotdt$pathway,focus$pathway[1:4])) else focus$pathway[1:4]
 rk<-mo$rank
 # fgsea1.38 prepareStats integerizes absolute rank weights internally; reproduce that kernel independently.
 scaleCoeff<-2^30/sum(abs(rk));if(scaleCoeff>=1)scaleCoeff<-floor(scaleCoeff);aw<-round(abs(rk)*scaleCoeff)
 for(id in paths){hits<-names(rk)%in%sets[[id]];nh<-sum(hits);step<-ifelse(hits,aw/sum(aw[hits]),-1/(length(rk)-nh));walk<-cumsum(step);es<-if(max(walk)> -min(walk))max(walk) else min(walk);ref<-g[scope==sc&pathway==id,ES];de<-abs(es-ref);ck(paste(sc,'ES',id),length(de)==1&&de<1e-10,de);k<-if(es>0)which.max(walk) else which.min(walk);le<-if(es>0)names(rk)[hits&seq_along(rk)<=k] else names(rk)[hits&seq_along(rk)>k];refle<-strsplit(g[scope==sc&pathway==id,leading_edge],';',fixed=TRUE)[[1]];ck(paste(sc,'leading_edge',id),setequal(le,refle));if(sc=='post_group'&&id%in%focus$pathway[1:4])walks[[id]]<-data.table(pathway=id,rank=seq_along(rk),running_ES=walk,hit=hits)}
}
wt(rbindlist(walks),'03_exact_running_ES.csv')
locks<-fread(file.path(out,'results/tables/00_input_manifest.csv'));ck('input_hashes_preserved',all(vapply(file.path(src,locks$path),sha,character(1))==locks$sha256));locks<-fread(file.path(out,'results/tables/00_protected14_before.csv'));ck('all_existing14_files_preserved',all(vapply(file.path(src,locks$path),sha,character(1))==locks$sha256))
wt(rbindlist(cks),'03_numerical_verification.csv');cat('VERIFICATION PASS',length(cks),'checks\n');print(exact[,.(label,NES,q_global,CAMERA_q,score_q)]);print(plotdt[,.(direction,pathway,NES,q_global,leading_edge_n)])
