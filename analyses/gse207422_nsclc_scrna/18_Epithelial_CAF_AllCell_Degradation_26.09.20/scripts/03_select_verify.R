#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])),'common.R'))
a<-fread(file.path(out,'results/tables/02_GSEA_all_scopes.csv'));pp<-readRDS(file.path(out,'resources/outcome_prepared.rds'));model<-readRDS(file.path(out,'resources/models/Epithelial_post_group.rds'));u<-names(model$rank)
a[,':='(absNES=abs(NES),bh_family_n=.N,status=ifelse(NES<=0,'LOW_DIRECTION',ifelse(q_full>=.05,'BH_Q_GE_0.05','CANDIDATE')),max_jaccard=as.numeric(NA))]
d<-a[collection%in%c('Hallmark','Reactome')];setorderv(d,c('q_full','absNES','pathway'),c(1L,-1L,1L));chosen<-character();selected<-list()
for(k in which(d$status=='CANDIDATE')){
 id<-d$pathway[k];genes<-intersect(pp$sets[[id]],u);ov<-if(length(chosen))vapply(chosen,function(x){prev<-intersect(pp$sets[[x]],u);length(intersect(genes,prev))/length(union(genes,prev))},numeric(1))else 0
 why<-if(length(chosen)>=8)'DISPLAY_LIMIT_8'else if(max(ov)>=.5)'GENE_OVERLAP_GE_0.50'else'SELECTED';d$status[k]<-why;d$max_jaccard[k]<-max(ov)
 if(why=='SELECTED'){chosen<-c(chosen,id);r<-copy(d[k]);r[,':='(display_order=length(chosen),measured_genes=paste(sort(genes),collapse=';'),neglog10_BH_q=-log10(q_full))];selected[[length(selected)+1L]]<-r}
}
sel<-if(length(selected))rbindlist(selected)else d[0];wt(sel,'03_HR_selected_plotdata.csv');wt(d,'03_HR_selection_audit.csv')
checks<-list();ck<-function(k,v,detail=''){checks[[length(checks)+1L]]<<-data.table(check=k,pass=isTRUE(v),detail=as.character(detail));if(!isTRUE(v))stop(k)}
# Direct weighted least squares, independently from lmFit.
X<-model$design;coefj<-model$coefficient
betas<-vapply(seq_len(nrow(model$voom$E)),function(i){w<-model$voom$weights[i,];solve(crossprod(X,w*X),crossprod(X,w*model$voom$E[i,]))[coefj]},numeric(1))
err<-max(abs(betas-model$fit$coefficients[,coefj]));ck('Independent weighted least squares coefficients',err<1e-8,err)
tref<-model$fit$coefficients[,coefj]/(model$fit$stdev.unscaled[,coefj]*sqrt(model$fit$s2.post));ck('Moderated t identity',max(abs(tref-model$fit$t[,coefj]))<1e-10)
# Independent BH implementation.
ord<-order(a$pval);qs<-pmin(1,rev(cummin(rev(a$pval[ord]*nrow(a)/seq_len(nrow(a))))));back<-numeric(nrow(a));back[ord]<-qs;ck('Independent full-family BH',max(abs(back-a$q_full))<1e-14,nrow(a))
# Reconstruct the full running enrichment statistic and leading edge for displayed positive sets.
for(id in chosen){
 hits<-names(model$rank)%in%pp$sets[[id]];steps<-ifelse(hits,abs(model$rank)/sum(abs(model$rank[hits])), -1/sum(!hits));walk<-cumsum(steps);es<-max(walk);at<-which.max(walk);le<-names(model$rank)[which(hits&seq_along(hits)<=at)]
 r<-a[pathway==id];ck(paste('Running ES',id),abs(es-r$ES)<1e-6,sprintf('Full double-precision walk versus fgsea batch ES: abs error %.3g; numerical tolerance 1e-6',abs(es-r$ES)));ck(paste('Leading edge',id),setequal(le,strsplit(r$leading_edge,';',fixed=TRUE)[[1]]))
}
ck('Selected sets positive and significant',all(sel$NES>0&sel$q_full<.05));ck('At most eight and no forced program priority',nrow(sel)<=8)
# Input hashes including the recorded plan are unchanged.
man<-fread(file.path(out,'results/tables/00_input_manifest.csv'));ck('Input hashes unchanged',all(vapply(man$path,sha,character(1))==man$sha256))
wt(rbindlist(checks),'04_independent_verification.csv');cat('VERIFIED\n');print(sel[,.(pathway,NES,q_full,leading_edge_n)]);print(a[,.(tested=.N,High_significant=sum(NES>0&q_full<.05),Low_significant=sum(NES<0&q_full<.05)),by=collection])
