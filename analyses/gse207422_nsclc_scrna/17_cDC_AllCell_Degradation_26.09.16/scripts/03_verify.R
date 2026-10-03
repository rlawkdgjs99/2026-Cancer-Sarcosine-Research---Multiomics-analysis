#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])),'common.R'))
checks<-list();ck<-function(k,v,d=''){checks[[length(checks)+1L]]<<-data.table(check=k,pass=isTRUE(v),detail=as.character(d));if(!isTRUE(v))stop(k,': ',d)}
a<-fread(file.path(out,'results/tables/02_GSEA_all_scopes.csv'));pp<-readRDS(file.path(out,'resources/outcome_prepared.rds'));md<-fread(file.path(legacy,'results/tables/01_FROZEN_allcell_patient_scores_groups.csv'));sets<-pp$sets
ck('unique complete hypotheses',!anyDuplicated(a[,.(scope,collection,pathway)])&&nrow(a)==4*3797)
bh<-function(p){o<-order(p);q<-numeric(length(p));q[o]<-pmin(1,rev(cummin(rev(p[o]*length(p)/seq_along(p)))));q}
for(sc in scopes){d<-a[scope==sc];ck(paste(sc,'independent global BH'),max(abs(bh(d$pval)-d$q_global_cDC))<1e-12);for(co in collections){r<-d[collection==co];ck(paste(sc,co,'collection BH'),max(abs(bh(r$pval)-r$q_collection))<1e-12)}}
X0<-pp$nuisance_design;h<-rowSums((X0%*%solve(crossprod(X0)))*X0);mn<-1/max(h);if(mn>10)mn<-10+(mn-10)*.7
keep<-rowSums(sweep(pp$raw,2,colSums(pp$raw),'/')*1e6>=10/median(colSums(pp$raw))*1e6)>=mn-1e-14&rowSums(pp$raw)>=15-1e-14
ck('independent nuisance gene filtering',identical(unname(keep),unname(pp$keep)))
xx<-pp$raw[keep,,drop=FALSE];libs<-colSums(xx);q75<-apply(xx,2,quantile,.75)/libs;ref<-if(median(q75)<1e-20)which.max(colSums(sqrt(xx)))else which.min(abs(q75-mean(q75)))
factors<-vapply(seq_len(ncol(xx)),function(i){o<-xx[,i];r<-xx[,ref];L<-log2((o/libs[i])/(r/libs[ref]));A<-(log2(o/libs[i])+log2(r/libs[ref]))/2;V<-(libs[i]-o)/libs[i]/o+(libs[ref]-r)/libs[ref]/r;ok<-is.finite(L)&is.finite(A)&A> -1e10;L<-L[ok];A<-A[ok];V<-V[ok];if(max(abs(L))<1e-6)return(1);n<-length(L);loL<-floor(n*.3)+1;loA<-floor(n*.05)+1;k<-rank(L)>=loL&rank(L)<=n+1-loL&rank(A)>=loA&rank(A)<=n+1-loA;2^(sum(L[k]/V[k])/sum(1/V[k]))},numeric(1));factors<-factors/exp(mean(log(factors)))
ck('independent TMM factors',max(abs(factors-pp$patients$outcome_TMM_factor))<1e-12);ck('outcome libraries unchanged',all(libs==pp$patients$outcome_retained_library)&&max(abs(libs*factors-pp$patients$outcome_effective_library))<1e-7)
plots<-fread(file.path(out,'results/tables/plotdata_GO_overview.csv'));exact<-fread(file.path(out,'results/tables/plotdata_Exact_Fig7_Pathways.csv'));ck('all four exact sets retained',nrow(exact)==4&&identical(exact$pathway,focus))
allchecks<-list();maxdel<-0
for(sc in scopes){b<-readRDS(file.path(out,'resources/models',paste0('cDC_',sc,'.rds')));X<-b$design;j<-b$coefficient;dt<-b$patients;E<-b$voom$E;W<-b$voom$weights
 ck(paste(sc,'patient design alignment'),identical(dt$Patient,rownames(X))&&identical(dt$Sample,colnames(E)))
 ck(paste(sc,'frozen exposure labels'),all(as.character(dt$group)==md$group[match(dt$Patient,md$Patient)])&&all(dt$degradation_mean_z==md$degradation_mean_z[match(dt$Patient,md$Patient)]))
 ck(paste(sc,'voom logCPM'),max(abs(E-log2(sweep(as.matrix(b$y$counts)+.5,2,dt$outcome_effective_library+1,'/')*1e6)))<1e-10)
 del<-0
 for(i in seq_len(nrow(E))){inv<-solve(crossprod(X,X*W[i,]));bt<-as.vector(inv%*%crossprod(X,E[i,]*W[i,]));s2<-sum(W[i,]*(E[i,]-as.vector(X%*%bt))^2)/(nrow(X)-ncol(X));prior<-b$fit$s2.prior;dfprior<-b$fit$df.prior;post<-if(is.infinite(dfprior))prior else (dfprior*prior+(nrow(X)-ncol(X))*s2)/(dfprior+nrow(X)-ncol(X));t<-bt[j]/sqrt(inv[j,j]*post);del<-max(del,max(abs(bt-b$fit$coefficients[i,])),abs(t-b$fit$t[i,j]),abs(s2-b$fit$sigma[i]^2),abs(post-b$fit$s2.post[i]))}
 ck(paste(sc,'all-gene WLS and moderated t'),del<1e-8,del);maxdel<-max(maxdel,del)
 rank<-b$rank;ck(paste(sc,'rank contract'),!any(exposure_genes%in%names(rank))&&!anyDuplicated(names(rank))&&all(diff(rank)<=0))
 sz<-vapply(sets,function(g)sum(names(rank)%in%g),integer(1));eligible<-names(sz)[sz>=15&sz<=500];r<-a[scope==sc];ck(paste(sc,'complete planned universe'),setequal(r$pathway,eligible)&&all(r$size==sz[r$pathway]))
 ids<-if(sc=='post_group')unique(c(focus,plots$pathway))else focus
 for(id in ids){hit<-names(rank)%in%sets[[id]];scalev<-2^30/sum(abs(rank));if(scalev>=1)scalev<-floor(scalev);w<-round(abs(rank)*scalev);walk<-cumsum(ifelse(hit,w/sum(w[hit]),-1/sum(!hit)));es<-if(abs(max(walk))>abs(min(walk)))max(walk)else min(walk);peak<-if(es>0)which.max(walk)else which.min(walk);le<-if(es>0)names(rank)[hit&seq_along(rank)<=peak]else names(rank)[hit&seq_along(rank)>peak];row<-r[pathway==id];allchecks[[paste(sc,id)]]<-data.table(scope=sc,pathway=id,error=abs(es-row$ES),leading_edge_match=setequal(le,strsplit(row$leading_edge,';',fixed=TRUE)[[1]]))}
 cat('Verified models and ES',sc,'\n');flush.console()
}
ees<-rbindlist(allchecks);ck('independent selected ES leading edges',max(ees$error)<1e-9&&all(ees$leading_edge_match),max(ees$error));wt(ees,'03_ES_checks.csv')
# Recreate display selection independently from complete primary results.
u<-names(readRDS(file.path(out,'resources/models/cDC_post_group.rds'))$rank)
for(dr in c('High','Low')){d<-a[scope=='post_group'&collection=='GO:BP'&q_global_cDC<.05&if(dr=='High')NES>0 else NES<0];d<-d[order(q_global_cDC,-abs(NES),pathway)];chosen<-character();for(id in d$pathway){genes<-intersect(sets[[id]],u);ov<-vapply(chosen,function(j){g<-intersect(sets[[j]],u);length(intersect(genes,g))/length(union(genes,g))},numeric(1));if(!length(ov)||all(ov<.5))chosen<-c(chosen,id);if(length(chosen)==8)break};ck(paste(dr,'independent GO selection'),identical(chosen,plots[direction==dr][order(display_order),pathway]))}
for(tab in list(plots,exact)){for(i in seq_len(nrow(tab))){row<-a[scope=='post_group'&pathway==tab$pathway[i]];ck(paste('plot matches source',tab$pathway[i]),all(unlist(tab[i,.(NES,q_global_cDC,leading_edge_n)])==unlist(row[,.(NES,q_global_cDC,leading_edge_n)])))}}
manifest<-fread(file.path(out,'results/tables/00_input_manifest.csv'));manifest[,sha_after:=vapply(file.path(src,path),sha,character(1))];ck('input and frozen plan unchanged',all(manifest$sha256==manifest$sha_after));wt(manifest,'03_input_preservation.csv')
prot<-fread(file.path(out,'results/tables/00_protected_before.csv'));prot[,sha_after:=vapply(file.path(src,path),sha,character(1))];ck('previous analyses unchanged',all(prot$sha256==prot$sha_after),nrow(prot));wt(prot,'03_previous_preservation.csv')
wt(rbindlist(checks),'03_verification_checks.csv');writeLines(c('PASS',paste('checks',length(checks)),paste('previous files preserved',nrow(prot)),paste('max WLS/t error',maxdel),paste('ES maximum error',max(ees$error))),file.path(out,'results/qa/verification_summary.txt'));cat('VERIFIED',length(checks),'checks\n')
