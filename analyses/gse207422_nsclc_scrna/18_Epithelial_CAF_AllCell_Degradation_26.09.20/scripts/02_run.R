#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])),'common.R'))
stopifnot(all(fread(file.path(out,'results/tables/01_preflight_checks.csv'))$pass))
pp<-readRDS(file.path(out,'resources/outcome_prepared.rds'));m<-pp$membership;sets<-pp$sets
results<-list();warns<-list();audits<-list()
scopes<-'post_group'
for(si in seq_along(scopes)){
 sc<-scopes[si];dt<-copy(pp$patients);if(startsWith(sc,'post'))dt<-dt[treatment=='Post'];ix<-match(dt$Patient,pp$patients$Patient);yy<-pp$y[,ix,keep.lib.sizes=TRUE]
 pred<-if(endsWith(sc,'continuous'))'exposure_z'else'group';X<-model.matrix(reformulate(c(covars(dt),pred)),dt);rownames(X)<-dt$Patient;j<-match(if(pred=='group')'groupHigh'else pred,colnames(X))
 stopifnot(nrow(dt)>=8,min(table(dt$group))>=3,qr(X)$rank==ncol(X),nrow(X)-ncol(X)>=5)
 vo<-voom(yy,design=X,lib.size=dt$outcome_effective_library,normalize.method='none',span=.5,adaptive.span=TRUE,save.plot=TRUE,plot=FALSE)
 fit<-eBayes(lmFit(vo,X),trend=FALSE,robust=FALSE)
 rk<-data.table(gene=rownames(fit$coefficients),beta=fit$coefficients[,j],t=fit$t[,j],p=fit$p.value[,j]);rk[,':='(excluded=gene%in%exposure_genes,q_gene=p.adjust(p,'BH'))];setorderv(rk,c('t','gene'),c(-1L,1L))
 stopifnot(all(is.finite(rk$t)),!anyDuplicated(rk$gene));rank<-setNames(rk[excluded==FALSE,t],rk[excluded==FALSE,gene]);audits[[sc]]<-data.table(scope=sc,ranked_genes=length(rank),ties=sum(duplicated(rank)),positive=sum(rank>0),negative=sum(rank<0),excluded_present=paste(intersect(exposure_genes,rk$gene),collapse=';'))
 fwrite(rk,file.path(out,'results/ranks',paste0('Epithelial_',sc,'.csv.gz')),compress='gzip')
 saveRDS(list(scope=sc,patients=dt,y=yy,design=X,coefficient=j,voom=vo,fit=fit,rank=rank),file.path(out,'resources/models',paste0('Epithelial_',sc,'.rds')))
 wt(data.table(Patient=dt$Patient,X),paste0('design_',sc,'.csv'))
 cat('\nMODEL',sc,'patients',nrow(dt),'genes',length(rank),'\n');flush.console()
 for(ci in seq_along(collections)){
  co<-collections[ci];paths<-lapply(split(m[collection==co,gene_symbol],m[collection==co,pathway]),unique);sz<-vapply(paths,function(g)sum(names(rank)%in%g),integer(1));paths<-paths[sz>=15&sz<=500]
  seed<-26092000+ci;set.seed(seed)
  chk<-file.path(out,'results/tables',paste0('checkpoint_',sc,'_',gsub(':','',co),'.csv'))
  cat(' GSEA',co,length(paths),'seed',seed,'\n');flush.console()
  if(file.exists(chk)){
   aa<-fread(chk);stopifnot(all(aa$seed==seed),setequal(aa$pathway,names(paths)))
  }else{
   aa<-withCallingHandlers(fgseaMultilevel(paths,rank,minSize=15,maxSize=500,eps=0,scoreType='std',gseaParam=1,sampleSize=101,nPermSimple=10000,nproc=1,BPPARAM=BiocParallel::SerialParam(progressbar=FALSE)),warning=function(w){warns[[length(warns)+1L]]<<-data.table(scope=sc,collection=co,message=conditionMessage(w));invokeRestart('muffleWarning')})
   aa<-as.data.table(aa);setnames(aa,'padj','q_collection');aa[,':='(scope=sc,collection=co,seed=seed,n=nrow(dt),High=sum(dt$group=='High'),Low=sum(dt$group=='Low'),n_cells=sum(dt$lineage_cells),leading_edge_n=lengths(leadingEdge),leading_edge=vapply(leadingEdge,paste,collapse=';',FUN.VALUE=character(1)))];aa[,leadingEdge:=NULL];wt(aa,basename(chk))
  }
  stopifnot(setequal(aa$pathway,names(paths)),all(is.finite(aa$NES)),all(is.finite(aa$pval)),all(aa$pval>0&aa$pval<=1))
  results[[paste(sc,co)]]<-aa
 }
}
ans<-rbindlist(results);ans[,q_full:=p.adjust(pval,'BH',n=.N),by=scope]
wt(ans,'02_GSEA_all_scopes.csv');wt(rbindlist(audits),'02_rank_audit.csv')
wt(if(length(warns))rbindlist(warns)else data.table(scope=character(),collection=character(),message=character()),'02_GSEA_warnings.csv')
exact<-ans[pathway%in%focus];exact[,display_order:=match(pathway,focus)];setorder(exact,scope,display_order);wt(exact,'02_exact_pathway_sensitivity.csv')
e<-exact[scope=='post_group'];e[,':='(significant=q_full<.05,q_label=ifelse(q_full<.001,formatC(q_full,format='e',digits=1),sprintf('%.3f',q_full)))];wt(e,'plotdata_Exact_Fig7_Pathways.csv')

wt(ans[,.(tested=.N,significant=sum(q_full<.05),positive_significant=sum(q_full<.05&NES>0),negative_significant=sum(q_full<.05&NES<0)),by=.(scope,collection)],'02_GSEA_summary.csv')
capture.output(sessionInfo(),file=file.path(out,'logs/sessionInfo_analysis.txt'))
cat('ANALYSIS COMPLETE\n');print(e[,.(pathway,NES,q_full,n,High,Low,n_cells)])
