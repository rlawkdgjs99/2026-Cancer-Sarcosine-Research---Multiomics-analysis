#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])),"common.R"))
logcon<-file(file.path(out,"logs/02_run.log"),"wt");sink(logcon,split=TRUE)
lock<-fread(file.path(out,"results/tables/01_exposure_lock.csv"));stopifnot(all(vapply(file.path(out,lock$file),sha,character(1))==lock$sha256))
stopifnot(all(fread(file.path(out,"results/tables/01_independent_exposure_checks.csv"))$pass))
stopifnot(!file.exists(file.path(out,"results/tables/02_GSEA_all_scopes.csv")))
a<-readRDS(file.path(out,"resources/allcell_prepared.rds"));md<-copy(a$metadata);md[,group:=factor(group,levels=c("Low","High"))]
m<-as.data.table(readRDS(member_path));stopifnot(all(m$db_version=="2026.1.Hs"))
sets<-lapply(split(m$gene_symbol,m$pathway),unique)
old<-fread(file.path(src,"04_Fig7_Immune_Programs_26.09.06/resources/Fig6_exact_MSigDB_2026.1.Hs_membership.csv"))
stopifnot(all(vapply(focus$pathway[1:4],function(id)setequal(sets[[id]],old[gs_name==id,gene_symbol]),logical(1))))
wt(focus,"00_fixed_focus.csv")
prep<-list()
for(ct in c("CD8","cDC1")){
 floor<-if(ct=="CD8")50L else 5L;dt<-copy(md[get(paste0(ct,"_cells"))>=floor]);raw<-a$counts[[ct]][,dt$Sample,drop=FALSE]
 yy<-DGEList(raw);keep<-filterByExpr(yy,design=model.matrix(reformulate(covars(dt)),dt),min.count=10,min.total.count=15,large.n=10,min.prop=.7)
 yy<-calcNormFactors(yy[keep,,keep.lib.sizes=FALSE],method="TMM");eff<-yy$samples$lib.size*yy$samples$norm.factors
 dt[,":="(outcome_effective_library=eff,outcome_TMM_factor=yy$samples$norm.factors,outcome_retained_library=yy$samples$lib.size)]
 em<-log2(1+sweep(as.matrix(yy$counts),2,eff,"/")*1e6)
 em<-em[!rownames(em)%in%exposure_genes & apply(em,1,sd)>0,,drop=FALSE];ez<-t(scale(t(em)))
 scores<-matrix(NA_real_,nrow(dt),4,dimnames=list(dt$Patient,focus$pathway[1:4]));coverage<-list()
 for(id in focus$pathway[1:4]){
  genes<-intersect(sets[[id]],rownames(ez))
  if(length(genes)>=15){val<-colMeans(ez[genes,,drop=FALSE]);if(sd(val)>0)scores[,id]<-z(val)}
  coverage[[id]]<-data.table(cell_type=ct,pathway=id,score_genes=length(genes),genes=paste(genes,collapse=";"),status=if(all(is.finite(scores[,id])))"ESTIMABLE" else "NE")
 }
 prep[[ct]]<-list(patients=dt,y=yy,raw=raw,keep=keep,scores=scores,expression=em,coverage=rbindlist(coverage))
 wt(dt,paste0("01_",ct,"_outcome_normalization.csv"));wt(data.table(gene=rownames(raw),retained=keep),paste0("01_",ct,"_outcome_gene_filter.csv"))
}
saveRDS(prep,file.path(out,"resources/outcome_prepared.rds"))
wt(rbindlist(lapply(prep,function(p)p$coverage)),"01_score_gene_coverage.csv")
ans<-list();cams<-list();scoreres<-list();scorevals<-list();statuses<-list();warnings<-list();cover<-list()
scopes<-c("post_group","post_continuous","all_group","all_continuous")
for(ti in seq_along(prep)){
 ct<-names(prep)[ti];pp<-prep[[ct]]
 for(si in seq_along(scopes)){
  sc<-scopes[si];dt<-copy(pp$patients);if(startsWith(sc,"post"))dt<-dt[treatment=="Post"]
  ix<-match(dt$Patient,pp$patients$Patient);yy<-pp$y[,ix,keep.lib.sizes=TRUE]
  pred<-if(endsWith(sc,"continuous"))"exposure_z" else "group"
  X<-model.matrix(reformulate(c(covars(dt),pred)),dt);rownames(X)<-dt$Patient;j<-match(if(pred=="group")"groupHigh" else pred,colnames(X))
  ok<-nrow(dt)>=8&&min(table(dt$group))>=3&&qr(X)$rank==ncol(X)&&nrow(X)-ncol(X)>=5&&(ct!="CD8"||nrow(pp$patients)>=10)
  statuses[[paste(ct,sc)]]<-data.table(cell_type=ct,scope=sc,n=nrow(dt),n_cells=sum(dt[[paste0(ct,"_cells")]]),High=sum(dt$group=="High"),Low=sum(dt$group=="Low"),rank=qr(X)$rank,df=nrow(X)-qr(X)$rank,status=if(ok)"TESTED" else "NE")
  if(!ok)next
  stopifnot(all(as.character(dt$group)==as.character(md$group[match(dt$Patient,md$Patient)])),all(dt$degradation_mean_z==md$degradation_mean_z[match(dt$Patient,md$Patient)]))
  cat("\nMODEL",ct,sc,"patients",nrow(dt),"cells",sum(dt[[paste0(ct,"_cells")]]),"genes",nrow(yy),"\n");flush.console()
  vo<-voom(yy,design=X,lib.size=dt$outcome_effective_library,normalize.method="none",span=.5,adaptive.span=TRUE,save.plot=TRUE,plot=FALSE)
  fit<-eBayes(lmFit(vo,X),trend=FALSE,robust=FALSE)
  rk<-data.table(gene=rownames(fit$coefficients),beta=fit$coefficients[,j],t=fit$t[,j],p=fit$p.value[,j])
  rk[,excluded:=gene%in%exposure_genes];setorderv(rk,c("t","gene"),c(-1L,1L))
  stopifnot(all(is.finite(rk$t)),!anyDuplicated(rk$gene))
  rank<-setNames(rk[excluded==FALSE,t],rk[excluded==FALSE,gene]);rk[,q_gene:=p.adjust(p,"BH")]
  fwrite(rk,file.path(out,"results/ranks",paste0(ct,"_",sc,".csv.gz")),compress="gzip")
  saveRDS(list(cell_type=ct,scope=sc,patients=dt,y=yy,design=X,coefficient=j,voom=vo,fit=fit,rank=rank),file.path(out,"resources/models",paste0(ct,"_",sc,".rds")))
  for(id in focus$pathway){
   ng<-sum(names(rank)%in%sets[[id]])
   cover[[paste(ct,sc,id)]]<-data.table(cell_type=ct,scope=sc,pathway=id,available_genes=ng,status=if(ng>=15&&ng<=500)"TESTED" else "NE_GENE_SET_SIZE")
  }
  for(ci in seq_along(collections)){
   co<-collections[ci];paths<-lapply(split(m[collection==co,gene_symbol],m[collection==co,pathway]),unique)
   sz<-vapply(paths,function(g)sum(names(rank)%in%g),integer(1));paths<-paths[sz>=15 & sz<=500]
   seed<-860908+ti*10000+si*100+ci;set.seed(seed);cat(" GSEA",co,length(paths),"seed",seed,"\n");flush.console()
   aa<-withCallingHandlers(fgseaMultilevel(paths,rank,minSize=15,maxSize=500,eps=0,scoreType="std",gseaParam=1,sampleSize=101,nPermSimple=10000,nproc=1,BPPARAM=BiocParallel::SerialParam(progressbar=FALSE)),
    warning=function(w){warnings[[length(warnings)+1L]]<<-data.table(cell_type=ct,scope=sc,collection=co,message=conditionMessage(w));invokeRestart("muffleWarning")})
   aa<-as.data.table(aa);setnames(aa,"padj","q_collection_cell")
   aa[,":="(cell_type=ct,scope=sc,collection=co,seed=seed,n=nrow(dt),High=sum(dt$group=="High"),Low=sum(dt$group=="Low"),leading_edge_n=lengths(leadingEdge),leading_edge=vapply(leadingEdge,paste,collapse=";",FUN.VALUE=character(1)))]
   aa[,leadingEdge:=NULL];stopifnot(setequal(aa$pathway,names(paths)));ans[[paste(ct,sc,co)]]<-aa
   # Per-collection checkpoint records the only executed seed, without rerunning completed estimates.
   wt(aa,paste0("checkpoint_",ct,"_",sc,"_",gsub(":","",co),".csv"))
   if(sc=="post_group"){
    vv<-vo[match(names(rank),rownames(vo$E)),];idx<-lapply(paths,function(g)which(rownames(vv$E)%in%g))
    ca<-as.data.table(camera(vv,index=idx,design=X,contrast=j,inter.gene.cor=NA,allow.neg.cor=FALSE,sort=FALSE),keep.rownames="pathway")
    ca[,":="(cell_type=ct,scope=sc,collection=co)];cams[[paste(ct,co)]]<-ca
   }
  }
  for(id in focus$pathway[1:4]){
   val<-pp$scores[ix,id];ss<-hc3(val,X,j);ss[,":="(cell_type=ct,scope=sc,pathway=id,n=nrow(dt),df=nrow(X)-ncol(X))]
   scoreres[[paste(ct,sc,id)]]<-ss
   scorevals[[paste(ct,sc,id)]]<-data.table(cell_type=ct,scope=sc,pathway=id,Patient=dt$Patient,Sample=dt$Sample,group=dt$group,treatment=dt$treatment,histology=dt$histology,n_cells=dt[[paste0(ct,"_cells")]],score=val,exposure=dt$degradation_mean_z)
  }
 }
}
aa<-rbindlist(ans);aa[,q_global:=p.adjust(pval,"BH",n=.N),by=scope];aa[,q_cell_all_collections:=p.adjust(pval,"BH",n=.N),by=.(scope,cell_type)]
ca<-rbindlist(cams);ca[,q_global:=p.adjust(PValue,"BH",n=.N),by=scope]
ss<-rbindlist(scoreres);ss[,q_eight_scores:=p.adjust(p,"BH",n=.N),by=scope]
cv<-rbindlist(cover);zz<-merge(cv,aa,by=c("cell_type","scope","pathway"),all.x=TRUE);zz<-merge(zz,focus,by="pathway",all.x=TRUE)
wt(aa,"02_GSEA_all_scopes.csv");wt(ca,"03_CAMERA_primary.csv");wt(ss,"04_exact_scores_HC3.csv")
wt(rbindlist(scorevals),"04_patient_score_values.csv");wt(rbindlist(statuses),"01_model_eligibility.csv");wt(cv,"05_pathway_coverage.csv");wt(zz,"05_focus_complete.csv")
wt(if(length(warnings))rbindlist(warnings)else data.table(cell_type=character(),scope=character(),collection=character(),message=character()),"02_warnings.csv")
wt(aa[,.(n_tests=.N,n_valid_p=sum(is.finite(pval)),positive_q05=sum(NES>0&q_global<.05,na.rm=TRUE),negative_q05=sum(NES<0&q_global<.05,na.rm=TRUE)),by=.(scope,cell_type,collection)],"02_collection_summary.csv")
stopifnot(all(vapply(file.path(out,lock$file),sha,character(1))==lock$sha256))
capture.output(sessionInfo(),file=file.path(out,"logs/sessionInfo_run.txt"))
print(zz[role=="Fig6 exact",.(cell_type,scope,label,available_genes,status,NES,pval,q_global)])
cat("\nOUTCOME ANALYSIS COMPLETE. Independent verification/visual QA pending.\n");sink();close(logcon)
