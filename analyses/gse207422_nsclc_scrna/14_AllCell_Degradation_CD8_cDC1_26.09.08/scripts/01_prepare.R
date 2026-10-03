#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])),"common.R"))
logcon<-file(file.path(out,"logs/01_prepare.log"),"wt");sink(logcon,split=TRUE)
plan<-file.path(out,"ANALYSIS_PLAN_FROZEN.md");lock<-file.path(out,"logs/frozen_plan_sha256.txt")
stopifnot(!file.exists(file.path(out,"resources/allcell_prepared.rds")))
if(file.exists(lock))stopifnot(identical(sha(plan),readLines(lock)))else writeLines(sha(plan),lock)
p02<-file.path(src,"02_lineage_reannotation");p11<-file.path(src,"11_cDC1_Identity_Audit_26.09.08")
paths<-c(raw=file.path(p02,"intermediate/01_full_counts_mt20_qc.rds"),
 cells=file.path(p02,"results/tables/09_final_cell_lineages_FROZEN.csv"),
 metadata=file.path(prior,"resources/verified_patient_pseudobulks.rds"),
 cDC1roster=file.path(p11,"results/tables/cDC1_candidate_cell_roster.csv"),
 oldCD8=file.path(p02,"intermediate/10_lineage_pseudobulk_counts.rds"),
 oldcDC1=file.path(src,"12_cDC1_CD8_Defined_Groups_Exploratory_26.09.08/resources/patient_pseudobulk.rds"),
 membership=member_path,focus=file.path(prior,"results/tables/00_prespecified_focus.csv"),
 exact=file.path(src,"04_Fig7_Immune_Programs_26.09.06/resources/Fig6_exact_MSigDB_2026.1.Hs_membership.csv"))
stopifnot(all(file.exists(paths)));cat("Hashing protected inputs\n")
wt(data.table(role=names(paths),path=substring(paths,nchar(src)+2),sha256=vapply(paths,sha,character(1))),"00_inputs.csv")
protected<-unlist(lapply(c("09_Fig7_Bulk_Aligned_MeanZ_26.09.08","10_cDC_Bulk_Aligned_MeanZ_Exploratory_26.09.08","11_cDC1_Identity_Audit_26.09.08","12_cDC1_CD8_Defined_Groups_Exploratory_26.09.08","13_cDC1_Abundance_CD8_Defined_Groups_26.09.08"),function(p)list.files(file.path(src,p),recursive=TRUE,full.names=TRUE)))
wt(data.table(path=substring(protected,nchar(src)+2),sha256=vapply(protected,sha,character(1))),"00_protected_before.csv")
f<-fread(paths["cells"]);r<-fread(paths["cDC1roster"]);a<-readRDS(paths["metadata"])
md<-as.data.table(a$meta);setorder(md,Sample);md[,group:=NULL]
stopifnot(nrow(md)==15,!anyDuplicated(md$Patient),!anyDuplicated(md$Sample),nrow(f)==92053,!anyDuplicated(f$cell_id),nrow(r)==129,!anyDuplicated(r$cell_id))
stopifnot(all(r$cell_id%in%f[final_lineage=="Conventional DC",cell_id]))
cat("Reading frozen full raw sparse matrix; no dense cell matrix is created\n");flush.console()
x<-readRDS(paths["raw"]);stopifnot(inherits(x,"dgCMatrix"),ncol(x)==92053,!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)))
f<-f[match(colnames(x),cell_id)];stopifnot(identical(colnames(x),f$cell_id),!anyNA(f$Patient),all(Matrix::colSums(x)==f$library_size))
for(st in seq.int(1L,length(x@x),by=5000000L)){ix<-st:min(length(x@x),st+4999999L);v<-x@x[ix];stopifnot(all(is.finite(v)),all(v>=0),all(v==floor(v)))}
keep<-f$final_lineage!="Excluded residual doublet";stopifnot(sum(keep)==91844)
pb<-list();celln<-list()
for(nm in c("All_cells","CD8","cDC1")){
 use<-switch(nm,All_cells=which(keep),CD8=which(f$final_lineage=="CD8 T cell"),cDC1=which(f$cell_id%in%r$cell_id))
 mem<-sparseMatrix(i=use,j=match(f$Sample[use],md$Sample),x=1,dims=c(ncol(x),nrow(md)),dimnames=list(colnames(x),md$Sample))
 pb[[nm]]<-as.matrix(x%*%mem)
 stopifnot(all(colSums(pb[[nm]])==as.numeric(tapply(f$library_size[use],factor(f$Sample[use],levels=md$Sample),sum,default=0))))
 celln[[nm]]<-as.integer(tabulate(match(f$Sample[use],md$Sample),nbins=15))
 cat(nm,"cells",length(use),"UMI",sum(pb[[nm]]),"\n")
}
old<-readRDS(paths["oldCD8"])[["CD8 T cell"]];stopifnot(all(pb$CD8==as.matrix(old[rownames(pb$CD8),md$Sample])))
oldc<-readRDS(paths["oldcDC1"]);oldcm<-oldc$counts
stopifnot(all(pb$cDC1[,match(colnames(oldcm),md$Patient),drop=FALSE]==oldcm[rownames(pb$cDC1),,drop=FALSE]))
unres<-which(f$final_lineage=="NK/gamma-delta T unresolved")
umat<-sparseMatrix(i=unres,j=match(f$Sample[unres],md$Sample),x=1,dims=c(ncol(x),15))
upb<-as.matrix(x%*%umat);oldall<-Reduce("+",readRDS(paths["oldCD8"]))
stopifnot(all(pb$All_cells==as.matrix(oldall[rownames(pb$All_cells),md$Sample])+upb))
target<-as.matrix(x[c("SARDH","PIPOX"),,drop=FALSE])
contrib<-rbindlist(lapply(1:2,function(i)data.table(Sample=f$Sample,Patient=f$Patient,lineage=f$final_lineage,gene=rownames(target)[i],raw=as.numeric(target[i,]))[keep,.(raw=sum(raw)),by=.(Sample,Patient,lineage,gene)]))
wt(contrib,"00_gene_lineage_contributions.csv")
composition<-f[keep,.(cells=.N),by=.(Patient,Sample,final_lineage)]
composition[,fraction:=cells/sum(cells),by=Patient];wt(composition,"00_patient_cell_composition.csv")
md[,":="(all_cells=celln$All_cells,CD8_cells=celln$CD8,cDC1_cells=celln$cDC1)]
wt(data.table(gene=rownames(upb),upb),"00_unresolved_immune_pseudobulk.csv.gz")
saveRDS(list(counts=pb,metadata=md,celln=celln),file.path(out,"resources/raw_patient_pseudobulks.rds"))
rm(x,target,mem,umat,a,old,oldall,upb);gc()
make_score<-function(raw,dt){
 yy<-DGEList(raw);X0<-model.matrix(reformulate(covars(dt)),dt)
 keepgene<-filterByExpr(yy,design=X0,min.count=10,min.total.count=15,large.n=10,min.prop=.7)
 yy<-calcNormFactors(yy[keepgene,,keep.lib.sizes=FALSE],method="TMM");eff<-yy$samples$lib.size*yy$samples$norm.factors
 xx<-log2(1+sweep(raw[c("SARDH","PIPOX"),,drop=FALSE],2,eff,"/")*1e6);zz<-t(apply(xx,1,z))
 D<-colMeans(zz);cut<-median(D);stopifnot(all(is.finite(D)),sd(D)>0)
 list(y=yy,keep=keepgene,x=xx,z=zz,D=D,D_z=z(D),cutoff=cut,group=ifelse(D>cut,"High","Low"),effective=eff)
}
e<-make_score(pb$All_cells,md)
md[,":="(SARDH_raw=as.numeric(pb$All_cells["SARDH",]),PIPOX_raw=as.numeric(pb$All_cells["PIPOX",]),
 SARDH_log2CPM=as.numeric(e$x["SARDH",]),PIPOX_log2CPM=as.numeric(e$x["PIPOX",]),SARDH_z=as.numeric(e$z["SARDH",]),PIPOX_z=as.numeric(e$z["PIPOX",]),
 degradation_mean_z=as.numeric(e$D),exposure_z=e$D_z,cutoff=e$cutoff,group=e$group,effective_library=e$effective,TMM_factor=e$y$samples$norm.factors)]
wt(md,"01_FROZEN_allcell_patient_scores_groups.csv")
wt(data.table(gene=rownames(e$x),mean_log=rowMeans(e$x),sd_log=apply(e$x,1,sd),nonzero=rowSums(pb$All_cells[rownames(e$x),]>0),ge5=rowSums(pb$All_cells[rownames(e$x),]>=5)),"01_exposure_component_reference.csv")
wt(data.table(gene=rownames(pb$All_cells),retained=e$keep),"01_exposure_normalization_gene_filter.csv")
comp<-merge(composition,md[,.(Patient,degradation_mean_z)],by="Patient")
wt(comp[,.(rho_spearman=cor(fraction,degradation_mean_z,method="spearman"),n=.N),by=final_lineage],"01_composition_descriptive_correlations.csv")
diags<-list()
for(nm in c("CD8","cDC1")){
 dd<-make_score(pb$All_cells-pb[[nm]],md);diags[[nm]]<-data.table(excluded_compartment=nm,Patient=md$Patient,diagnostic_D=dd$D,diagnostic_group=dd$group,primary_D=md$degradation_mean_z,primary_group=md$group)
}
wt(rbindlist(diags),"01_leave_compartment_out_DIAGNOSTIC_ONLY.csv")
saveRDS(list(counts=pb,metadata=md,exposure=e),file.path(out,"resources/allcell_prepared.rds"))
wt(data.table(scope=c("All15","Post12"),n=c(15,sum(md$treatment=="Post")),High=c(sum(md$group=="High"),sum(md$group=="High"&md$treatment=="Post")),Low=c(sum(md$group=="Low"),sum(md$group=="Low"&md$treatment=="Post"))),"01_group_counts.csv")
lockpaths<-file.path(out,c("ANALYSIS_PLAN_FROZEN.md","results/tables/01_FROZEN_allcell_patient_scores_groups.csv","resources/allcell_prepared.rds"))
wt(data.table(file=substring(lockpaths,nchar(out)+2),sha256=vapply(lockpaths,sha,character(1))),"01_exposure_lock.csv")
stopifnot(identical(sha(plan),readLines(lock)))
capture.output(sessionInfo(),file=file.path(out,"logs/sessionInfo_prepare.txt"))
print(md[,.(Patient,treatment,histology,all_cells,CD8_cells,cDC1_cells,SARDH_raw,PIPOX_raw,degradation_mean_z,group)])
cat("PREPARATION COMPLETE. No outcome fit performed.\n");sink();close(logcon)
