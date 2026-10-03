#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])),'common.R'))
checks<-list();ck<-function(k,v,detail=''){checks[[length(checks)+1L]]<<-data.table(check=k,pass=isTRUE(v),detail=as.character(detail));cat(k,isTRUE(v),detail,'\n');flush.console();if(!isTRUE(v))stop(k)}
paths<-c(raw=file.path(src,'02_lineage_reannotation/intermediate/01_full_counts_mt20_qc.rds'),cells=file.path(src,'02_lineage_reannotation/results/tables/09_final_cell_lineages_FROZEN.csv'),stored=file.path(src,'02_lineage_reannotation/intermediate/10_lineage_pseudobulk_counts.rds'),groups=file.path(legacy,'results/tables/01_FROZEN_allcell_patient_scores_groups.csv'),exposure=file.path(legacy,'resources/allcell_prepared.rds'),membership=member_path,plan=file.path(out,'ANALYSIS_PLAN.md'))
wt(data.table(role=names(paths),path=paths,sha256=vapply(paths,sha,character(1))),'00_input_manifest.csv')
f<-fread(paths['cells']);md<-fread(paths['groups']);a<-readRDS(paths['exposure']);md[,group:=factor(group,levels=c('Low','High'))]
ck('unique patient sample cell IDs',nrow(md)==15&&!anyDuplicated(md$Patient)&&!anyDuplicated(md$Sample)&&nrow(f)==92053&&!anyDuplicated(f$cell_id))
ck('frozen all-cell labels',identical(as.character(md$group),as.character(a$metadata$group))&&max(abs(md$degradation_mean_z-a$metadata$degradation_mean_z))<1e-14)
xlog<-log2(1+sweep(a$counts$All_cells[c('SARDH','PIPOX'),],2,a$exposure$effective,'/')*1e6);D<-colMeans(t(scale(t(xlog))))
ck('exposure reconstructed without regrouping',max(abs(D-md$degradation_mean_z))<1e-12&&all(ifelse(D>median(D),'High','Low')==md$group))
x<-readRDS(paths['raw']);f<-f[match(colnames(x),cell_id)];ck('raw cell alignment',identical(colnames(x),f$cell_id)&&nrow(x)==24292&&ncol(x)==92053&&all(Matrix::colSums(x)==f$library_size))
stored<-readRDS(paths['stored']);elig<-list();pbs<-list()
for(lineage in c('Epithelial','CAF')){
 use<-which(f$final_lineage==lineage);xc<-x[,use,drop=FALSE];ck(paste(lineage,'raw integer nonnegative'),all(is.finite(xc@x))&&all(xc@x>=0)&&all(xc@x==floor(xc@x)))
 mem<-sparseMatrix(i=seq_along(use),j=match(f$Sample[use],md$Sample),x=1,dims=c(length(use),15),dimnames=list(f$cell_id[use],md$Sample));pb<-as.matrix(xc%*%mem)
 ck(paste(lineage,'all genes/patients equal stored pseudobulk'),all(pb==as.matrix(stored[[lineage]][rownames(pb),md$Sample])))
 d<-copy(md);d[,lineage_cells:=tabulate(match(f$Sample[use],Sample),nbins=15)];d[,lineage:=lineage]
 ck(paste(lineage,'library sums'),all(colSums(pb)==tapply(f$library_size[use],factor(f$Sample[use],levels=md$Sample),sum)))
 d[,eligible_50:=lineage_cells>=50];elig[[lineage]]<-d;pbs[[lineage]]<-pb
 wt(f[use,.(cell_id,Patient,Sample,final_lineage,library_size)],paste0('00_',lineage,'_cell_roster.csv'))
}
eligibility<-rbindlist(elig);wt(eligibility,'01_patient_eligibility.csv');saveRDS(list(counts=pbs,metadata=eligibility),file.path(out,'resources/raw_verified.rds'))
rm(x,xc,mem,stored,a);gc()
dt<-elig[['Epithelial']][eligible_50==TRUE];ck('Epithelial all14/post11',nrow(dt)==14&&nrow(dt[treatment=='Post'])==11&&sum(dt[treatment=='Post']$group=='High')==6)
raw<-pbs[['Epithelial']][,dt$Sample,drop=FALSE];X0<-model.matrix(reformulate(covars(dt)),dt);yy<-DGEList(raw);keep<-filterByExpr(yy,design=X0,min.count=10,min.total.count=15,large.n=10,min.prop=.7);yy<-calcNormFactors(yy[keep,,keep.lib.sizes=FALSE],method='TMM');eff<-yy$samples$lib.size*yy$samples$norm.factors
ck('filtered normalized genes',sum(keep)>1000&&all(is.finite(eff))&&all(eff>0));dt[,':='(outcome_effective_library=eff,outcome_TMM_factor=yy$samples$norm.factors,outcome_retained_library=yy$samples$lib.size)]
wt(dt,'01_Epithelial_normalization.csv');wt(data.table(gene=rownames(raw),retained=keep),'01_gene_filter.csv')
d<-dt[treatment=='Post'];X<-model.matrix(reformulate(c(covars(d),'group')),d);ck('Epithelial primary model eligible',nrow(d)>=8&&min(table(d$group))>=3&&qr(X)$rank==ncol(X)&&nrow(X)-ncol(X)>=5)
wt(eligibility[treatment=='Post',.(patients=.N,cells=sum(lineage_cells),eligible_50=sum(eligible_50)),by=.(lineage,group)],'01_eligibility_summary.csv')
m<-as.data.table(readRDS(member_path));sets<-lapply(split(m$gene_symbol,m$pathway),unique);ck('MSigDB frozen version',all(m$db_version=='2026.1.Hs'))
saveRDS(list(patients=dt,raw=raw,y=yy,keep=keep,nuisance_design=X0,sets=sets,membership=m),file.path(out,'resources/outcome_prepared.rds'))
wt(rbindlist(checks),'01_preflight_checks.csv');capture.output(sessionInfo(),file=file.path(out,'logs/sessionInfo_prepare.txt'));cat('PREPARATION PASSED\n')
