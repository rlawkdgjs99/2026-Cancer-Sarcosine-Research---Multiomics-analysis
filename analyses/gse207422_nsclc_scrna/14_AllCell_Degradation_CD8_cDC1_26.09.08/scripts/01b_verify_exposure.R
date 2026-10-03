#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])),"common.R"))
a<-readRDS(file.path(out,"resources/allcell_prepared.rds"));m<-fread(file.path(out,"results/tables/01_FROZEN_allcell_patient_scores_groups.csv"))
checks<-list();add<-function(n,b){checks[[n]]<<-data.table(check=n,pass=isTRUE(b));stopifnot(isTRUE(b))}
eff<-a$exposure$y$samples$lib.size*a$exposure$y$samples$norm.factors
xx<-t(log2(1+t(a$counts$All_cells[c("SARDH","PIPOX"),])/eff*1e6))
zz<-(xx-rowMeans(xx))/apply(xx,1,sd);dd<-apply(zz,2,mean)
add("independent_all15_formula",max(abs(dd-m$degradation_mean_z))<1e-12)
add("one_median_groups",identical(ifelse(dd>median(dd),"High","Low"),setNames(m$group,colnames(xx))))
add("cutoff_fixed",max(abs(m$cutoff-median(dd)))<1e-12)
add("15_independent_samples",nrow(m)==15&&!anyDuplicated(m$Patient)&&!anyDuplicated(m$Sample))
add("both_components_detected",all(m$SARDH_raw>0)&all(m$PIPOX_raw>0))
add("all_cell_count",sum(m$all_cells)==91844)
add("CD8_cell_count",sum(m$CD8_cells)==17614)
add("cDC1_cell_count",sum(m$cDC1_cells)==129)
model_rows<-list()
for(ct in c("CD8","cDC1")){
 d<-copy(m[get(paste0(ct,"_cells"))>=if(ct=="CD8")50L else 5L]);d[,group:=factor(group,levels=c("Low","High"))]
 for(sc in c("post_group","post_continuous","all_group","all_continuous")){
  s<-copy(d);if(startsWith(sc,"post"))s<-s[treatment=="Post"]
  pred<-if(endsWith(sc,"continuous"))"exposure_z" else "group"
  X<-model.matrix(reformulate(c(covars(s),pred)),s)
  ok<-nrow(s)>=8 && min(table(s$group))>=3 && qr(X)$rank==ncol(X)&&nrow(s)-ncol(X)>=5 && (ct!="CD8"||nrow(d)>=10)
  model_rows[[paste(ct,sc)]]<-data.table(cell_type=ct,scope=sc,n=nrow(s),n_cells=sum(s[[paste0(ct,"_cells")]]),High=sum(s$group=="High"),Low=sum(s$group=="Low"),rank=qr(X)$rank,df=nrow(s)-qr(X)$rank,status=if(ok)"ELIGIBLE" else "NE",max_leverage=if(qr(X)$rank==ncol(X))max(rowSums((X%*%solve(crossprod(X)))*X))else NA_real_)
 }
}
wt(rbindlist(checks),"01_independent_exposure_checks.csv");wt(rbindlist(model_rows),"01_PREFIT_eligibility.csv")
lock<-fread(file.path(out,"results/tables/01_exposure_lock.csv"));stopifnot(all(vapply(file.path(out,lock$file),sha,character(1))==lock$sha256))
print(rbindlist(model_rows));cat("Exposure checks passed; no outcome fitted.\n")
