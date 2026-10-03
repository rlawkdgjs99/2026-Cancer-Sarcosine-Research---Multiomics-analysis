#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])),"common.R"))
checks<-list();record<-function(k,ok,detail=""){checks[[length(checks)+1L]]<<-data.table(check=k,passed=isTRUE(ok),detail=as.character(detail));if(!isTRUE(ok))stop(k,": ",detail)}
aa<-fread(file.path(out,"results/tables/02_GSEA_all_scopes.csv"));ca<-fread(file.path(out,"results/tables/03_CAMERA_primary.csv"))
ss<-fread(file.path(out,"results/tables/04_exact_scores_HC3.csv"));prep<-readRDS(file.path(out,"resources/outcome_prepared.rds"));pbo<-readRDS(file.path(out,"resources/allcell_prepared.rds"))
md<-pbo$metadata;m<-as.data.table(readRDS(member_path));sets<-lapply(split(m$gene_symbol,m$pathway),unique)
manifest<-fread(file.path(out,"results/tables/00_inputs.csv"));manifest[,sha_after:=vapply(file.path(src,path),sha,character(1))];manifest[,unchanged:=sha256==sha_after]
record("Input checksums unchanged",all(manifest$unchanged));wt(manifest,"06_input_preservation.csv")
prot<-fread(file.path(out,"results/tables/00_protected_before.csv"));prot[,sha_after:=vapply(file.path(src,path),sha,character(1))];prot[,unchanged:=sha256==sha_after]
record("Previous09_to13 artifacts unchanged",all(prot$unchanged),nrow(prot));wt(prot,"06_previous_artifact_preservation.csv")
lock<-fread(file.path(out,"results/tables/01_exposure_lock.csv"));record("Frozen plan scores and preparation unchanged",all(vapply(file.path(out,lock$file),sha,character(1))==lock$sha256))
eff<-pbo$exposure$y$samples$lib.size*pbo$exposure$y$samples$norm.factors
xx<-log2(1+sweep(pbo$counts$All_cells[c("SARDH","PIPOX"),],2,eff,"/")*1e6)
zx<-sweep(sweep(xx,1,rowMeans(xx),"-"),1,apply(xx,1,sd),"/");refD<-colMeans(zx)
record("Independent all-cell D",max(abs(refD-md$degradation_mean_z))<1e-12)
record("Independent global median groups",all(ifelse(refD>median(refD),"High","Low")==md$group))
record("All15 7High 8Low",nrow(md)==15&&sum(md$group=="High")==7&&sum(md$group=="Low")==8)
record("Post12 fixed 7High 5Low",nrow(md[treatment=="Post"])==12&&sum(md$group=="High"&md$treatment=="Post")==7&&sum(md$group=="Low"&md$treatment=="Post")==5)
cc<-fread(file.path(out,"results/tables/00_patient_cell_composition.csv"))
record("Patient composition totals match retained singlets",all(cc[,.(cells=sum(cells)),by=Patient][match(md$Patient,Patient),cells]==md$all_cells)&&sum(md$all_cells)==91844)
bh<-function(p){r<-rep(NA_real_,length(p));ok<-which(is.finite(p));ii<-ok[order(p[ok])];r[ii]<-pmin(1,rev(cummin(rev(p[ii]*length(p)/seq_along(ii)))));r}
err<-0
for(sc in unique(aa$scope)){
 d<-aa[scope==sc];record(paste(sc,"GSEA BH missingness"),identical(is.na(d$q_global),is.na(bh(d$pval))))
 err<-max(err,max(abs(d$q_global-bh(d$pval)),na.rm=TRUE))
 for(ct in unique(d$cell_type)){
  e<-d[cell_type==ct];err<-max(err,max(abs(e$q_cell_all_collections-bh(e$pval)),na.rm=TRUE))
  for(co in unique(e$collection)){t<-e[collection==co];err<-max(err,max(abs(t$q_collection_cell-bh(t$pval)),na.rm=TRUE))}
 }
}
err<-max(err,max(abs(ca$q_global-bh(ca$PValue)),na.rm=TRUE))
for(sc in unique(ss$scope)){d<-ss[scope==sc];ref<-bh(d$p);record(paste(sc,"score BH missingness"),identical(is.na(ref),is.na(d$q_eight_scores)));if(any(is.finite(ref)))err<-max(err,max(abs(ref-d$q_eight_scores),na.rm=TRUE))}
record("Independent BH all families",err<1e-12,err)
record("Unique GSEA hypotheses",!anyDuplicated(aa[,.(scope,cell_type,collection,pathway)]))
record("GSEA P range",all(is.na(aa$pval)|(aa$pval>=0&aa$pval<=1)))
scoreerr<-0
for(ct in names(prep)){
 nm<-prep[[ct]];dt<-nm$patients;X0<-model.matrix(reformulate(covars(dt)),dt)
 h<-rowSums((X0%*%solve(crossprod(X0)))*X0);minn<-1/max(h);if(minn>10)minn<-10+(minn-10)*.7
 medlib<-median(colSums(nm$raw));cut<-10/medlib*1e6
 manualkeep<-rowSums(sweep(nm$raw,2,colSums(nm$raw),"/")*1e6>=cut)>=minn-1e-14 & rowSums(nm$raw)>=15-1e-14
 record(paste(ct,"independent nuisance filter"),identical(unname(manualkeep),unname(nm$keep)),sum(nm$keep))
 record(paste(ct,"outcome library sizes"),all(dt$outcome_retained_library==colSums(nm$y$counts))&&max(abs(dt$outcome_effective_library-dt$outcome_retained_library*dt$outcome_TMM_factor))<1e-8)
 raw<-pbo$counts[[ct]][,dt$Sample,drop=FALSE];record(paste(ct,"original outcome count alignment"),all(raw==nm$raw))
 ee<-log2(1+sweep(nm$y$counts,2,dt$outcome_effective_library,"/")*1e6)
 ee<-ee[!rownames(ee)%in%exposure_genes&apply(ee,1,sd)>0,,drop=FALSE]
 ez<-sweep(sweep(ee,1,rowMeans(ee),"-"),1,apply(ee,1,sd),"/")
 for(id in focus$pathway[1:4]){
  genes<-intersect(sets[[id]],rownames(ez))
  if(length(genes)<15){record(paste(ct,id,"score NE gene gate"),all(is.na(nm$scores[,id])));next}
  val<-colMeans(ez[genes,,drop=FALSE]);val<-(val-mean(val))/sd(val)
  scoreerr<-max(scoreerr,max(abs(val-nm$scores[,id])))
 }
}
record("Independent patient score reconstruction",scoreerr<1e-10,scoreerr)
curve<-list();eschecks<-list();hcerr<-0;maxwls<-0
for(f in sort(list.files(file.path(out,"resources/models"),full.names=TRUE,pattern="\\.rds$"))){
 b<-readRDS(f);sc<-b$scope;ct<-b$cell_type;X<-b$design;j<-b$coefficient;dt<-b$patients;nm<-prep[[ct]];ii<-match(dt$Patient,nm$patients$Patient);tag<-paste(ct,sc)
 record(paste(tag,"sample design alignment"),identical(dt$Patient,rownames(X))&&identical(dt$Sample,colnames(b$voom$E))&&!anyDuplicated(dt$Patient))
 record(paste(tag,"rank df group gate"),qr(X)$rank==ncol(X)&&nrow(X)-ncol(X)>=5&&min(table(dt$group))>=3)
 record(paste(tag,"unchanged allcell patient exposure"),all(as.character(dt$group)==as.character(md$group[match(dt$Patient,md$Patient)]))&&all(dt$degradation_mean_z==md$degradation_mean_z[match(dt$Patient,md$Patient)]))
 E<-b$voom$E;W<-b$voom$weights;del<-0
 record(paste(tag,"voom transformation"),max(abs(E-log2(sweep(as.matrix(b$y$counts)+.5,2,dt$outcome_effective_library+1,"/")*1e6)))<1e-10)
 for(i in seq_len(nrow(E))){
  inv<-solve(crossprod(X,X*W[i,]));bt<-as.vector(inv%*%crossprod(X,E[i,]*W[i,]))
  tt<-bt[j]/sqrt(inv[j,j]*b$fit$s2.post[i]);del<-max(del,max(abs(bt-b$fit$coefficients[i,])),abs(tt-b$fit$t[i,j]))
 }
 record(paste(tag,"independent all-gene WLS and moderated t"),del<1e-8,del);maxwls<-max(maxwls,del)
 rank<-b$rank;record(paste(tag,"rank gene contract"),!any(exposure_genes%in%names(rank))&&all(diff(rank)<=0)&&!anyDuplicated(names(rank)))
 sz<-vapply(sets,function(g)sum(names(rank)%in%g),integer(1));eligible<-names(sz)[sz>=15&sz<=500]
 data<-aa[cell_type==ct&scope==sc]
 record(paste(tag,"complete set universe"),setequal(data$pathway,eligible)&&all(data$size==sz[data$pathway]),length(eligible))
 for(id in focus$pathway){
  if(!(id%in%eligible))next
  hit<-names(rank)%in%sets[[id]];scalev<-2^30/sum(abs(rank));if(scalev>=1)scalev<-floor(scalev);w<-round(abs(rank)*scalev)
  walk<-cumsum(ifelse(hit,w/sum(w[hit]),-1/sum(!hit)))
  es<-if(abs(max(walk))>abs(min(walk)))max(walk)else min(walk);peak<-if(es>0)which.max(walk)else which.min(walk)
  le<-if(es>0)names(rank)[hit&seq_along(rank)<=peak]else names(rank)[hit&seq_along(rank)>peak]
  row<-data[pathway==id]
  eschecks[[paste(ct,sc,id)]]<-data.table(cell_type=ct,scope=sc,pathway=id,error=abs(es-row$ES),LE_equal=setequal(le,strsplit(row$leading_edge,";",fixed=TRUE)[[1]]))
  if(sc=="post_group"&&id%in%focus$pathway[1:4])curve[[paste(ct,id)]]<-data.table(cell_type=ct,pathway=id,rank=c(0L,seq_along(rank)),running_ES=c(0,walk),gene=c("",names(rank)),hit=c(FALSE,hit))
 }
 for(id in focus$pathway[1:4]){
  row<-ss[cell_type==ct&scope==sc&pathway==id];val<-nm$scores[ii,id]
  if(anyNA(val)){record(paste(tag,id,"HC3 NE low coverage"),row$status=="NE_GENE_COVERAGE_OR_SCORE");next}
  if(sd(val)<1e-12){record(paste(tag,id,"HC3 NE constant"),row$status=="NE_CONSTANT_SCORE");next}
  dd<-copy(dt);dd[,response:=val];pred<-if(endsWith(sc,"continuous"))"exposure_z"else"group"
  ff<-lm(reformulate(c(covars(dd),pred),response="response"),data=dd);h<-hatvalues(ff)
  record(paste(tag,id,"HC3 leverage eligibility"),(row$status=="TESTED")==all(1-h>=1e-8))
  if(row$status!="TESTED")next
  xm<-model.matrix(ff);inv<-solve(crossprod(xm));meat<-matrix(0,ncol(xm),ncol(xm))
  for(k in seq_len(nrow(xm)))meat<-meat+tcrossprod(xm[k,])*residuals(ff)[k]^2/(1-h[k])^2
  se<-sqrt((inv%*%meat%*%inv)[j,j]);bt<-coef(ff)[j];pv<-2*pt(-abs(bt/se),df.residual(ff))
  ci<-bt+c(-1,1)*qt(.975,df.residual(ff))*se
  hcerr<-max(hcerr,abs(bt-row$beta),abs(se-row$se),abs(pv-row$p),max(abs(ci-c(row$ci_low,row$ci_high))))
 }
 cat("Verified",tag,"\n")
}
es<-rbindlist(eschecks);record("Independent focus ES and leading edges",all(es$LE_equal)&&max(es$error)<1e-9,max(es$error))
record("Independent HC3 beta SE P CI",hcerr<1e-10,hcerr)
fv<-fread(file.path(out,"results/tables/05_focus_complete.csv"))
record("NE not assigned a P or q",all(is.na(fv[status!="TESTED",q_global]))&&all(is.na(fv[status!="TESTED",pval])))
pf<-fread(file.path(out,"results/tables/01_PREFIT_eligibility.csv"));fin<-fread(file.path(out,"results/tables/01_model_eligibility.csv"))
ii<-match(paste(fin$cell_type,fin$scope),paste(pf$cell_type,pf$scope))
record("Eligibility preflight equals fitted scopes",all(fin$n==pf$n[ii])&&all(fin$High==pf$High[ii])&&all(fin$Low==pf$Low[ii])&&all(fin$rank==pf$rank[ii]))
wt(es,"06_ES_checks.csv");wt(rbindlist(curve),"06_primary_running_ES.csv");wt(rbindlist(checks),"06_verification_checks.csv")
writeLines(c("PASS",paste("checks",length(checks)),paste("BH max error",err),paste("score max error",scoreerr),paste("all-gene WLS/t max error",maxwls),paste("ES max error",max(es$error)),paste("HC3 max error",hcerr),paste("prior artifacts unchanged",nrow(prot))),file.path(out,"results/qa/verification_summary.txt"))
cat("INDEPENDENT VERIFICATION PASS:",length(checks),"checks\n")
