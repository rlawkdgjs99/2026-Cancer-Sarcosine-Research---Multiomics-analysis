suppressPackageStartupMessages({library(data.table);library(Matrix);library(digest)})
arg<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);out<-normalizePath(file.path(dirname(arg),'..'));src<-dirname(out)
roots<-c(CD8='14_AllCell_Degradation_CD8_cDC1_26.09.08',cDC='17_cDC_AllCell_Degradation_26.09.16',Epithelial=basename(out))
mods<-lapply(names(roots),function(ct)readRDS(file.path(src,roots[ct],'resources/models',paste0(ct,'_post_group.rds'))));names(mods)<-names(roots)
frozen<-fread(file.path(src,roots['CD8'],'results/tables/01_FROZEN_allcell_patient_scores_groups.csv'))
checks<-list();add<-function(label,pass,detail=''){checks[[length(checks)+1L]]<<-data.table(check=label,pass=isTRUE(pass),detail=detail);stopifnot(pass)}
summary<-list()
for(ct in names(mods)){
 m<-mods[[ct]];dt<-as.data.table(m$patients);ref<-frozen[match(dt$Patient,Patient)]
 add(paste(ct,'frozen groups and scores'),identical(as.character(dt$group),as.character(ref$group))&&max(abs(dt$degradation_mean_z-ref$degradation_mean_z))<1e-12)
 add(paste(ct,'post only independent patient units'),all(dt$treatment=='Post')&&!anyDuplicated(dt$Patient)&&!anyDuplicated(dt$Sample))
 nc<-dt[[switch(ct,CD8='CD8_cells',cDC='cDC_cells',Epithelial='lineage_cells')]]
 add(paste(ct,'same >=50 lineage-cell rule'),all(nc>=50))
 add(paste(ct,'same histology group design'),identical(colnames(m$design),colnames(mods$CD8$design)))
 add(paste(ct,'positive High coefficient'),colnames(m$design)[m$coefficient]=='groupHigh')
 add(paste(ct,'4 exposure genes excluded'),!any(c('GNMT','DMGDH','SARDH','PIPOX')%in%names(m$rank)))
 summary[[ct]]<-data.table(compartment=ct,patients=nrow(dt),High=sum(dt$group=='High'),Low=sum(dt$group=='Low'),cells=sum(nc),ranked_genes=length(m$rank),min_cells=min(nc),model=paste(colnames(m$design),collapse=' + '),patient_ids=paste(dt$Patient,collapse=';'))
}
# Compare the exact calls used in source, independently of prose.
getcalls<-function(path,name){ans<-list();walk<-function(e){if(missing(e))return(invisible(NULL));if(is.call(e)){if(is.symbol(e[[1]])&&as.character(e[[1]])==name)ans[[length(ans)+1L]]<<-e;for(x in as.list(e)[-1])walk(x)}else if(is.expression(e))for(x in e)walk(x)};walk(parse(path));ans}
for(fn in c('voom','eBayes','fgseaMultilevel')){
 calls<-lapply(roots,function(r)getcalls(file.path(src,r,'scripts/02_run.R'),fn));add(paste('Exact source call parity:',fn),all(vapply(calls,function(x)length(x)==1,logical(1)))&&all(vapply(calls,function(x)identical(x[[1]],calls[[1]][[1]]),logical(1))))
}
# Source inputs remained unchanged since the original run.
man<-fread(file.path(out,'results/tables/00_input_manifest.csv'));add('Recorded analysis inputs unchanged',all(vapply(man$path,function(p)digest(file=p,algo='sha256',serialize=FALSE),character(1))==man$sha256))
fwrite(rbindlist(summary),file.path(out,'results/qa/method_consistency_cohorts.csv'));fwrite(rbindlist(checks),file.path(out,'results/qa/method_consistency_checks.csv'));print(rbindlist(summary)[,patient_ids:=NULL]);cat('All',length(checks),'consistency checks passed. BH family distinction remains explicitly documented.\n')
