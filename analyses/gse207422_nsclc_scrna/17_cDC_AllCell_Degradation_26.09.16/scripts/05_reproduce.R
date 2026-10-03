#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])),'common.R'))
m<-as.data.table(readRDS(member_path));original<-fread(file.path(out,'results/tables/02_GSEA_all_scopes.csv'));checks<-list()
for(si in seq_along(scopes)){sc<-scopes[si];b<-readRDS(file.path(out,'resources/models',paste0('cDC_',sc,'.rds')));rk<-b$rank
 for(ci in seq_along(collections)){co<-collections[ci];paths<-lapply(split(m[collection==co,gene_symbol],m[collection==co,pathway]),unique);n<-vapply(paths,function(g)sum(names(rk)%in%g),integer(1));paths<-paths[n>=15&n<=500];set.seed(26091600+si*100+ci)
  a<-fgseaMultilevel(paths,rk,minSize=15,maxSize=500,eps=0,scoreType='std',gseaParam=1,sampleSize=101,nPermSimple=10000,nproc=1,BPPARAM=BiocParallel::SerialParam(progressbar=FALSE));r<-original[scope==sc&collection==co];setorder(a,pathway);setorder(r,pathway)
  for(k in c('ES','NES','pval','log2err','size')){ok<-identical(a[[k]],r[[k]])||isTRUE(all.equal(a[[k]],r[[k]],tolerance=1e-13));checks[[length(checks)+1L]]<-data.table(scope=sc,collection=co,check=k,pass=ok);stopifnot(ok)}
  ok<-identical(vapply(a$leadingEdge,paste,collapse=';',FUN.VALUE=character(1)),r$leading_edge);checks[[length(checks)+1L]]<-data.table(scope=sc,collection=co,check='leading_edges',pass=ok);stopifnot(ok);cat('Reproduced',sc,co,'\n');flush.console()
 }
}
wt(rbindlist(checks),'05_fixed_seed_reproducibility.csv');capture.output(sessionInfo(),file=file.path(out,'logs/sessionInfo_reproduce.txt'));cat('PASS',length(checks),'reproducibility checks\n')
