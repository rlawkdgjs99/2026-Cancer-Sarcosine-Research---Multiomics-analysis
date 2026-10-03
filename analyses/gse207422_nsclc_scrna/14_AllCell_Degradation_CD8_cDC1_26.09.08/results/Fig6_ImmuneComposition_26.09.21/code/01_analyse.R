#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]); base<-normalizePath(file.path(dirname(script),'..'))
plan<-fromJSON(file.path(base,'ANALYSIS_PLAN.json')); tabs<-file.path(base,'tables')
md<-fread(plan$inputs$patients$path); comp<-fread(plan$inputs$composition$path)
post<-md[treatment=='Post'][order(Patient)]
stopifnot(nrow(post)==12L,!anyDuplicated(post$Patient),sum(post$group=='High')==7L,sum(post$group=='Low')==5L,sum(post$all_cells)==78192L)
stopifnot(all(ifelse(md$degradation_mean_z>md$cutoff,'High','Low')==md$group))
post[,histology:=factor(histology)];post[,High:=as.integer(group=='High')]
counts<-comp[Patient %in% post$Patient];stopifnot(!anyDuplicated(counts[,.(Patient,final_lineage)]))
ct<-counts[,.(total=sum(cells)),by=Patient];stopifnot(all(ct$total[match(post$Patient,ct$Patient)]==post$all_cells))
stopifnot(setequal(setdiff(unique(counts$final_lineage),c('Epithelial','CAF')),plan$targets))
ids<-split(seq_len(nrow(post)),post$histology)
choices<-lapply(ids,function(ix)combn(ix,sum(post$High[ix]),simplify=FALSE))
perms<-do.call(rbind,unlist(lapply(choices[[1]],function(a)lapply(choices[[2]],function(z){h<-integer(nrow(post));h[c(a,z)]<-1L;h})),recursive=FALSE))
stopifnot(nrow(perms)==300L,!anyDuplicated(as.data.frame(perms)),sum(apply(perms,1,function(x)all(x==post$High)))==1L)
model<-function(y,h){m<-lm(y~histology+h,data=data.frame(y=y,histology=post$histology,h=h));c(beta=unname(coef(m)['h']),t=unname(coef(summary(m))['h','t value']))}
rows<-list();stats<-list();permutation_rows<-list()
for(k in seq_along(plan$targets)){
 cell<-plan$targets[k];num<-counts[final_lineage==cell][match(post$Patient,Patient),cells];den<-post$all_cells
 stopifnot(length(num)==12L,all(num>0),all(num<den),all(is.finite(num)))
 pct<-100*num/den;y<-log(num/(den-num));ob<-model(y,post$High)
 ts<-apply(perms,1,function(h)model(y,h)[2]);ne<-sum(abs(ts)>=abs(ob[2])-1e-12)
 lo<-pct[post$group=='Low'];hi<-pct[post$group=='High']
 rows[[k]]<-data.table(Patient=post$Patient,Sample=post$Sample,histology=as.character(post$histology),group=post$group,degradation_mean_z=post$degradation_mean_z,cell_type=cell,target_cells=num,denominator_cells=den,percent=pct,log_ratio=y,row_z=as.numeric(scale(pct)))
 stats[[k]]<-data.table(cell_type=cell,n_Low=5L,n_High=7L,mean_Low_percent=mean(lo),mean_High_percent=mean(hi),SD_Low_percent=sd(lo),SD_High_percent=sd(hi),difference_pp=mean(hi)-mean(lo),log_ratio_beta=ob[1],t_statistic=ob[2],extreme_assignments=ne,total_assignments=300L,raw_P=ne/300)
 permutation_rows[[k]]<-data.table(cell_type=cell,assignment=seq_along(ts),t_statistic=ts)
}
res<-rbindlist(stats);res[,BH13_q:=p.adjust(raw_P,'BH')]
values<-rbindlist(rows);stopifnot(nrow(values)==156L,all(is.finite(values$row_z)))
old<-fread(plan$inputs$old_stats$path)[scope=='all_singlets']
for(cell in old$cell_type){stopifnot(abs(res[cell_type==cell,raw_P]-old[cell_type==cell,raw_P])<1e-12,abs(res[cell_type==cell,difference_pp]-old[cell_type==cell,raw_mean_difference_pp])<1e-10)}
fwrite(values,file.path(tabs,'patient_values.csv'));fwrite(res,file.path(tabs,'statistics.csv'));fwrite(rbindlist(permutation_rows),file.path(tabs,'permutation_statistics.csv'));fwrite(post,file.path(tabs,'frozen_patient_metadata.csv'))
for(value in c('percent','row_z','target_cells'))fwrite(dcast(values,cell_type~Patient,value.var=value),file.path(tabs,paste0(value,'_matrix.csv')))
pd<-as.data.table(perms);setnames(pd,post$Patient);fwrite(pd,file.path(tabs,'permutation_assignments.csv'))
sink(file.path(base,'qa/sessionInfo_R.txt'));print(sessionInfo());sink()
print(res[,.(cell_type,mean_Low_percent,mean_High_percent,difference_pp,raw_P,BH13_q)])
