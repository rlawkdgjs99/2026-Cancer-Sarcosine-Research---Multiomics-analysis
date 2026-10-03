options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly=TRUE); o<-normalizePath(if(length(args)) args[1] else '.'); b<-dirname(dirname(o));.libPaths(c(file.path(b,'R_libs'),.libPaths()))
library(data.table);library(IOBR)
options(IOBR.cache_dir=file.path(b,'reference_data/IOBR_v2.2.3_data-v1.0'))
m<-fread(file.path(o,'inputs/matched16_source.csv'))
sids<-c(m$sample_PRE,m$sample_EDT)
e<-fread(readLines(file.path(o,'inputs/expression_path.txt'),warn=FALSE),check.names=FALSE)
f<-as.matrix(e[,..sids]); storage.mode(f)<-'double';rownames(f)<-trimws(e[[1]])
stopifnot(ncol(f)==32,!anyDuplicated(rownames(f)),all(is.finite(f)),all(f>=0))
q<-deconvo_quantiseq(f,tumor=TRUE,arrays=FALSE,scale_mrna=TRUE)
stopifnot(nrow(q)==32,setequal(q$ID,sids))
fwrite(q,file.path(o,'tables/quanTIseq_matched32.csv'))
p<-fread(file.path(o,'inputs/prior_quanTIseq_PRE73.csv'))
qa<-merge(as.data.table(q),p,by='ID',suffixes=c('_new','_old'))
errs<-sapply(setdiff(names(p),'ID'),function(n) max(abs(qa[[paste0(n,'_new')]]-qa[[paste0(n,'_old')]])))
print(errs);stopifnot(max(errs)<1e-8)
cyt<-log2(sqrt((f['GZMA',]/colSums(f)*1e6)*(f['PRF1',]/colSums(f)*1e6)))
src<-fread(file.path(o,'inputs/CYT_all91.csv'))
stopifnot(max(abs(cyt-src$CYT_score[match(names(cyt),src$sample_id)]))<1e-10)
long<-rbindlist(lapply(c('PRE','EDT'),function(tp){
 ids<-m[[paste0('sample_',tp)]]
 z<-data.table(patient_id=m$patient_id,sample_id=ids,response=m$response_group,therapy=m$therapy,age=m$age,sex=m$sex,timepoint=tp)
 for(g in c('SARDH','PIPOX','GNMT','DMGDH')){
  vals<-log2(f[g,ids]+1);stopifnot(max(abs(vals-m[[paste0(tp,'_',g)]]))<1e-10);set(z,j=g,value=as.numeric(vals))
 }
 z[,CD8_fraction:=q$T_cells_CD8_quantiseq[match(ids,q$ID)]];z[,CD8_percent:=100*CD8_fraction];z[,CYT:=as.numeric(cyt[ids])];z
}))
# Source clinical table has an additional row-index column, consumed by read.delim row.names=1.
clin<-read.delim(file.path(b,'Melanoma-PRJEB23709_ClinicalData.tsv'),row.names=1,check.names=FALSE)
c<-clin[match(long$sample_id,clin$sample_id),]
stopifnot(identical(c$patient_name,long$patient_id),identical(c$Treatment,long$timepoint),identical(ifelse(c$response_NR=='R','R','NR'),long$response))
stopifnot(all(c$Gender==long$sex),all(c$age_start==long$age),all(ifelse(c$Therapy=='anti-PD-1','antiPD1','combo')==long$therapy))
fwrite(long,file.path(o,'tables/patient_data_long32.csv'))
write.csv(data.frame(feature=names(errs),max_abs_PRE_difference=errs),file.path(o,'qa/PRE_quanTIseq_reproduction.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(o,'qa/sessionInfo.txt'))
cat('PASS: gene32 raw expressions, original clinical pair IDs, CYT32 reproduction and PRE16 quanTIseq reproduction\n')
