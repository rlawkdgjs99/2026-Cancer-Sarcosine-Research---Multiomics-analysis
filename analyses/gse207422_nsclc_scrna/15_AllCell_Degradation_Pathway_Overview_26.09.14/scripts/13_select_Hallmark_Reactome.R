# Display-only selection from frozen scRNA patient-pseudobulk GSEA results.
# Same q/NES/Jaccard/cap rule as original GO:BP figures; collection changes only.
suppressPackageStartupMessages({library(data.table);library(digest);library(jsonlite)})
arg <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
out <- normalizePath(file.path(dirname(arg),'..'))
base <- dirname(out)
roots <- c(Whole_tumour=out,CD8=file.path(base,'14_AllCell_Degradation_CD8_cDC1_26.09.08'),cDC=file.path(base,'17_cDC_AllCell_Degradation_26.09.16'))
modelnames <- c(Whole_tumour='All_cells',CD8='CD8',cDC='cDC')
member <- file.path(base,'07_Fig6_Fig7_Three_Collection_Link_26.09.07/resources/MSigDB_2026.1.Hs_three_collections.rds')
tab <- file.path(out,'results/tables'); qa <- file.path(out,'results/qa')
heatfile <- file.path(tab,'Fig6_scRNA_GSEA_heatmap.csv')
sourcefiles <- c(member,heatfile,file.path(out,'results/figures/SC.png'),unlist(lapply(names(roots),function(k)c(file.path(roots[k],'resources/models',paste0(modelnames[k],'_post_group.rds')),file.path(roots[k],'results/tables/02_GSEA_all_scopes.csv')))))
sha <- function(p)digest(file=p,algo='sha256',serialize=FALSE)
hashes <- data.table(path=sourcefiles,sha256=vapply(sourcefiles,sha,character(1)))
members <- as.data.table(readRDS(member)); stopifnot(identical(unique(members$db_version),'2026.1.Hs'))
sets <- lapply(split(members$gene_symbol,members$pathway),unique)
heat <- fread(heatfile); focus <- unique(heat$pathway)
allrows <- selected <- counts <- rosters <- universes <- checks <- list()
for (ct in names(roots)) {
 model <- readRDS(file.path(roots[ct],'resources/models',paste0(modelnames[ct],'_post_group.rds')))
 rank <- model$rank; universe <- names(rank); patients <- as.data.table(model$patients)
 stopifnot(all(is.finite(rank)),!anyDuplicated(universe),all(diff(rank)<=0),!any(c('GNMT','DMGDH','SARDH','PIPOX')%in%universe),!anyDuplicated(patients$Patient))
 hh <- heat[compartment==ct]
 stopifnot(nrow(hh)==4,nrow(patients)==hh$n[1],sum(patients$group=='High')==hh$High[1],sum(patients$group=='Low')==hh$Low[1])
 rosters[[ct]] <- patients[,.(compartment=ct,Patient,group,histology)]
 universes[[ct]] <- data.table(compartment=ct,gene=universe)
 raw <- fread(file.path(roots[ct],'results/tables/02_GSEA_all_scopes.csv'))[scope=='post_group']
 qcol <- if(ct=='cDC') 'q_global_cDC' else 'q_global'
 family_n <- nrow(raw)
 stopifnot(family_n==hh$bh_family_n[1],max(abs(p.adjust(raw$pval,'BH')-raw[[qcol]]))<1e-12)
 if('cell_type'%in%names(raw))raw <- raw[cell_type==modelnames[ct]]
 stopifnot(!anyDuplicated(raw$pathway))
 for(id in focus){ref<-raw[pathway==id];h<-hh[pathway==id];stopifnot(nrow(ref)==1,ref$NES==h$NES,abs(ref[[qcol]]-h$BH_q)<1e-12*h$BH_q)}
 d <- raw[collection%in%c('Hallmark','Reactome'),.(pathway,collection,NES,ES,BH_q=get(qcol),gene_set_size=size,leading_edge_n,leading_edge,n,High,Low)]
 d[,':='(compartment=ct,bh_family_n=family_n,q_source_column=qcol,neglog10_BH_q=-log10(BH_q),absNES=abs(NES),display_order=NA_integer_,max_Jaccard_to_selected=NA_real_,overlap_with='',display_status=ifelse(NES<=0,'LOW_DIRECTION',ifelse(BH_q>=.05,'BH_Q_GE_0.05','CANDIDATE')),fixed_program=pathway%in%focus)]
 stopifnot(all(is.finite(d$NES)),all(d$BH_q>0&d$BH_q<=1),all(d$n==nrow(patients)),all(d$High==hh$High[1]),all(d$Low==hh$Low[1]))
 d[,measured_genes:=vapply(pathway,function(id)paste(sort(intersect(sets[[id]],universe)),collapse=';'),character(1))]
 stopifnot(all(lengths(strsplit(d$measured_genes,';',fixed=TRUE))==d$gene_set_size),all(lengths(strsplit(d$leading_edge,';',fixed=TRUE))==d$leading_edge_n))
 setorderv(d,c('BH_q','absNES','pathway'),c(1L,-1L,1L))
 chosen <- character()
 for(i in which(d$display_status=='CANDIDATE')) {
   id <- d$pathway[i]; a <- intersect(sets[[id]],universe)
   overlaps <- if(length(chosen)) vapply(chosen,function(prev){b<-intersect(sets[[prev]],universe);length(intersect(a,b))/length(union(a,b))},numeric(1)) else numeric()
   mx <- if(length(overlaps))max(overlaps) else 0
   reason <- if(length(chosen)>=8) 'DISPLAY_LIMIT_8' else if(mx>=.5) 'GENE_OVERLAP_GE_0.50' else 'SELECTED'
   d[i,':='(max_Jaccard_to_selected=mx,overlap_with=if(length(overlaps))names(which.max(overlaps))else'',display_status=reason)]
   if(reason=='SELECTED'){chosen<-c(chosen,id);d[i,display_order:=length(chosen)]}
 }
 sel <- d[display_status=='SELECTED'];stopifnot(nrow(sel)<=8,all(sel$NES>0),all(sel$BH_q<.05),all(sel$max_Jaccard_to_selected<.5))
 allrows[[ct]] <- d; selected[[ct]] <- sel
 counts[[ct]] <- d[,.(tested=.N,significant_High=sum(NES>0&BH_q<.05),displayed=sum(display_status=='SELECTED')),by=.(compartment,collection)]
 checks[[ct]] <- data.table(compartment=ct,n=nrow(patients),High=sum(patients$group=='High'),Low=sum(patients$group=='Low'),ranked_genes=length(universe),bh_family_n=family_n,source_BH_reproduced=TRUE,all_four_original_program_values_match=TRUE)
}
aa <- rbindlist(allrows); ss <- rbindlist(selected)
fwrite(aa,file.path(tab,'Fig6_HR_all_results_selection_audit.csv'),quote=TRUE,bom=TRUE)
fwrite(ss,file.path(tab,'Fig6_HR_selected_plotdata.csv'),quote=TRUE,bom=TRUE)
fwrite(aa[fixed_program==TRUE],file.path(tab,'Fig6_HR_four_program_status.csv'),quote=TRUE,bom=TRUE)
fwrite(rbindlist(counts),file.path(tab,'Fig6_HR_display_counts.csv'),quote=TRUE,bom=TRUE)
fwrite(rbindlist(rosters),file.path(qa,'HR_patient_rosters.csv'),quote=TRUE)
fwrite(rbindlist(universes),file.path(qa,'HR_ranked_genes.csv'),quote=TRUE)
fwrite(members[collection%in%c('Hallmark','Reactome'),.(pathway,gene_symbol,collection)],file.path(qa,'HR_membership.csv'),quote=TRUE)
stopifnot(identical(hashes$sha256,unname(vapply(hashes$path,sha,character(1)))))
fwrite(hashes,file.path(qa,'HR_input_hashes.csv'),quote=TRUE)
fwrite(rbindlist(checks),file.path(qa,'HR_build_checks.csv'),quote=TRUE)
write_json(list(type='display selection only; no GSEA rerun',collection=c('Hallmark','Reactome'),scope='post_group',direction='NES > 0',BH_q_threshold=0.05,ordering=c('BH q ascending','absolute NES descending','pathway ID ascending'),gene_overlap='Jaccard on gene-set members intersected with each original ranked universe; reject >= 0.50',maximum_per_compartment=8,fixed_program_priority=FALSE,BH_families=c(Whole_tumour=4958,CD8=5959,cDC=3797),database='MSigDB 2026.1.Hs'),file.path(qa,'HR_selection_rules.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,'logs/sessionInfo_HR_selection.txt'))
print(ss[,.(compartment,display_order,pathway,NES,BH_q,leading_edge_n)])
print(aa[fixed_program==TRUE,.(compartment,pathway,BH_q,display_status)])
