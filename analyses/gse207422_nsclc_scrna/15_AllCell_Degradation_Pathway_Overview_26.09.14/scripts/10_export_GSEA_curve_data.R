# Reconstruct exact GSEA walks from frozen ranks. No model or permutation rerun.
suppressPackageStartupMessages({library(data.table);library(digest);library(jsonlite)})
arg<-sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
out<-normalizePath(file.path(dirname(arg),'..'));base<-normalizePath(file.path(out,'..'))
.libPaths(c(.libPaths(),file.path(base,'03_lee_fig3_style_FINAL/hallmark_GSEA_Q4_vs_Q1_26.08.26/R_libs')))
suppressPackageStartupMessages(library(fgsea))
sha<-function(p)digest(file=p,algo='sha256',serialize=FALSE)
roots<-c(Whole_tumour=out,CD8=file.path(base,'14_AllCell_Degradation_CD8_cDC1_26.09.08'),cDC=file.path(base,'17_cDC_AllCell_Degradation_26.09.16'))
models<-c(Whole_tumour='All_cells',CD8='CD8',cDC='cDC')
member<-file.path(base,'07_Fig6_Fig7_Three_Collection_Link_26.09.07/resources/MSigDB_2026.1.Hs_three_collections.rds')
heatfile<-file.path(out,'results/tables/Fig6_scRNA_GSEA_heatmap.csv')
protected<-c(member,heatfile,file.path(out,'results/figures/SC.png'),
 unlist(lapply(names(roots),function(k)c(file.path(roots[k],'resources/models',paste0(models[k],'_post_group.rds')),file.path(roots[k],'results/tables/02_GSEA_all_scopes.csv')))))
before<-data.table(path=protected,sha256=vapply(protected,sha,character(1)))
h<-fread(heatfile);focus<-unique(h$pathway);stopifnot(length(focus)==4,nrow(h)==12)
m<-as.data.table(readRDS(member));sets<-lapply(split(m$gene_symbol,m$pathway),unique)
curves<-hits<-ranks<-patients<-checks<-annotations<-list()
for(ct in names(roots)){
 b<-readRDS(file.path(roots[ct],'resources/models',paste0(models[ct],'_post_group.rds')))
 rk<-b$rank;stopifnot(all(is.finite(rk)),all(diff(rk)<=0),!anyDuplicated(names(rk)),!any(c('GNMT','DMGDH','SARDH','PIPOX')%in%names(rk)))
 p<-as.data.table(b$patients);hh<-h[compartment==ct];stopifnot(nrow(p)==hh$n[1],sum(p$group=='High')==hh$High[1],sum(p$group=='Low')==hh$Low[1])
 patients[[ct]]<-p[,.(compartment=ct,Patient,group,histology)]
 ranks[[ct]]<-data.table(compartment=ct,rank=seq_along(rk),gene=names(rk),statistic=unname(rk))
 allg<-fread(file.path(roots[ct],'results/tables/02_GSEA_all_scopes.csv'))
 g<-allg[scope=='post_group'];if('cell_type'%in%names(g))g<-g[cell_type==models[ct]]
 # Preserve the archived fgsea 1.38.0 prepareStats integerized weighting kernel.
 factor<-2^30/sum(abs(rk));if(factor>=1)factor<-floor(factor)
 weights<-round(abs(rk)*factor)
 for(id in focus){
  ref<-g[pathway==id];label<-hh[pathway==id];qcol<-if(ct=='cDC')'q_global_cDC'else 'q_global'
  stopifnot(nrow(ref)==1,nrow(label)==1,identical(ref$NES,label$NES),abs(ref[[qcol]]-label$BH_q)<=1e-12*abs(label$BH_q))
  hit<-names(rk)%in%sets[[id]];idx<-which(hit);N<-length(rk);K<-sum(hit)
  stopifnot(K==ref$size,K>=15,K<=500)
  walk<-c(0,cumsum(ifelse(hit,weights/sum(weights[hit]),-1/(N-K))))
  es<-if(max(walk)> -min(walk))max(walk)else min(walk)
  peak<-if(es>0)which.max(walk)-1L else which.min(walk)-1L
  le<-if(es>0)names(rk)[hit&seq_along(rk)<=peak]else names(rk)[hit&seq_along(rk)>peak]
  native<-fgsea::calcGseaStat(stats=setNames(sign(rk)*weights,names(rk)),selectedStats=idx,gseaParam=1,scoreType='std')
  stopifnot(abs(es-ref$ES)<1e-10,abs(es-native)<1e-10,abs(tail(walk,1))<1e-10,
    setequal(le,strsplit(ref$leading_edge,';',fixed=TRUE)[[1]]),length(le)==ref$leading_edge_n)
  # All hit jumps and adjacent vertices retained; omitted miss-only runs are linear.
  nodes<-sort(unique(c(0L,N,idx,idx-1L)))
  curves[[paste(ct,id)]]<-data.table(compartment=ct,pathway=id,rank=nodes,rank_percent=nodes/N*100,running_ES=walk[nodes+1L])
  hits[[paste(ct,id)]]<-data.table(compartment=ct,pathway=id,rank=idx,rank_percent=idx/N*100,gene=names(rk)[idx],leading_edge=names(rk)[idx]%in%le)
  annotations[[paste(ct,id)]]<-data.table(compartment=ct,pathway=id,program=label$program,collection=ref$collection,n=label$n,High=label$High,Low=label$Low,
    ranked_genes=N,gene_set_size=K,ES=ref$ES,NES=ref$NES,BH_q=label$BH_q,bh_family_n=label$bh_family_n,
    leading_edge_n=ref$leading_edge_n,peak_rank=peak,weight_scale=factor)
  checks[[paste(ct,id)]]<-data.table(compartment=ct,pathway=id,ES_error=abs(es-ref$ES),fgsea_kernel_error=abs(es-native),leading_edge_match=TRUE,
    source_NES_q_match=TRUE,endpoint_error=abs(tail(walk,1)))
 }
}
dest<-file.path(out,'results/tables')
fwrite(rbindlist(curves),file.path(dest,'Fig6_exact_GSEA_curve_vertices.csv'),quote=TRUE)
fwrite(rbindlist(hits),file.path(dest,'Fig6_exact_GSEA_gene_hits.csv'),quote=TRUE)
fwrite(rbindlist(annotations),file.path(dest,'Fig6_exact_GSEA_curve_annotations.csv'),quote=TRUE)
qa<-file.path(out,'results/qa')
fwrite(rbindlist(ranks),file.path(qa,'GSEA_curve_frozen_ranks.csv'),quote=TRUE)
fwrite(rbindlist(patients),file.path(qa,'GSEA_curve_patient_rosters.csv'),quote=TRUE)
fwrite(rbindlist(checks),file.path(qa,'GSEA_curve_ES_checks.csv'),quote=TRUE)
stopifnot(identical(before$sha256,unname(vapply(before$path,sha,character(1)))))
fwrite(before,file.path(qa,'GSEA_curve_input_hashes.csv'),quote=TRUE)
capture.output(sessionInfo(),file=file.path(out,'logs/sessionInfo_GSEA_curves.txt'))
print(rbindlist(annotations)[,.(compartment,program,ES,NES,BH_q,gene_set_size)])
cat('12 original contrasts reconstructed and verified without inference rerun.\n')
