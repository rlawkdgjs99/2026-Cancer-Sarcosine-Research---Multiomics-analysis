#!/usr/bin/env Rscript
# Display-only replacement candidate for Figure8a. Frozen UMAP, labels and
# patient-level degradation scores are retained; no new clustering or inference.
options(warn=1, device=function(...) grDevices::pdf(file=NULL, ...))
args <- commandArgs(FALSE)
script <- sub('^--file=','',args[grepl('^--file=',args)])
base <- normalizePath(file.path(dirname(script),'..'))
final <- dirname(base)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(grid);library(jsonlite);library(digest)})
inputs <- c(
 cells=file.path(final,'02_lineage_reannotation/results/tables/09_final_cell_lineages_FROZEN.csv'),
 patients=file.path(base,'results/tables/01_FROZEN_allcell_patient_scores_groups.csv'),
 composition=file.path(base,'results/tables/00_patient_cell_composition.csv'),
 reference=file.path(base,'results/tables/01_exposure_component_reference.csv'),
 lock=file.path(base,'results/tables/01_exposure_lock.csv')
)
sha <- function(p)digest(file=p,algo='sha256')
input_hashes <- vapply(inputs,sha,character(1))
cells <- fread(inputs['cells']); md <- fread(inputs['patients']);comp<-fread(inputs['composition']);ref<-fread(inputs['reference']);lock<-fread(inputs['lock'])
stopifnot(nrow(cells)==92053,!anyDuplicated(cells$cell_id),nrow(md)==15,!anyDuplicated(md$Patient),all(is.finite(cells$umap_1)),all(is.finite(cells$umap_2)))
stopifnot(sha(inputs['patients'])==lock[file=='results/tables/01_FROZEN_allcell_patient_scores_groups.csv',sha256])
cutoff<-unique(md$cutoff);stopifnot(length(cutoff)==1,abs(median(md$degradation_mean_z)-cutoff)<1e-14)
# Independent arithmetic check of existing scores using preserved raw counts and TMM library sizes.
for(g in c('SARDH','PIPOX')){
 xx<-log2(1+md[[paste0(g,'_raw')]]/md$effective_library*1e6)
 rr<-ref[gene==g];zz<-(xx-rr$mean_log)/rr$sd_log
 stopifnot(max(abs(xx-md[[paste0(g,'_log2CPM')]]))<1e-12,max(abs(zz-md[[paste0(g,'_z')]]))<1e-12)
}
stopifnot(max(abs((md$SARDH_z+md$PIPOX_z)/2-md$degradation_mean_z))<1e-12)
stopifnot(all(ifelse(md$degradation_mean_z>cutoff,'High','Low')==md$group))
keep <- cells[final_lineage!='Excluded residual doublet']
stopifnot(nrow(keep)==91844,sum(keep$final_lineage=='NK/gamma-delta T unresolved')==1332)
patient_n<-keep[,.N,by=Patient];stopifnot(all(patient_n$N==md$all_cells[match(patient_n$Patient,md$Patient)]))
pp<-md[treatment=='Post'][order(degradation_mean_z,Patient)]
stopifnot(nrow(pp)==12,sum(pp$group=='High')==7,sum(pp$group=='Low')==5)
plotcells<-keep[Patient%in%pp$Patient];stopifnot(nrow(plotcells)==78192,uniqueN(plotcells$Patient)==12)
cc<-plotcells[,.N,by=final_lineage]
cc[,fraction:=N/sum(N)]
refcomp<-comp[Patient%in%pp$Patient,.(N=sum(cells)),by=final_lineage]
stopifnot(all(cc$N==refcomp$N[match(cc$final_lineage,refcomp$final_lineage)]))
lineage_order<-c('Epithelial','CAF','B cell','Plasma cell','CD4 T cell','CD8 T cell','Cycling T cell','NK cell','Mast cell','Neutrophil','Monocyte','Macrophage','Conventional DC','pDC','NK/gamma-delta T unresolved')
stopifnot(setequal(cc$final_lineage,lineage_order))
# Original lineage palette, with a neutral colour for retained unresolved singlets.
lineage_cols<-c(setNames(hcl.colors(14,palette='Dark 3'),lineage_order[1:14]),'NK/gamma-delta T unresolved'='#929BA3')
lineage_labels<-setNames(lineage_order,lineage_order)
lineage_labels['CD4 T cell']<-'CD4⁺ T cells';lineage_labels['CD8 T cell']<-'CD8⁺ T cells';lineage_labels['NK/gamma-delta T unresolved']<-'Unresolved NK/γδ T'
plotcells[,final_lineage:=factor(final_lineage,levels=lineage_order)]
set.seed(260917)
plotcells<-plotcells[sample(.N)] # draw order only; all retained Post cells plotted
ink<-'#223843';muted<-'#637784';group_cols<-c('Low'='#2D8195','High'='#C44770')
p_umap<-ggplot(plotcells,aes(umap_1,umap_2,colour=final_lineage))+
 geom_point(size=.17,alpha=1,stroke=0)+
 scale_colour_manual(values=lineage_cols,labels=lineage_labels,drop=FALSE)+
 coord_equal(expand=FALSE)+theme_void(base_family='Arial',base_size=11)+
 labs(colour=NULL)+theme(legend.position='right',legend.text=element_text(size=10.4,colour=ink),
 legend.key.height=unit(3.55,'mm'),legend.key.width=unit(4,'mm'),legend.margin=margin(0,0,0,3),plot.margin=margin(5,3,3,2))+
 guides(colour=guide_legend(override.aes=list(size=2.8,alpha=1),ncol=1))
pp[,Patient_plot:=factor(Patient,levels=Patient)]
p_score<-ggplot(pp,aes(Patient_plot,degradation_mean_z,fill=group))+
 geom_col(width=.65)+geom_hline(yintercept=cutoff,linewidth=.45,linetype='dashed',colour='#667680')+
 scale_fill_manual(values=group_cols,guide='none')+
 scale_y_continuous(breaks=c(-.5,0,.5,1),limits=c(-.85,1.22),expand=expansion(mult=c(0,0)))+
 labs(x=NULL,y='Mean gene z score')+theme_classic(base_family='Arial',base_size=11)+
 theme(axis.text=element_text(colour=ink,size=10.6),axis.title.y=element_text(size=10.8,margin=margin(r=7)),
 axis.ticks=element_line(linewidth=.35,colour='#8797A0'),axis.line=element_line(linewidth=.4,colour='#8797A0'),
 axis.ticks.length=unit(1.2,'mm'),plot.margin=margin(3,10,3,6))

tx<-function(s,x,y,size=11,bold=FALSE,col=ink,just='left')grid.text(s,x,y,just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col,lineheight=1.16))
draw<-function(){
 grid.newpage();grid.rect(gp=gpar(fill='white',col=NA))
 tx('Tumour cell types and degradation groups',.025,.965,16,TRUE)
 tx('Post-treatment NSCLC · GSE207422 · 12 patients',.025,.921,10.8,col=muted)
 tx('Cell types',.035,.871,12.5,TRUE)
 tx('78,192 captured singlets',.96,.871,10.3,col=muted,just='right')
 pushViewport(viewport(x=.5,y=.647,width=.965,height=.415));grid.draw(ggplotGrob(p_umap));popViewport()
 tx('Whole-tumour sarcosine-degradation score',.035,.403,12.5,TRUE)
 tx('Low  n = 5',.035,.363,10.7,TRUE,col=group_cols['Low'])
 tx('High  n = 7',.29,.363,10.7,TRUE,col=group_cols['High'])
 tx('One bar per patient',.96,.363,10.4,col=muted,just='right')
 pushViewport(viewport(x=.5,y=.235,width=.965,height=.216));grid.draw(ggplotGrob(p_score));popViewport()
 tx('Score = mean of SARDH and PIPOX gene z scores in whole-tumour pseudobulk.',.025,.087,9.1,col=muted)
 tx('Dashed line: original 15-patient median (0.0021). High > median; Low ≤ median.',.025,.053,9.1,col=muted)
 tx('Fixed groups retained for cell-type analyses, with cell-type eligibility applied.',.025,.019,9.1,col=muted)
}
figdir<-file.path(base,'results/figures');tabdir<-file.path(base,'results/tables');qadir<-file.path(base,'results/qa')
stopifnot(dir.exists(figdir),dir.exists(tabdir),dir.exists(qadir))
fwrite(pp[,.(Patient,Sample,treatment,all_cells,SARDH_z,PIPOX_z,degradation_mean_z,cutoff,group)],file.path(tabdir,'Fig8a_Patient_Groups.csv'),scipen=0)
fwrite(cc[order(match(as.character(final_lineage),lineage_order))],file.path(tabdir,'Fig8a_Cell_Composition.csv'))
pdf<-file.path(figdir,'A.pdf');master<-file.path(figdir,'Fig8a_Cell_Context.png')
cairo_pdf(pdf,width=6.4,height=6.3,family='Arial',onefile=TRUE);draw();dev.off()
ragg::agg_png(master,width=6.4,height=6.3,units='in',res=600,background='white');draw();dev.off()
stopifnot(identical(vapply(inputs,sha,character(1)),input_hashes))
qa<-list(scope='Display only: original UMAP subset to primary Post patients and frozen patient scores/groups. No new clustering, score fit, median split or statistical test.',
 inputs=setNames(lapply(seq_along(inputs),function(i)list(path=inputs[i],sha256=input_hashes[i])),names(inputs)),
 patient_n=12,cell_n=78192,High=7,Low=5,original_median=unname(cutoff),all_reference_patients=15,
 score_validation='Gene log2(1+TMMCPM), original mean/sd z, mean-two-gene score and median labels agree within 1e-12.',
 cell_validation='All Post singlets, including1202 unresolved NK/gamma-delta cells, match frozen lineage and patient composition counts. 209 residual doublets excluded from all15reference cells.',
 UMAP='Frozen coordinates retained, no reduction/clustering recomputed. Draw-order seed260917 only.',
 q_or_p_values='None: group definition/display, not an inferential comparison.',
 RNA_proxy='Expression-derived score, not sarcosine concentration or flux.',
 palettes=list(lineage=as.list(lineage_cols),group=as.list(group_cols)),software=list(R=R.version.string,ggplot2=as.character(packageVersion('ggplot2')),ragg=as.character(packageVersion('ragg'))),
 transport_QA='pending',visual_QA='pending')
write_json(qa,file.path(qadir,'Fig8a_Cell_Context_QA.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat('Written:',master,'\n',pdf,'\n')
