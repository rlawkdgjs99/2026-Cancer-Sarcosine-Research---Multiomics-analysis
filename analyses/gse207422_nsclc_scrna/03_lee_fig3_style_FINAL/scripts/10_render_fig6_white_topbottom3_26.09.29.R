#!/usr/bin/env Rscript
# Figure 6 candidates: frozen UMAP; display-only, no inferential tests.
a <- commandArgs(FALSE); script <- sub("^--file=", "", a[grepl("^--file=", a)])
b <- normalizePath(file.path(dirname(script), ".."))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(scales);library(ragg);library(jsonlite)})
out <- file.path(b,"results/figures_publication/Fig6_UMAP_Redraw_26.09.29")
for(x in c("figures","tables","qa"))dir.create(file.path(out,x),recursive=TRUE,showWarnings=FALSE)
s <- as.data.table(readRDS(file.path(b,"intermediate/01_cell_paper_style_scores.rds")))
lin <- fread(file.path(b,"../02_lineage_reannotation/results/tables/09_final_cell_lineages_FROZEN.csv"))
caps <- fread(file.path(b,"results/tables/plotdata_01_UMAP_display_caps.csv"))
meta_path <- file.path(b,"../27_Fig6_Current_Degradation_CSV_Export_26.09.23/staging/scRNA15_Metadata_Sarcosine.csv")
m <- fread(meta_path)
stopifnot(nrow(m)==15,!anyDuplicated(m$Patient),nrow(s)==92053,!anyDuplicated(s$cell_id))
i <- match(s$cell_id,lin$cell_id)
stopifnot(!anyNA(i),identical(s$analysis_eligible,lin$analysis_eligible[i]),identical(s$final_lineage,lin$final_lineage[i]),identical(s$umap_1,lin$umap_1[i]),identical(s$umap_2,lin$umap_2[i]))
# Independently reconstruct the stored 15-reference production score, sample SD z.
pcheck <- (as.numeric(scale(m$GNMT_log2CPM_ref15))+as.numeric(scale(m$DMGDH_log2CPM_ref15)))/2
stopifnot(max(abs(pcheck-m$Production_score_ref15))<1e-12)
setorder(m,-Production_score_ref15,Patient)
m[,production_rank_desc:=seq_len(.N)]
stopifnot(m$Production_score_ref15[3]>m$Production_score_ref15[4],m$Production_score_ref15[12]>m$Production_score_ref15[13])
m[,selection:=fifelse(production_rank_desc<=3,"High",fifelse(production_rank_desc>=13,"Low","Not selected"))]
d <- copy(s[analysis_eligible==TRUE]);stopifnot(nrow(d)==90512,uniqueN(d$Patient)==15)
d[,selection:=m$selection[match(Patient,m$Patient)]]
stopifnot(!anyNA(d$selection))
levels <- c("Epithelial","CAF","B cell","Plasma cell","CD4 T cell","CD8 T cell","Cycling T cell","NK cell","Mast cell","Neutrophil","Monocyte","Macrophage","Conventional DC","pDC")
cols <- setNames(c("#E64B6A","#CB6132","#AF8100","#817A00","#75AE00","#109F48","#00A88B","#0098B8","#2674D9","#5145BD","#9570D5","#CA65D7","#AC288F","#EE77A8"),levels)
d[,final_lineage:=factor(final_lineage,levels=levels)];stopifnot(!anyNA(d$final_lineage))
xr <- range(d$umap_1);yr <- range(d$umap_2)
th <- theme_void(base_family="Arial",base_size=12)+theme(plot.title=element_text(size=17,face="bold",colour="#193340",margin=margin(b=5)),plot.subtitle=element_text(size=10.5,colour="#586D77",lineheight=1.1,margin=margin(b=8)),plot.margin=margin(10,10,10,10),plot.background=element_rect(fill="white",colour=NA),legend.title=element_text(size=11),legend.text=element_text(size=10),legend.position="bottom")
set.seed(260929)
sel <- d[selection!="Not selected"];sel <- sel[sample(.N)]
lineage_plot <- function(z,title,subtitle){ggplot(z,aes(umap_1,umap_2,colour=final_lineage))+geom_point(size=.28,stroke=0,alpha=1)+scale_colour_manual(values=cols,drop=FALSE)+coord_equal(xlim=xr,ylim=yr)+th+labs(title=title,subtitle=subtitle,colour=NULL)+guides(colour=guide_legend(ncol=7,byrow=TRUE,override.aes=list(size=2.6,alpha=1)))+theme(legend.key.width=grid::unit(3.5,"mm"),legend.key.height=grid::unit(5,"mm"))}
pp <- list()
for(gr in c("Low","High")){
 z<-sel[selection==gr];ids<-m[selection==gr,Patient]
 pp[[gr]]<-lineage_plot(z,paste("Production",gr),paste0(paste(ids,collapse=" / "),"  |  ",format(nrow(z),big.mark=",")," cells"))
}
D <- (pp$Low|pp$High)+plot_layout(guides="collect")+plot_annotation(title="Cell lineages by patient production score",caption="Bottom 3 vs top 3 of 15 patients; whole-tumour pseudobulk mean z(GNMT, DMGDH).\nFrozen UMAP coordinates; Low: 1 pretreatment + 2 post-treatment; High: 3 post-treatment.",theme=theme(plot.title=element_text(size=19,face="bold",family="Arial",colour="#193340"),plot.caption=element_text(size=9,family="Arial",colour="#586D77",hjust=0))) & theme(legend.position="bottom")
L <- lineage_plot(d[sample(.N)],"Cell lineages","90,512 annotated cells | 15 patients")+guides(colour=guide_legend(ncol=2,byrow=TRUE,override.aes=list(size=2.8,alpha=1)))
pal <- list(Production=c("#FFFFFF","#FFF3DD","#FFBF55","#E56A0A","#943100"),Degradation=c("#FFFFFF","#E4F4FF","#65B9EB","#1776C6","#073B8C"))
plots <- list();audits<-list()
for(ax in c("Production","Degradation")){
 col<-paste0(tolower(ax),"_module_shifted");value<-d[[col]];cap<-caps[axis==ax,display_cap_99_5]
 stopifnot(all(is.finite(value)),min(s[[col]])==0,abs(as.numeric(quantile(s[[col]],.995))-cap)<1e-12)
 dc<-copy(d);dc[,display_value:=pmin(value,cap)];setorder(dc,display_value,cell_id)
 plots[[ax]]<-ggplot(dc,aes(umap_1,umap_2,colour=display_value))+geom_point(size=.30,stroke=0,alpha=1)+scale_colour_gradientn(colours=pal[[ax]],values=c(0,.18,.4,.7,1),limits=c(0,cap),oob=squish,breaks=c(0,cap/2,cap),labels=c("0",sprintf("%.3f",cap/2),paste0("\u2265",sprintf("%.3f",cap))),name="Shifted module score")+coord_equal(xlim=xr,ylim=yr)+th+labs(title=paste("Sarcosine",tolower(ax)),subtitle=if(ax=="Production")"GNMT + DMGDH" else "SARDH + PIPOX")+guides(colour=guide_colourbar(title.position="top",barwidth=grid::unit(57,"mm"),barheight=grid::unit(4,"mm"),frame.colour="#CED6DC",ticks.colour="#586D77"))
 audits[[ax]]<-data.table(axis=ax,scale_min=0,cap=cap,cap_reference_cells=nrow(s),plotted_cells=nrow(d),above_cap=sum(value>cap),palette=paste(pal[[ax]],collapse="|"),positions="0|0.18|0.4|0.7|1",point_size=.30)
}
E <- (plots$Production|plots$Degradation)+plot_annotation(caption="90,512 annotated cells | 15 patients | Original UMAP and scores retained\nWhite = 0; colour saturation at the original 99.5th-percentile cap (separate scales).",theme=theme(plot.caption=element_text(size=9,family="Arial",colour="#586D77",hjust=0)))
all<-list(Fig6d_Prod_TopBottom3=D,Fig6e_White_Modules=E,Production_White=plots$Production,Degradation_White=plots$Degradation,Lineages_All15=L)
dims<-list(c(12,6.7),c(11,5.5),c(5.5,5.1),c(5.5,5.1),c(6,6.6))
for(i in seq_along(all)){
 nm<-names(all)[i];wh<-dims[[i]]
 ggsave(file.path(out,"figures",paste0(nm,".png")),all[[i]],width=wh[1],height=wh[2],units="in",dpi=400,device=ragg::agg_png,bg="white")
 ggsave(file.path(out,"figures",paste0(nm,".pdf")),all[[i]],width=wh[1],height=wh[2],units="in",device=cairo_pdf,bg="white")
}
fwrite(m,file.path(out,"tables/Patient_Production_Ranking.csv"),bom=TRUE)
fwrite(sel[,.(cell_id,Patient,Resource,selection,final_lineage,umap_1,umap_2)],file.path(out,"tables/Fig6d_Plotted_Cells.csv.gz"))
fwrite(d[,.(cell_id,Patient,final_lineage,umap_1,umap_2,production_module_shifted,degradation_module_shifted)],file.path(out,"tables/Fig6e_Plotted_Cells.csv.gz"))
fwrite(sel[,.(cells=.N),by=.(selection,Patient,Resource,final_lineage)],file.path(out,"tables/Selected_Patient_Cell_Counts.csv"),bom=TRUE)
fwrite(rbindlist(audits),file.path(out,"tables/Colour_Scales.csv"))
fwrite(data.table(lineage=levels,colour=unname(cols)),file.path(out,"tables/Lineage_Colours.csv"))
write_json(list(all_cells=nrow(s),annotated_cells=nrow(d),low_patients=m[selection=="Low",Patient],high_patients=m[selection=="High",Patient],selected_cells=nrow(sel),group_counts=sel[,.(patients=uniqueN(Patient),cells=.N),by=selection],production_score_reconstruction_max_error=max(abs(pcheck-fread(meta_path)$Production_score_ref15)),coordinates_lineages_exact_match=TRUE,new_statistics=FALSE,UMAP_recomputed=FALSE),file.path(out,"qa/Checks.json"),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"qa/sessionInfo.txt"))
print(sel[,.(patients=uniqueN(Patient),cells=.N),by=selection]);cat("OUTPUT",out,"\n")
