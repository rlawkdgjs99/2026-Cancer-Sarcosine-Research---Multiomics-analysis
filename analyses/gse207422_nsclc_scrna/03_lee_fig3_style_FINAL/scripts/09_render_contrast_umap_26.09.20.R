#!/usr/bin/env Rscript
# Display-only redraw from frozen scores/coordinates; no reanalysis.
a <- commandArgs(FALSE); script <- sub("^--file=", "", a[grepl("^--file=", a)])
b <- normalizePath(file.path(dirname(script), ".."))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(scales);library(ragg);library(jsonlite)})
out <- file.path(b,"results/figures_publication/UMAP_Contrast_26.09.20")
s <- as.data.table(readRDS(file.path(b,"intermediate/01_cell_paper_style_scores.rds")))
lin <- fread(file.path(b,"../02_lineage_reannotation/results/tables/09_final_cell_lineages_FROZEN.csv"))
caps <- fread(file.path(b,"results/tables/plotdata_01_UMAP_display_caps.csv"))
stopifnot(nrow(s)==92053, !anyDuplicated(s$cell_id), all(is.finite(s$umap_1)),all(is.finite(s$umap_2)))
i <- match(s$cell_id,lin$cell_id)
stopifnot(!anyNA(i),identical(s$analysis_eligible,lin$analysis_eligible[i]),identical(s$final_lineage,lin$final_lineage[i]),identical(s$umap_1,lin$umap_1[i]),identical(s$umap_2,lin$umap_2[i]))
d <- copy(s[analysis_eligible==TRUE])
stopifnot(nrow(d)==90512,uniqueN(d$Patient)==15)
axes <- c("Production","Degradation")
for(ax in axes){
 col <- paste0(tolower(ax),"_module_shifted")
 stopifnot(all(is.finite(s[[col]])),min(s[[col]])==0,
 abs(quantile(s[[col]],.995)-caps[axis==ax,display_cap_99_5])<1e-12)
}
keep <- c("cell_id","Sample","Patient","Resource","final_lineage","analysis_eligible","umap_1","umap_2","production_module_raw","degradation_module_raw","production_module_shifted","degradation_module_shifted")
fwrite(s[,..keep],file.path(out,"tables/source_cells.csv.gz"))
fwrite(d[,..keep],file.path(out,"tables/plotted_cells.csv.gz"))
fwrite(d[,.(cells=.N),by=.(Patient,Resource,final_lineage)],file.path(out,"tables/cell_counts.csv"))
level_order <- c("Epithelial","CAF","B cell","Plasma cell","CD4 T cell","CD8 T cell","Cycling T cell","NK cell","Mast cell","Neutrophil","Monocyte","Macrophage","Conventional DC","pDC")
cols <- setNames(hcl.colors(length(level_order),"Dark 3"),level_order)
d[,final_lineage:=factor(final_lineage,levels=level_order)]
stopifnot(!anyNA(d$final_lineage))
xr <- range(d$umap_1); yr <- range(d$umap_2)
base_theme <- theme_void(base_family="Arial",base_size=10)+theme(
 plot.title=element_text(size=13,face="bold",colour="#24343C",margin=margin(b=4)),
 plot.subtitle=element_text(size=9,colour="#53636B",margin=margin(b=8)),
 legend.position="bottom",legend.title=element_text(size=8.5),legend.text=element_text(size=8),
 plot.margin=margin(8,8,8,8),plot.background=element_rect(fill="white",colour=NA))
set.seed(260920); dl <- d[sample(.N)]
L <- ggplot(dl,aes(umap_1,umap_2,colour=final_lineage))+geom_point(size=.18,stroke=0,alpha=1)+
 scale_colour_manual(values=cols,drop=FALSE)+coord_equal(xlim=xr,ylim=yr)+base_theme+
 labs(title="Cell lineages",subtitle="90,512 annotated cells | 15 patients",colour=NULL)+
 guides(colour=guide_legend(ncol=2,byrow=TRUE,override.aes=list(size=2,alpha=1)))+
 theme(legend.key.height=grid::unit(2.8,"mm"),legend.key.width=grid::unit(3.5,"mm"))
# Same global scale per module for all lineages. Colour changes only;
# original cap is from all 92,053 frozen cells, not recalculated after subset.
pal <- list(Production=c("#EFF0F2","#DADDE1","#F49C37","#CC4C02","#7F2704"),Degradation=c("#EFF0F2","#DADDE1","#4292C6","#08519C","#041F49"))
plots <- list(); audits <- list()
for(ax in axes){
 dc <- copy(d); value <- dc[[paste0(tolower(ax),"_module_shifted")]]
 cap <- caps[axis==ax,display_cap_99_5];dc[,value_display:=pmin(value,cap)]
 setorder(dc,value_display,cell_id)
 breaks <- c(0,cap/2,cap)
 labels <- c("0",sprintf("%.3f",cap/2),paste0("≥",sprintf("%.3f",cap)))
 plots[[ax]] <- ggplot(dc,aes(umap_1,umap_2,colour=value_display))+
  geom_point(size=.30,stroke=0,alpha=1)+
  scale_colour_gradientn(colours=pal[[ax]],values=c(0,.18,.4,.7,1),limits=c(0,cap),oob=squish,breaks=breaks,labels=labels,name="Shifted module score")+
  coord_equal(xlim=xr,ylim=yr)+base_theme+
  labs(title=paste("Sarcosine",tolower(ax)),subtitle=if(ax=="Production")"GNMT + DMGDH" else "SARDH + PIPOX")+
  guides(colour=guide_colourbar(title.position="top",barwidth=grid::unit(48,"mm"),barheight=grid::unit(3,"mm")))
 audits[[ax]] <- data.table(axis=ax,minimum=0,display_cap=cap,cap_source_cells=nrow(s),plotted_cells=nrow(d),above_cap=sum(value>cap),at_or_above_cap=sum(value>=cap),palette=paste(pal[[ax]],collapse="|"),colour_positions="0|0.18|0.4|0.7|1")
}
cap_note <- "Original scores and UMAP coordinates retained; colours saturate at the original 99.5th percentile."
DP <- (plots$Production | plots$Degradation)+plot_annotation(caption=cap_note,theme=theme(plot.caption=element_text(size=8,colour="#53636B",hjust=0)))
M <- (L | plots$Production | plots$Degradation)+plot_annotation(caption=cap_note,theme=theme(plot.caption=element_text(size=8,colour="#53636B",hjust=0)))
all <- list(L=L,P=plots$Production,D=plots$Degradation,DP=DP,M=M)
dims <- list(L=c(4.4,5.1),P=c(4.4,4.4),D=c(4.4,4.4),DP=c(8.8,4.6),M=c(13.2,5.1))
for(nm in names(all)){
 wh <- dims[[nm]]
 ggsave(file.path(out,"figures/PPT_insert",paste0(nm,".png")),all[[nm]],width=wh[1],height=wh[2],units="in",dpi=600,device=ragg::agg_png,bg="white")
 ggsave(file.path(out,"figures",paste0(nm,".pdf")),all[[nm]],width=wh[1],height=wh[2],units="in",device=cairo_pdf,bg="white")
}
fwrite(rbindlist(audits),file.path(out,"tables/display_specification.csv"))
write_json(list(source_cells=nrow(s),plotted_cells=nrow(d),excluded_unassigned_cells=nrow(s)-nrow(d),patients=uniqueN(d$Patient),lineages=uniqueN(d$final_lineage),new_tests=FALSE,source_coordinates_and_labels_verified=TRUE,score_recalculation=FALSE,lineage_point_size=.18,score_point_size=.30,alpha=1,draw_order="ascending displayed score, then cell_id; lineage shuffled seed260920",production_detected=sum(s[analysis_eligible==TRUE,GNMT_log1pCP10k>0|DMGDH_log1pCP10k>0]),degradation_detected=sum(s[analysis_eligible==TRUE,SARDH_log1pCP10k>0|PIPOX_log1pCP10k>0])),file.path(out,"qa/render_checks.json"),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"qa/sessionInfo.txt"))
cat("COMPLETE",out,"\n")
