#!/usr/bin/env Rscript
# Author correction: combine the same selected six patients on ONE lineage UMAP.
a<-commandArgs(FALSE);f<-sub("^--file=","",a[grepl("^--file=",a)]);b<-normalizePath(file.path(dirname(f),".."))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg);library(jsonlite)})
out<-file.path(b,"results/figures_publication/Fig6_UMAP_Redraw_26.09.29")
d<-fread(file.path(out,"tables/Fig6d_Plotted_Cells.csv.gz"))
s<-as.data.table(readRDS(file.path(b,"intermediate/01_cell_paper_style_scores.rds")))
i<-match(d$cell_id,s$cell_id)
stopifnot(nrow(d)==41685,uniqueN(d$Patient)==6,!anyDuplicated(d$cell_id),!anyNA(i),all(s$analysis_eligible[i]),all(d$umap_1==s$umap_1[i]),all(d$umap_2==s$umap_2[i]),all(d$final_lineage==s$final_lineage[i]),setequal(unique(d$Patient),c("P12","P07","P04","P08","P09","P15")))
lv<-c("Epithelial","CAF","B cell","Plasma cell","CD4 T cell","CD8 T cell","Cycling T cell","NK cell","Mast cell","Neutrophil","Monocyte","Macrophage","Conventional DC","pDC")
co<-setNames(c("#F02C55","#D65720","#C28B00","#878600","#78B300","#009940","#00AA85","#009DC9","#0071D9","#3E36CA","#8153D4","#B14DD5","#AD158F","#EE54B1"),lv)
d[,final_lineage:=factor(final_lineage,levels=lv)];stopifnot(!anyNA(d$final_lineage))
set.seed(260929);d<-d[sample(.N)]
p<-ggplot(d,aes(umap_1,umap_2,colour=final_lineage))+geom_point(size=.48,stroke=0,alpha=1)+scale_colour_manual(values=co,drop=FALSE)+coord_equal(xlim=range(s[analysis_eligible==TRUE]$umap_1),ylim=range(s[analysis_eligible==TRUE]$umap_2))+theme_void(base_family="Arial",base_size=12)+labs(title="Cell lineages",colour=NULL)+theme(plot.title=element_text(size=18,face="bold",colour="black",margin=margin(b=8)),legend.position="right",legend.text=element_text(size=11.5,colour="black"),legend.key.height=grid::unit(5.3,"mm"),legend.key.width=grid::unit(4,"mm"),legend.margin=margin(0,0,0,4),plot.margin=margin(10,10,10,10),plot.background=element_rect(fill="white",colour=NA))+guides(colour=guide_legend(ncol=1,override.aes=list(size=3,alpha=1)))
ggsave(file.path(out,"figures/Fig6d_Single_Vivid.png"),p,width=8,height=5.8,dpi=500,device=ragg::agg_png,bg="white")
ggsave(file.path(out,"figures/Fig6d_Single_Vivid.pdf"),p,width=8,height=5.8,device=cairo_pdf,bg="white")
fwrite(data.table(lineage=lv,colour=unname(co)),file.path(out,"tables/Single_Vivid_Lineage_Colours.csv"))
write_json(list(cells=nrow(d),patients=sort(unique(d$Patient)),one_combined_UMAP=TRUE,coordinates_and_lineages_verified=TRUE,point_size=.48,alpha=1,seed=260929,selection_unchanged=TRUE),file.path(out,"qa/Single_Vivid_Checks.json"),pretty=TRUE,auto_unbox=TRUE)
cat(out,"\n")
