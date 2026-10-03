#!/usr/bin/env Rscript
# Author correction: combine the same selected six patients on ONE lineage UMAP.
a<-commandArgs(FALSE);f<-sub("^--file=","",a[grepl("^--file=",a)]);b<-normalizePath(file.path(dirname(f),".."))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ragg);library(jsonlite)})
out<-file.path(b,"results/figures_publication/Fig6_Post12_Group_UMAP_26.09.29")
for(x in c("figures","tables","qa"))dir.create(file.path(out,x),recursive=TRUE,showWarnings=FALSE)
s<-as.data.table(readRDS(file.path(b,"intermediate/01_cell_paper_style_scores.rds")))
meta_path<-file.path(b,"../14_AllCell_Degradation_CD8_cDC1_26.09.08/results/Fig6_ImmuneComposition_26.09.21/tables/frozen_patient_metadata.csv")
m<-fread(meta_path)
stopifnot(nrow(m)==12,!anyDuplicated(m$Patient),all(m$treatment=="Post"),sum(m$group=="Low")==5,sum(m$group=="High")==7)
stopifnot(all(ifelse(m$degradation_mean_z>m$cutoff,"High","Low")==m$group))
d<-copy(s[analysis_eligible==TRUE & Patient %in% m$Patient])
d[,group:=m$group[match(Patient,m$Patient)]]
stopifnot(!anyNA(d$group),!anyDuplicated(d$cell_id),uniqueN(d$Patient)==12)
lv<-c("Epithelial","CAF","B cell","Plasma cell","CD4 T cell","CD8 T cell","Cycling T cell","NK cell","Mast cell","Neutrophil","Monocyte","Macrophage","Conventional DC","pDC")
co<-setNames(c("#F02C55","#D65720","#C28B00","#878600","#78B300","#009940","#00AA85","#009DC9","#0071D9","#3E36CA","#8153D4","#B14DD5","#AD158F","#EE54B1"),lv)
d[,final_lineage:=factor(final_lineage,levels=lv)];stopifnot(!anyNA(d$final_lineage))
make_plot<-function(d,title){
set.seed(260929);d<-d[sample(.N)]
p<-ggplot(d,aes(umap_1,umap_2,colour=final_lineage))+geom_point(size=.48,stroke=0,alpha=1)+scale_colour_manual(values=co,drop=FALSE)+coord_equal(xlim=range(s[analysis_eligible==TRUE]$umap_1),ylim=range(s[analysis_eligible==TRUE]$umap_2))+theme_void(base_family="Arial",base_size=12)+labs(title=title,colour=NULL)+theme(plot.title=element_text(size=18,face="bold",colour="black",margin=margin(b=8)),legend.position="right",legend.text=element_text(size=11.5,colour="black"),legend.key.height=grid::unit(5.3,"mm"),legend.key.width=grid::unit(4,"mm"),legend.margin=margin(0,0,0,4),plot.margin=margin(10,10,10,10),plot.background=element_rect(fill="white",colour=NA))+guides(colour=guide_legend(ncol=1,override.aes=list(size=3,alpha=1)))

return(p)
}

for(gr in c("Low","High")){
 dd<-d[group==gr];n<-uniqueN(dd$Patient)
 p<-make_plot(dd,paste0("Degradation ",gr," (n = ",n,")"))
 ggsave(file.path(out,"figures",paste0(gr,".png")),p,width=8,height=5.8,dpi=500,device=ragg::agg_png,bg="white")
 ggsave(file.path(out,"figures",paste0(gr,".pdf")),p,width=8,height=5.8,device=cairo_pdf,bg="white")
}
fwrite(m,file.path(out,"tables/Patient_Groups.csv"),bom=TRUE)
fwrite(d[,.(cell_id,Sample,Patient,Resource,group,final_lineage,umap_1,umap_2)],file.path(out,"tables/Plotted_Cells.csv.gz"))
fwrite(d[,.(cells=.N),by=.(group,Patient,final_lineage)],file.path(out,"tables/Cell_Counts.csv"),bom=TRUE)
# Panel e clean display: same all15 cells, original modules and original colour caps.
suppressPackageStartupMessages({library(patchwork);library(scales)})
caps<-fread(file.path(b,"results/tables/plotdata_01_UMAP_display_caps.csv"))
z<-copy(s[analysis_eligible==TRUE]);stopifnot(nrow(z)==90512)
pal<-list(Production=c("#FFFFFF","#FFF3DD","#FFBF55","#E56A0A","#943100"),Degradation=c("#FFFFFF","#E4F4FF","#65B9EB","#1776C6","#073B8C"))
plots<-list()
for(ax in names(pal)){
 dc<-copy(z);col<-paste0(tolower(ax),"_module_shifted");cap<-caps[axis==ax,display_cap_99_5]
 dc[,display_value:=pmin(get(col),cap)];setorder(dc,display_value,cell_id)
 plots[[ax]]<-ggplot(dc,aes(umap_1,umap_2,colour=display_value))+geom_point(size=.30,stroke=0,alpha=1)+scale_colour_gradientn(colours=pal[[ax]],values=c(0,.18,.4,.7,1),limits=c(0,cap),oob=squish,breaks=c(0,cap/2,cap),labels=c("0",sprintf("%.3f",cap/2),paste0("\u2265",sprintf("%.3f",cap))),name="Shifted module score")+coord_equal(xlim=range(z$umap_1),ylim=range(z$umap_2))+theme_void(base_family="Arial",base_size=12)+labs(title=paste("Sarcosine",tolower(ax)),subtitle=if(ax=="Production")"GNMT + DMGDH" else "SARDH + PIPOX")+theme(plot.title=element_text(size=17,face="bold",colour="#193340",margin=margin(b=5)),plot.subtitle=element_text(size=10.5,colour="#586D77",margin=margin(b=8)),plot.margin=margin(10,10,10,10),plot.background=element_rect(fill="white",colour=NA),legend.title=element_text(size=11),legend.text=element_text(size=10),legend.position="bottom")+guides(colour=guide_colourbar(title.position="top",barwidth=grid::unit(57,"mm"),barheight=grid::unit(4,"mm"),frame.colour="#CED6DC",ticks.colour="#586D77"))
}
E<-plots$Production|plots$Degradation
for(ext in c("png","pdf")){
 path<-file.path(out,"figures",paste0("E_clean.",ext))
 if(ext=="png")ggsave(path,E,width=11,height=5.1,dpi=400,device=ragg::agg_png,bg="white") else ggsave(path,E,width=11,height=5.1,device=cairo_pdf,bg="white")
}
write_json(list(group_counts=d[,.(patients=uniqueN(Patient),cells=.N),by=group],metadata_source=meta_path,groups_redefined=FALSE,umap_recomputed=FALSE,downsampled=FALSE,point_size=.48,alpha=1,shared_axis_limits=TRUE,footnotes=FALSE,e_cells=nrow(z)),file.path(out,"qa/Checks.json"),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"qa/sessionInfo.txt"))
print(d[,.(patients=uniqueN(Patient),cells=.N),by=group]);cat(out,"\n")
