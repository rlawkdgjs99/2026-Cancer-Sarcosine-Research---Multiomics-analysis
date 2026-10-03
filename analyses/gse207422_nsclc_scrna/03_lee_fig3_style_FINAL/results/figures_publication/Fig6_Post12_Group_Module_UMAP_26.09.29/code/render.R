args<-commandArgs(trailingOnly=TRUE); b<-normalizePath(args[1]);out<-normalizePath(args[2])
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(scales);library(ragg);library(jsonlite)})
for(x in c("figures","tables","qa"))dir.create(file.path(out,x),showWarnings=FALSE)
sources<-c(file.path(b,"intermediate/01_cell_paper_style_scores.rds"),file.path(b,"results/tables/plotdata_01_UMAP_display_caps.csv"),file.path(b,"results/figures_publication/Fig6_Post12_Group_UMAP_26.09.29/tables/Patient_Groups.csv"),file.path(b,"results/figures_publication/Fig6_Post12_Group_UMAP_26.09.29/tables/Plotted_Cells.csv.gz"))
hashes<-tools::md5sum(sources)
s<-as.data.table(readRDS(sources[1]));caps<-fread(sources[2]);m<-fread(sources[3]);old<-fread(sources[4])
stopifnot(nrow(m)==12,all(m$treatment=="Post"),sum(m$group=="Low")==5,sum(m$group=="High")==7,!anyDuplicated(m$Patient))
d<-copy(s[analysis_eligible==TRUE & Patient %in% m$Patient]);d[,group:=m$group[match(Patient,m$Patient)]]
stopifnot(nrow(d)==76990,!anyNA(d$group),!anyDuplicated(d$cell_id),setequal(d$cell_id,old$cell_id))
i<-match(d$cell_id,old$cell_id)
for(nm in c("Patient","group","umap_1","umap_2"))stopifnot(identical(d[[nm]],old[[nm]][i]))
xr<-range(s[analysis_eligible==TRUE]$umap_1);yr<-range(s[analysis_eligible==TRUE]$umap_2)
pal<-list(Production=c("#FFFFFF","#FFF3DD","#FFBF55","#E56A0A","#943100"),Degradation=c("#FFFFFF","#E4F4FF","#65B9EB","#1776C6","#073B8C"))
plots<-list();audit<-list()
for(ax in names(pal))for(gr in c("Low","High")){
 dc<-copy(d[group==gr]);col<-paste0(tolower(ax),"_module_shifted");cap<-caps[axis==ax,display_cap_99_5]
 stopifnot(length(cap)==1,all(is.finite(dc[[col]])),min(s[[col]])==0,abs(as.numeric(quantile(s[[col]],.995))-cap)<1e-12)
 dc[,display_value:=pmin(get(col),cap)];setorder(dc,display_value,cell_id)
 nm<-paste(ax,gr,sep="_")
 p<-ggplot(dc,aes(umap_1,umap_2,colour=display_value))+geom_point(size=.30,stroke=0,alpha=1)+
 scale_colour_gradientn(colours=pal[[ax]],values=c(0,.18,.4,.7,1),limits=c(0,cap),oob=squish,breaks=c(0,cap/2,cap),labels=c("0",sprintf("%.3f",cap/2),paste0("\u2265",sprintf("%.3f",cap))),name="Shifted module score")+
 coord_equal(xlim=xr,ylim=yr)+theme_void(base_family="Arial",base_size=12)+
 labs(title=paste("Sarcosine",tolower(ax)),subtitle=paste0(if(ax=="Production")"GNMT + DMGDH" else "SARDH + PIPOX","\nDegradation ",gr," (n = ",uniqueN(dc$Patient),")"))+
 theme(plot.title=element_text(size=19,face="bold",colour="#193340",margin=margin(b=6)),plot.subtitle=element_text(size=13,colour="#586D77",lineheight=1.15,margin=margin(b=10)),plot.margin=margin(12,12,12,12),plot.background=element_rect(fill="white",colour=NA),legend.title=element_text(size=12),legend.text=element_text(size=11),legend.position="bottom")+
 guides(colour=guide_colourbar(title.position="top",barwidth=grid::unit(64,"mm"),barheight=grid::unit(4,"mm"),frame.colour="#CED6DC",ticks.colour="#586D77"))
 plots[[nm]]<-p
 ggsave(file.path(out,"figures",paste0(nm,".png")),p,width=6,height=6,dpi=400,device=ragg::agg_png,bg="white")
 ggsave(file.path(out,"figures",paste0(nm,".pdf")),p,width=6,height=6,device=cairo_pdf,bg="white")
 audit[[nm]]<-data.table(metric=ax,group=gr,patients=uniqueN(dc$Patient),cells=nrow(dc),scale_min=0,cap=cap,above_cap=sum(dc[[col]]>cap),point_size=.30,alpha=1,palette=paste(pal[[ax]],collapse="|"))
}
combined<-wrap_plots(plots,ncol=2)
ggsave(file.path(out,"figures/Modules_Low_High.png"),combined,width=12,height=12,dpi=400,device=ragg::agg_png,bg="white")
ggsave(file.path(out,"figures/Modules_Low_High.pdf"),combined,width=12,height=12,device=cairo_pdf,bg="white")
fwrite(d[,.(cell_id,Patient,group,final_lineage,umap_1,umap_2,production_module_shifted,degradation_module_shifted)],file.path(out,"tables/Plotted_Cells.csv.gz"))
fwrite(m,file.path(out,"tables/Patient_Groups.csv"),bom=TRUE)
fwrite(rbindlist(audit),file.path(out,"tables/Plot_Audit.csv"))
stopifnot(identical(hashes,tools::md5sum(sources)))
fwrite(data.table(source=sources,md5=unname(hashes)),file.path(out,"qa/Source_Hashes.csv"))
writeLines(capture.output(sessionInfo()),file.path(out,"qa/sessionInfo.txt"))
print(rbindlist(audit)[,.(metric,group,patients,cells,cap)]);cat("PASS: original cell roster, coordinates, groups, scores and source immutability checks\n")
