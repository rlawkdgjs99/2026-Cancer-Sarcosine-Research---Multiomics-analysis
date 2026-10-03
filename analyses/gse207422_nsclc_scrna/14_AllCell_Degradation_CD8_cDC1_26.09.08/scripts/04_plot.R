#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])),"common.R"))
dir.create(file.path(out,"resources/font_cache"),showWarnings=FALSE)
Sys.setenv(XDG_CACHE_HOME=file.path(out,"resources/font_cache"))
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
stopifnot(readLines(file.path(out,"results/qa/verification_summary.txt"))[1]=="PASS")
theme_set(theme_classic(base_size=12,base_family="Arial")+theme(plot.title=element_text(face="bold",size=15),
 plot.subtitle=element_text(size=11,colour="#505660"),axis.text=element_text(colour="#26333B"),
 plot.caption=element_text(size=9,colour="#505660",hjust=0),legend.position="bottom",plot.margin=margin(12,16,12,12)))
hi<-"#C94B60";lo<-"#268F95";ink<-"#233B46"
fmt<-function(x)ifelse(is.na(x),"NE",ifelse(x<.001,formatC(x,format="e",digits=1),formatC(x,format="f",digits=3)))
nlab<-c(CD8="CD8 T cells",cDC1="cDC1-like cells")
focus[,display:=c("APC cross-presentation","TCR signaling","TNFR2-related\nnoncanonical NF-kB","IFNG response",
 "Antigen processing / presentation","Regulation of T-cell activation","Type-II IFN production","TNFA signaling via NF-kB","Inflammatory response")]
savefig<-function(p,name,w,h){
 ggsave(file.path(out,"results/figures",paste0(name,".png")),p,width=w,height=h,dpi=400,device=ragg::agg_png,bg="white")
 ggsave(file.path(out,"results/figures",paste0(name,".pdf")),p,width=w,height=h,device=cairo_pdf,bg="white")
 ggsave(file.path(out,"results/figures",paste0(name,".svg")),p,width=w,height=h,device=grDevices::svg,bg="white")
}
md<-fread(file.path(out,"results/tables/01_FROZEN_allcell_patient_scores_groups.csv"))
md[,group:=factor(group,levels=c("Low","High"))];md[,Patient_order:=factor(Patient,levels=Patient[order(degradation_mean_z)])]
md[,treatment_label:=factor(treatment,levels=c("Pre","Post"),labels=c("Pre-treatment biopsy","Post-treatment surgery"))]
cols<-scale_colour_manual(values=c(Low=lo,High=hi),labels=c(Low="Degradation-Low",High="Degradation-High"),name="Whole-sample group")
fills<-scale_fill_manual(values=c(Low=lo,High=hi),labels=c(Low="Degradation-Low",High="Degradation-High"),name="Whole-sample group")
p1a<-ggplot(md,aes(Patient_order,degradation_mean_z,colour=group))+geom_segment(aes(xend=Patient_order,y=0,yend=degradation_mean_z),linewidth=1.4,alpha=.65)+
 geom_hline(yintercept=unique(md$cutoff),linetype=2,colour="#717D84",linewidth=.5)+
 geom_point(aes(shape=treatment_label),size=3.5,stroke=1.1)+cols+
 scale_shape_manual(values=c("Pre-treatment biopsy"=2,"Post-treatment surgery"=16),name="Specimen")+
 guides(colour=guide_legend(order=1,nrow=1,override.aes=list(shape=16)),shape=guide_legend(order=2,nrow=1))+
 labs(title="One patient-level exposure",subtitle="All retained cell types pooled within each patient",x="Patient",y="SARDH/PIPOX mean-z score")+
 theme(axis.text.x=element_text(angle=45,hjust=1))
p1b<-ggplot(md,aes(SARDH_log2CPM,PIPOX_log2CPM,colour=group,shape=treatment_label))+geom_point(size=3.4,stroke=1.1)+
 ggrepel::geom_text_repel(aes(label=Patient),size=3,show.legend=FALSE,seed=260908,max.overlaps=Inf,min.segment.length=0,box.padding=.4)+cols+
 scale_shape_manual(values=c("Pre-treatment biopsy"=2,"Post-treatment surgery"=16),name="Specimen")+
 scale_y_continuous(expand=expansion(mult=c(.08,.16)))+
 labs(title="Both score components are observed",subtitle="Same two-gene formula; no CD8- or cDC1-specific regrouping",x="SARDH  log2(1 + TMM CPM)",y="PIPOX  log2(1 + TMM CPM)")+guides(colour="none",shape="none")
p1<-(p1a+p1b)+plot_layout(guides="collect")+plot_annotation(title="A shared tumour-sample degradation grouping for cell-specific analyses",
 subtitle="15 independent patients | 91,844 retained singlets | High 7 / Low 8",
 caption="Median cutoff fixed across all 15 patients. Primary Post-only comparisons retain these labels: CD8 High 7 / Low 5; cDC1 High 6 / Low 4.\nThis is an all-cell RNA proxy, not measured sarcosine concentration. Dashed line: the fixed cohort median.")
p1<-p1 & theme(legend.position="bottom",legend.box="vertical",legend.text=element_text(size=10),legend.title=element_text(size=10))
savefig(p1,"Fig_01_AllCell_Exposure_and_Fixed_Patient_Groups",12.6,6.1);wt(md,"plotdata_01_exposure.csv")
a<-fread(file.path(out,"results/tables/05_focus_complete.csv"));walks<-fread(file.path(out,"results/tables/06_primary_running_ES.csv"))
elig<-fread(file.path(out,"results/tables/01_model_eligibility.csv"));allg<-fread(file.path(out,"results/tables/02_GSEA_all_scopes.csv"))
nf<-allg[scope=="post_group",.N]
for(ct in c("CD8","cDC1")){
 d<-a[cell_type==ct&scope=="post_group"&role=="Fig6 exact"];nn<-elig[cell_type==ct&scope=="post_group"]
 plots<-lapply(focus$pathway[1:4],function(id){
  rr<-d[pathway==id];ttl<-focus$display[match(id,focus$pathway)]
  if(rr$status!="TESTED")return(ggplot()+annotate("text",x=.5,y=.5,label=paste0("Not estimable\n",rr$available_genes," genes; minimum 15"),size=5)+labs(title=ttl)+theme_void())
  w<-walks[cell_type==ct&pathway==id];rg<-range(w$running_ES);tick<-rg[1]-.07
  ggplot(w,aes(rank,running_ES))+geom_hline(yintercept=0,colour="#AEB8BE",linetype=2,linewidth=.45)+
   geom_line(colour=if(rr$NES>0)hi else lo,linewidth=.9)+
   geom_segment(data=w[hit==TRUE],aes(x=rank,xend=rank,y=tick,yend=tick+.025),inherit.aes=FALSE,linewidth=.23,alpha=.7,colour=ink)+
   scale_x_continuous(expand=expansion(mult=c(.01,.015)))+
   labs(title=ttl,subtitle=paste0("NES ",sprintf("%.2f",rr$NES),"  |  BH q ",fmt(rr$q_global)),
    x="Gene rank: Degradation-High to Degradation-Low",y="Running enrichment score")+
   theme(plot.title=element_text(size=13),plot.subtitle=element_text(size=11),axis.title=element_text(size=9.5),axis.text=element_text(size=9))
 })
 p<-wrap_plots(plots,ncol=2)+plot_annotation(title=paste0(nlab[ct],": immune programs by whole-sample degradation state"),
  subtitle=paste0("Post-treatment patients: High ",nn$High," / Low ",nn$Low," | ",format(nn$n_cells,big.mark=",")," cells | Patient-level pseudobulk"),
  caption=paste0("Histology-adjusted High-minus-Low ranks; the same all-cell patient grouping is used throughout.\nBH q covers all ",format(nf,big.mark=",")," tested sets across both cell types and GO:BP, Hallmark, Reactome.\n",
   if(ct=="cDC1")"Exploratory: >=5 cDC1-like cells/patient. "else">=50 CD8 cells/patient. ",
   "CAMERA and patient-score comparisons do not pass BH q < 0.05."),
  theme=theme(plot.title=element_text(face="bold",size=16),plot.subtitle=element_text(size=11),plot.caption=element_text(size=9,hjust=0)))
 savefig(p,paste0(if(ct=="CD8")"Fig_02"else"Fig_03","_",ct,"_AllCell_Groups_Exact_GSEA"),11.8,7.7)
 wt(d,paste0("plotdata_",ct,"_primary_exact.csv"))
}
d<-copy(a[scope=="post_group"]);d[,display:=factor(focus$display[match(pathway,focus$pathway)],levels=rev(focus$display))]
d[,panel:=factor(cell_type,levels=c("CD8","cDC1"),labels=c("CD8 T cells\nHigh 7 / Low 5","cDC1-like cells\nHigh 6 / Low 4"))]
d[,txt:=ifelse(status=="TESTED",paste0(sprintf("%.2f",NES),"\nq ",fmt(q_global)),"NE")]
lim<-max(3,max(abs(d$NES),na.rm=TRUE));lim<-ceiling(lim*2)/2
p4<-ggplot(d,aes(panel,display))+geom_tile(aes(fill=NES),colour="white",linewidth=.8)+
 geom_tile(data=d[!is.na(q_global)&q_global<.05],fill=NA,colour=ink,linewidth=1)+
 geom_text(aes(label=txt),size=3.4,lineheight=1.1)+geom_hline(yintercept=5.5,colour=ink,linewidth=.7)+
 scale_fill_gradient2(low=lo,mid="white",high=hi,limits=c(-lim,lim),na.value="#EDF0F2",name="NES")+
 scale_x_discrete(position="top")+
 labs(title="Cell-specific immune programs under one patient grouping",
  subtitle="Whole-sample Degradation-High versus Low | Post-treatment primary analysis",
  x=NULL,y=NULL,caption=paste0("Dark border: BH q < 0.05 across all ",format(nf,big.mark=",")," tested sets. Above divider: four exact Figure 6 programs.\nBelow divider: five prespecified supporting themes. Positive NES denotes Degradation-High.\ncDC1 is a low-cell exploratory analysis; enrichment is not a demonstrated causal sequence."))+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(face="bold",size=11),axis.text.y=element_text(size=11))
savefig(p4,"Fig_04_CD8_cDC1_Shared_Group_Immune_Themes",10.7,7.8);wt(d,"plotdata_04_primary_themes.csv")
sv<-fread(file.path(out,"results/tables/04_patient_score_values.csv"))[scope=="post_group"]
sr<-fread(file.path(out,"results/tables/04_exact_scores_HC3.csv"))[scope=="post_group"]
ylim<-range(sv$score,na.rm=TRUE);ylim<-c(floor(ylim[1]*2)/2-.25,ceiling(ylim[2]*2)/2+.3)
p5rows<-lapply(c("CD8","cDC1"),function(ct){
 panels<-lapply(focus$pathway[1:4],function(id){
  v<-copy(sv[cell_type==ct&pathway==id]);v[,group:=factor(group,levels=c("Low","High"))];rs<-sr[cell_type==ct&pathway==id];nn<-elig[cell_type==ct&scope=="post_group"]
  ggplot(v,aes(group,score,colour=group))+geom_hline(yintercept=0,colour="#CDD4D9",linetype=2)+
   geom_boxplot(width=.48,outlier.shape=NA,linewidth=.65,fill="white")+
   geom_point(position=position_jitter(width=.065,height=0,seed=260908),size=2.6,alpha=.85)+
   scale_colour_manual(values=c(Low=lo,High=hi),guide="none")+
   scale_x_discrete(labels=c(Low=paste0("Low\nn=",nn$Low),High=paste0("High\nn=",nn$High)))+
   coord_cartesian(ylim=ylim)+
   labs(title=focus$display[match(id,focus$pathway)],subtitle=paste0("HC3 BH q ",fmt(rs$q_eight_scores)),x=NULL,y=if(id==focus$pathway[1])"Standardized pathway score" else NULL)+
   theme(plot.title=element_text(size=11),plot.subtitle=element_text(size=10),axis.title=element_text(size=10),axis.text=element_text(size=9))
 })
 wrap_plots(panels,ncol=4)+plot_annotation(title=nlab[ct],theme=theme(plot.title=element_text(face="bold",size=14)))
})
p5<-(wrap_elements(full=p5rows[[1]])/wrap_elements(full=p5rows[[2]]))+plot_annotation(title="Patient-level pathway scores by whole-sample degradation group",
 subtitle="Post-treatment primary comparisons | One point per patient",
 caption="Low / High always refer to the all-cell SARDH/PIPOX grouping. Outcome scores are standardized within each cell type.\nDisplayed q values: histology-adjusted HC3 comparisons, BH across eight score tests; these are distinct from ranked-list GSEA.",
 theme=theme(plot.title=element_text(face="bold",size=16),plot.subtitle=element_text(size=11),plot.caption=element_text(hjust=0,size=9)))
savefig(p5,"Fig_05_CD8_cDC1_Primary_Patient_Pathway_Scores",13.8,8.7);wt(merge(sv,sr[,.(cell_type,pathway,beta,p,q_eight_scores)],by=c("cell_type","pathway")),"plotdata_05_primary_scores.csv")
s<-copy(a[role=="Fig6 exact"]);s[,display:=factor(focus$display[match(pathway,focus$pathway)],levels=rev(focus$display[1:4]))]
s[,contrast:=factor(scope,levels=c("post_group","post_continuous","all_group","all_continuous"),labels=c("Post High-Low\nPRIMARY","Post continuous\nSensitivity","All High-Low\nSensitivity","All continuous\nSensitivity"))]
s[,cell_label:=factor(cell_type,levels=c("CD8","cDC1"),labels=unname(nlab))]
s[,txt:=ifelse(status=="TESTED",paste0(sprintf("%.2f",NES),"\nq ",fmt(q_global)),"NE")]
lim<-ceiling(max(3,abs(s$NES),na.rm=TRUE)*2)/2
p6<-ggplot(s,aes(contrast,display))+geom_tile(aes(fill=NES),colour="white",linewidth=.75)+
 geom_tile(data=s[!is.na(q_global)&q_global<.05],fill=NA,colour=ink,linewidth=1)+geom_text(aes(label=txt),size=3.25,lineheight=1.15)+
 facet_grid(cell_label~.,switch="y")+scale_x_discrete(position="top")+
 scale_fill_gradient2(low=lo,mid="white",high=hi,limits=c(-lim,lim),name="NES",na.value="#EDF0F2")+
 labs(title="Exact immune programs across fixed analysis checks",
 subtitle="One frozen all-cell exposure; no patient regrouping across cell types or scopes",x=NULL,y=NULL,
 caption="Border: global BH q < 0.05 within the indicated scope. The primary result remains Post High-Low regardless of sensitivity results.\nAll scopes use the same15-patient score reference. Post: histology adjustment; all patients: treatment and histology adjustment.\nThese checks reuse the same dataset and do not constitute independent validation.")+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(face="bold",size=10),strip.background=element_blank(),strip.text=element_text(face="bold",size=12))
savefig(p6,"Fig_S01_Fixed_Contrast_Sensitivity_Checks",12,8);wt(s,"plotdata_S01_sensitivity.csv")
ca<-fread(file.path(out,"results/tables/03_CAMERA_primary.csv"))
cross<-merge(merge(a[scope=="post_group"&role=="Fig6 exact"],ca[,.(cell_type,pathway,CAMERA_p=PValue,CAMERA_q=q_global)],by=c("cell_type","pathway")),sr[,.(cell_type,pathway,score_beta=beta,score_ci_low=ci_low,score_ci_high=ci_high,score_p=p,score_q=q_eight_scores)],by=c("cell_type","pathway"))
wt(cross,"07_primary_cross_method_summary.csv")
writeLines("Six analysis plots exported as PNG/PDF/SVG, including separate CD8 and cDC1 GSEA figures; visual QA pending.",file.path(out,"logs/plot_status.txt"))
cat("FIGURE EXPORT COMPLETE; visual QA pending\n")
