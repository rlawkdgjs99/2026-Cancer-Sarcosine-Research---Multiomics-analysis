suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
arg <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
root <- normalizePath(file.path(dirname(arg),'..'))
tab <- file.path(root,'results/tables'); out <- file.path(root,'results/figures')
d <- fread(file.path(tab,'03_HR_selected_plotdata.csv'))
labels <- c(
 REACTOME_NEUTROPHIL_DEGRANULATION='Neutrophil degranulation',
 HALLMARK_INTERFERON_GAMMA_RESPONSE='Interferon-γ response',
 REACTOME_SIGNALING_BY_INTERLEUKINS='Signaling by interleukins',
 HALLMARK_COMPLEMENT='Complement',
 HALLMARK_INFLAMMATORY_RESPONSE='Inflammatory response',
 REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION='Asparagine N-linked\nglycosylation',
 HALLMARK_INTERFERON_ALPHA_RESPONSE='Interferon-α response',
 REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION='Antigen processing:\ncross-presentation',
 HALLMARK_OXIDATIVE_PHOSPHORYLATION='Oxidative phosphorylation',
 REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION='Platelet activation, signaling\nand aggregation',
 REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX='Degradation of the\nextracellular matrix',
 HALLMARK_TNFA_SIGNALING_VIA_NFKB='TNF-α signaling via NF-κB',
 REACTOME_MITOCHONDRIAL_TRANSLATION='Mitochondrial translation',
 REACTOME_CLASS_A_1_RHODOPSIN_LIKE_RECEPTORS='Class A/1\n(rhodopsin-like receptors)',
 REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE='Plasma lipoprotein assembly,\nremodeling and clearance',
 REACTOME_COMPLEX_I_BIOGENESIS='Complex I biogenesis'
)
missing<-setdiff(d$pathway,names(labels)); labels[missing]<-vapply(missing,function(id)paste(strwrap(tools::toTitleCase(tolower(gsub('_',' ',sub('^(HALLMARK|REACTOME)_','',id)))),width=34),collapse='\n'),character(1)); d[,compartment:='Epithelial']
d[,display_label:=paste0(unname(labels[pathway]),ifelse(collection=='Hallmark',' [H]',' [R]'))]
fwrite(d,file.path(tab,'03_HR_selected_plotdata.csv'),quote=TRUE,bom=TRUE)
names_ct<-c(Epithelial='Epithelial cells')
stems<-c(Epithelial='E')
qmax<-c(Epithelial=max(2,ceiling(max(d$neglog10_BH_q)/2)*2))
qbreaks<-list(Epithelial=pretty(c(0,qmax[1]),n=4)); qbreaks[[1]]<-qbreaks[[1]][qbreaks[[1]]<=qmax[1]]
ink <- '#233743'; muted <- '#627581'; red <- '#B32E46'
makeplot <- function(ct) {
 dd <- d[compartment==ct][order(display_order)]
 stopifnot(nrow(dd)>0,nrow(dd)<=8,max(dd$neglog10_BH_q)<=qmax[ct])
 dd[,y:=9-display_order]
 caption <- sprintf('High n = %d, Low n = %d | %s\nFull-collection BH: %s tests\n[H] Hallmark; [R] Reactome | Up to 8 terms after overlap filtering',dd$High[1],dd$Low[1],'Post-treatment NSCLC; ≥50 cells/patient',format(dd$bh_family_n[1],big.mark=','))
 ggplot(dd,aes(x=NES,y=y))+
  geom_point(aes(fill=neglog10_BH_q,size=leading_edge_n),shape=21,colour='white',stroke=.35)+
  scale_y_continuous(breaks=dd$y,labels=dd$display_label,limits=c(.45,8.55),expand=c(0,0))+
  scale_x_continuous(breaks=0:max(4,ceiling(max(d$NES))),limits=c(0,max(4,ceiling(max(d$NES)))),expand=c(0,0))+
  scale_fill_gradientn(colours=c('#F8DDD5','#E89B9A','#C65367','#A11836'),limits=c(0,qmax[ct]),breaks=qbreaks[[ct]],name=expression(-log[10]("BH q")),guide=guide_colourbar(display='rectangles',order=1,barwidth=unit(2.8,'mm'),barheight=unit(27,'mm'),title.position='top',ticks=TRUE,frame.colour='#65737B',frame.linewidth=.25))+
  scale_size_area(max_size=7.3,limits=c(0,max(300,ceiling(max(d$leading_edge_n)/50)*50)),breaks=c(50,150,250),name='Leading-edge\ngenes',guide=guide_legend(order=2,title.position='top',override.aes=list(fill='#344D58',colour='white')))+
  labs(title=names_ct[ct],subtitle='Degradation-High | Hallmark + Reactome',x='Normalized enrichment score (NES)',y=NULL,caption=caption)+
  theme_classic(base_size=11,base_family='Arial')+
  theme(plot.title.position='plot',plot.caption.position='plot',
    plot.title=element_text(face='bold',colour=ink,size=15,margin=margin(b=5)),
    plot.subtitle=element_text(face='bold',colour=red,size=10.5,margin=margin(b=15)),
    axis.title.x=element_text(size=10.5,colour=ink,margin=margin(t=8)),
    axis.text.x=element_text(size=10,colour=ink),axis.text.y=element_text(size=10.3,colour=ink,lineheight=1.02,margin=margin(r=8)),
    axis.ticks.y=element_blank(),axis.line.y=element_blank(),axis.line.x=element_line(colour='#72818A',linewidth=.35),
    axis.ticks.x=element_line(colour='#72818A',linewidth=.3),
    legend.position='right',legend.box='vertical',legend.box.just='left',legend.justification=c(0,1),
    legend.title=element_text(size=9,colour=ink),legend.text=element_text(size=9,colour=ink),
    legend.key.height=unit(5,'mm'),legend.key.width=unit(6,'mm'),legend.spacing.y=unit(5,'mm'),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,8),
    plot.caption=element_text(size=8.1,colour=muted,hjust=0,lineheight=1.18,margin=margin(t=14)),plot.margin=margin(14,12,10,12))
}
saveplot <- function(p,stem,w,h){
 ggsave(file.path(out,paste0(stem,'.png')),p,device=ragg::agg_png,width=w,height=h,units='in',dpi=600,bg='white')
 grDevices::cairo_pdf(file=file.path(out,paste0(stem,'.pdf')),width=w,height=h,bg='white',family='Arial')
 print(p);invisible(dev.off())
}
plots <- lapply(names(names_ct),makeplot)
for(j in seq_along(plots))saveplot(plots[[j]],stems[j],6.8,5.9)
capture.output(sessionInfo(),file=file.path(root,'logs/sessionInfo_plots.txt'))
