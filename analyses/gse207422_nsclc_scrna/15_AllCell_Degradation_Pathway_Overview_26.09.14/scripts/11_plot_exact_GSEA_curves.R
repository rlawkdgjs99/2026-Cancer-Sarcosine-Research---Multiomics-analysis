suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
a <- commandArgs(trailingOnly=TRUE)
root <- if(length(a)) normalizePath(a[1],mustWork=TRUE) else dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]))))
tab <- file.path(root,'results','tables'); out <- file.path(root,'results','figures')
v <- fread(file.path(tab,'Fig6_exact_GSEA_curve_vertices.csv'))
h <- fread(file.path(tab,'Fig6_exact_GSEA_gene_hits.csv'))
a <- fread(file.path(tab,'Fig6_exact_GSEA_curve_annotations.csv'))
ids <- c('REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION','REACTOME_TCR_SIGNALING','REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY','HALLMARK_INTERFERON_GAMMA_RESPONSE')
labels <- c('Antigen cross-presentation\n(MHC-I)','TCR signaling','TNFR2-related\nnoncanonical NF-κB','IFN-γ response')
comp <- c('Whole_tumour','CD8','cDC')
compnames <- c('Whole tumour','CD8⁺ T cells','Conventional DCs')
short <- c('W','T','D')
ink <- '#253C49'; muted <- '#667A86'; high <- '#BB4A57'; low <- '#3579AB'
super <- setNames(strsplit('⁰¹²³⁴⁵⁶⁷⁸⁹⁻','')[[1]],strsplit('0123456789-','')[[1]])
qstr <- function(q) {
 if(q>=.001) return(sprintf('%.3f',q))
 s <- strsplit(sprintf('%.1e',q),'e',fixed=TRUE)[[1]]
 paste0(s[1],' × 10',paste0(super[strsplit(as.character(as.integer(s[2])),'')[[1]]],collapse=''))
}
make <- function(cc, ii, compact=FALSE) {
 id <- ids[ii]; vv <- v[compartment==cc & pathway==id]; hh <- h[compartment==cc & pathway==id]; aa <- a[compartment==cc & pathway==id]
 stopifnot(nrow(aa)==1)
 colour <- if(aa$NES>0) high else low
 fs <- if(compact) 10 else 11
 ggplot(vv,aes(rank_percent,running_ES))+
  geom_hline(yintercept=0,colour='#AAB7BF',linewidth=.35)+
  geom_line(colour=colour,linewidth=.72,lineend='round')+
  geom_segment(data=hh,aes(x=rank_percent,xend=rank_percent,y=-.335,yend=-.295),inherit.aes=FALSE,colour=colour,alpha=.50,linewidth=.22)+
  annotate('text',x=99,y=.775,label=paste0('NES ',sprintf('%.2f',aa$NES),'\nBH q = ',qstr(aa$BH_q)),hjust=1,vjust=1,size=fs/ggplot2::.pt,family='Arial',lineheight=1.12,colour=ink)+
  scale_x_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),expand=c(0,0))+
  scale_y_continuous(limits=c(-.35,.80),breaks=c(-.2,0,.2,.4,.6),expand=c(0,0))+
  labs(title=labels[ii],x='Ranked genes (%)',y='Running enrichment score')+
  theme_classic(base_size=fs,base_family='Arial')+
  theme(plot.title=element_text(size=fs+1,face='bold',colour=ink,margin=margin(b=9)),
    axis.title=element_text(colour=ink),axis.text=element_text(colour=muted),
    axis.line=element_line(linewidth=.35,colour='#748794'),axis.ticks=element_line(linewidth=.3,colour='#748794'),
    plot.margin=margin(8,12,8,8),axis.title.x=element_text(margin=margin(t=5)),axis.title.y=element_text(margin=margin(r=6)))
}
writeplot <- function(p,name,w,hh) {
 ggsave(file.path(out,paste0(name,'.png')),p,device=ragg::agg_png,width=w,height=hh,units='in',dpi=600,bg='white')
 ggsave(file.path(out,paste0(name,'.pdf')),p,device=grDevices::cairo_pdf,width=w,height=hh,units='in',bg='white')
}
for(j in seq_along(comp)) {
 aa <- a[compartment==comp[j]][1]
 sub <- sprintf('Post-treatment NSCLC  |  Patient-level pseudobulk  |  Degradation High n = %d, Low n = %d',aa$High,aa$Low)
 if(j==3) sub <- paste0(sub,'  |  Exploratory')
 caption <- paste0('Genes ranked by histology-adjusted High − Low contrast: left, High-associated; right, Low-associated.\n',
 'Ticks mark gene-set members. Original NES and BH q retained; ',format(aa$bh_family_n,big.mark=','),' tests in the original BH family.')
 p <- wrap_plots(lapply(1:4,function(i)make(comp[j],i)),ncol=2)+plot_annotation(title=paste0('Four immune programs | ',compnames[j]),subtitle=sub,caption=caption,
 theme=theme(plot.title=element_text(family='Arial',size=18,face='bold',colour=ink),plot.subtitle=element_text(family='Arial',size=9.5,colour=muted,margin=margin(b=8)),plot.caption=element_text(family='Arial',size=8.5,colour=muted,hjust=0,lineheight=1.2),plot.margin=margin(14,14,12,14)))
 writeplot(p,short[j],9.5,7.4)
}
plots <- list()
for(i in 1:4) for(j in 1:3) {
 p <- make(comp[j],i,TRUE)
 if(i==1) p <- p+labs(subtitle=sprintf('%s  |  High %d / Low %d',compnames[j],a[compartment==comp[j]][1]$High,a[compartment==comp[j]][1]$Low))+theme(plot.subtitle=element_text(size=10,face='bold',colour=ink,margin=margin(b=8)))
 plots[[length(plots)+1]] <- p
}
p <- wrap_plots(plots,ncol=3)+plot_annotation(title='Four immune programs across whole tumour, CD8⁺ T cells and DCs',
 subtitle='Same fixed whole-tumour Degradation groups | Post-treatment NSCLC | Patient-level pseudobulk',
 caption='Left of each rank axis: High-associated genes; right: Low-associated genes. Ticks mark gene-set members.\nOriginal NES and BH q retained, including nonsignificant results. BH families: whole tumour 4,958; CD8⁺ 5,959; cDC 3,797 tests.\nWhole tumour and CD8⁺: 12 patients (High 7 / Low 5). cDC: 8 patients (High 5 / Low 3), exploratory.',
 theme=theme(plot.title=element_text(family='Arial',size=18,face='bold',colour=ink),plot.subtitle=element_text(family='Arial',size=10,colour=muted,margin=margin(b=10)),plot.caption=element_text(family='Arial',size=9,colour=muted,hjust=0,lineheight=1.2),plot.margin=margin(16,16,12,16)))
writeplot(p,'P',13,12.2)
writeLines(capture.output(sessionInfo()),file.path(root,'logs','sessionInfo_GSEA_curve_plots.txt'))
cat('Saved W/T/D/P PNG + PDF\n')
