#!/usr/bin/env Rscript
# Author-requested High/Low paired display of the VERIFIED 12 GSEA contrasts.
# Each contrast is shown in two orientations. Reverse NES = -original NES;
# original BH q is copied. These are not 24 tests or group-wise expression values.
# No GSEA rerun, patient regrouping, gene selection or statistical testing.
# Workflow: validate upstream hashes/12 source rows -> expand display -> render.
set.seed(260917)
suppressPackageStartupMessages({library(data.table);library(grid)})
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
out <- dirname(dirname(script)); root <- out
while(!file.exists(file.path(root,'PROJECT_HANDOFF.md'))) {
  stopifnot(dirname(root)!=root); root<-dirname(root)
}
rel <- function(p) substring(normalizePath(p),nchar(root)+2L)
hash <- function(p) digest::digest(file=p,algo='sha256')
options(device=function(...) grDevices::pdf(file=file.path(tempdir(),'sc6_fallback.pdf')))
checks <- list()
check <- function(name,ok) {
  checks[[name]] <<- isTRUE(ok)
  if(!isTRUE(ok)) stop('Validation failed: ',name)
}
previous <- jsonlite::read_json(file.path(out,'results/qa/SC_build.json'),simplifyVector=TRUE)
check('All nine original input hashes unchanged',all(vapply(file.path(root,previous$inputs$path),hash,'')==previous$inputs$sha256))
source <- file.path(out,'results/tables/Fig6_scRNA_GSEA_heatmap.csv')
protected <- c(source,file.path(out,'results/figures',c('SC.png','SC.pdf')))
before <- vapply(protected,hash,'')
s <- fread(source)
ids <- c('REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION','REACTOME_TCR_SIGNALING',
         'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY','HALLMARK_INTERFERON_GAMMA_RESPONSE')
cts <- c('Whole_tumour','CD8','cDC')
labels <- c('Antigen cross-presentation\n(MHC-I)','TCR signaling',
            'TNFR2-related\nnoncanonical NF-κB','IFN-γ response')
check('12 unique source contrasts',nrow(s)==12L && nrow(unique(s[,.(pathway,compartment)]))==12L)
check('Expected programs and compartments',setequal(s$pathway,ids) && setequal(s$compartment,cts))
for(ct in cts) {
  x<-s[compartment==ct]; original<-fread(file.path(root,x$source_table[1]))
  if(ct=='CD8') original<-original[cell_type=='CD8']
  original<-original[match(x$pathway,pathway)]
  check(paste(ct,'source NES/q match'),all(x$NES==original$NES) && all(x$BH_q==original[[x$source_q_field[1]]]))
}
data <- rbindlist(lapply(1:3,function(j) rbindlist(lapply(1:2,function(k) {
  x<-copy(s[compartment==cts[j]])[match(ids,pathway)]
  x[,':='(original_NES=NES,contrast=if(k==1L)'High vs Low' else 'Low vs High',
           numerator=if(k==1L)'High' else 'Low',denominator=if(k==1L)'Low' else 'High',
           multiplier=if(k==1L)1L else -1L,column_index=2L*j-2L+k,
           row_index=seq_len(.N),original_contrast_id=paste(compartment,pathway,sep='::'))]
  x[,NES:=original_NES*multiplier]
  x
}))))
check('24 display cells, 12 original source contrasts',nrow(data)==24L && uniqueN(data$original_contrast_id)==12L)
for(ct in cts) {
  h<-data[compartment==ct & multiplier==1L]; l<-data[compartment==ct & multiplier== -1L]
  check(paste(ct,'mirrored NES'),all(h$NES == -l$NES))
  check(paste(ct,'identical BH q'),all(h$BH_q==l$BH_q))
}
check('No clipping or invalid values',all(is.finite(data$NES)) && all(abs(data$NES)<3.5) && all(data$BH_q>=0 & data$BH_q<=1))
qfmt <- function(q) if(q<.001) formatC(q,format='e',digits=1) else sprintf('%.3f',q)
data[,nes_label:=paste0(sprintf('%+.2f',NES),ifelse(BH_q<.05,'*',''))]
data[,q_label:=paste0('q = ',vapply(BH_q,qfmt,''))]
palette <- colorRamp(c('#326995','#FBFAF8','#B6453E'),space='Lab')
colour <- function(x) rgb(palette((x+3.5)/7)/255)
data[,fill:=colour(NES)]
text_colour <- function(hex) {
  v<-col2rgb(hex)/255
  lin<-ifelse(v<=.04045,v/12.92,((v+.055)/1.055)^2.4)
  l<-as.numeric(c(.2126,.7152,.0722)%*%lin)
  ifelse(1.05/(l+.05)>(l+.05)/.05,'#FFFFFF','#000000')
}
data[,text_colour:=text_colour(fill)]
csv<-file.path(out,'results/tables/Fig6_scRNA_GSEA_two_directions.csv')
fwrite(data,csv)
ink<-'#202C35';muted<-'#64717B';width<-13.2;height<-6.7
left<-.31;right<-.975;gap<-.017
pairwidth<-(right-left-2*gap)/3;cw<-pairwidth/2
pairleft<-left+(0:2)*(pairwidth+gap)
xs<-as.vector(rbind(pairleft+cw/2,pairleft+cw*1.5))
top<-.702;bottom<-.276;rh<-(top-bottom)/4;ys<-top-rh*(.5+0:3)
txt<-function(label,x,y,size=15,font=1,colour=ink,just='centre') {
  grid.text(label,x=x,y=y,just=just,gp=gpar(fontfamily='Arial',fontsize=size,
    fontface=font,col=colour,lineheight=1.08))
}
draw<-function() {
  grid.newpage();grid.rect(gp=gpar(fill='white',col=NA))
  txt('Shared immune programs in NSCLC',.026,.957,24,2,just='left')
  txt('Whole-tumour sarcosine-degradation score groups · Post-treatment NSCLC (GSE207422)',
      .026,.906,13.5,colour=muted,just='left')
  for(j in 1:3) {
    head<-if(j==1L)'Whole-tumour\npseudobulk' else if(j==2L) expression(bold(CD8^'+'~'T cells')) else 'Conventional\nDCs'
    centre<-pairleft[j]+pairwidth/2
    txt(head,centre,.841,17,2)
    z<-s[compartment==cts[j]][1]
    txt(sprintf('High %d / Low %d patients',z$High,z$Low),centre,.782,11.2,colour=muted)
    grid.lines(x=c(pairleft[j]+.005,pairleft[j]+pairwidth-.005),y=c(.760,.760),gp=gpar(col='#D7DEE3',lwd=.8))
  }
  for(j in 1:6) txt(if(j%%2==1)'High vs Low' else 'Low vs High',xs[j],.730,11.8,2)
  for(i in 1:4) {
    txt(labels[i],left-.018,ys[i],16.5,just='right')
    for(j in 1:6) {
      z<-data[row_index==i & column_index==j]
      grid.rect(x=xs[j],y=ys[i],width=cw-.004,height=rh-.007,gp=gpar(fill=z$fill,col=NA))
      txt(z$nes_label,xs[j],ys[i]+.016,17.8,2,colour=z$text_colour)
      txt(z$q_label,xs[j],ys[i]-.022,10.4,colour=z$text_colour)
    }
  }
  txt('NES for the stated contrast',(left+right)/2,.222,13.5)
  nx<-512L;bl<-left+.090;br<-right-.090
  grid.rect(x=bl+(seq_len(nx)-.5)/nx*(br-bl),y=.181,width=(br-bl)/nx+.00003,
    height=.020,gp=gpar(col=NA,fill=colour(seq(-3.5,3.5,length.out=nx))))
  for(v in c(-3.5,0,3.5)) {
    xp<-bl+(v+3.5)/7*(br-bl)
    grid.lines(x=c(xp,xp),y=c(.167,.172),gp=gpar(col=muted,lwd=.6))
    txt(if(v>0)paste0('+',v) else as.character(v),xp,.148,11.5)
  }
  txt('Second group enriched',bl,.114,11.3,colour='#326995',just='left')
  txt('First group enriched',br,.114,11.3,colour='#A63B35',just='right')
  txt('* BH q < 0.05; original correction families retained.',.026,.063,11.3,colour=muted,just='left')
  txt('Paired columns show the same contrast in opposite directions: mirrored NES, identical q.',
      .026,.029,11.1,colour=muted,just='left')
}
png<-file.path(out,'results/figures/SC6.png');pdf<-file.path(out,'results/figures/SC6.pdf')
ragg::agg_png(png,width=width,height=height,units='in',res=600,background='white')
draw();invisible(dev.off())
if(isTRUE(capabilities('aqua'))) {
  grDevices::quartz(type='pdf',file=pdf,width=width,height=height,family='Arial',bg='white')
} else {
  grDevices::cairo_pdf(pdf,width=width,height=height,family='Arial',bg='white')
}
draw();invisible(dev.off())
check('Original three-column files preserved',identical(before,vapply(protected,hash,'')))
check('Original nine inputs preserved',all(vapply(file.path(root,previous$inputs$path),hash,'')==previous$inputs$sha256))
capture.output(sessionInfo(),file=file.path(out,'logs/sessionInfo_SC6.txt'))
receipt<-list(checks=checks,original_inputs=previous$inputs,
  preserved=data.table(path=vapply(protected,rel,''),sha256=unname(before)),
  source_csv=rel(source),display_csv=rel(csv),seed=260917L,
  plot=list(width_inches=width,height_inches=height,dpi=600,x_centres=xs,y_centres=ys,
            cell_width=cw,cell_height=rh,colour_limits=c(-3.5,3.5)),
  statement='12 original comparisons in two orientations; no additional hypothesis tests or group-level expression measurements')
jsonlite::write_json(receipt,file.path(out,'results/qa/SC6_build.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
cat('PASS:',length(checks),'source/display checks. SC6.png and SC6.pdf created.\n')
