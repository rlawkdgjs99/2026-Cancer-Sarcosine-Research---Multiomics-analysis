#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(grid)})
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);base<-normalizePath(file.path(dirname(script),'..'))
plan<-fromJSON(file.path(base,'ANALYSIS_PLAN.json'));v<-fread(file.path(base,'tables/patient_values.csv'));s<-fread(file.path(base,'tables/statistics.csv'))
m<-fread(file.path(base,'tables/heatmap_column_annotations.csv'));patients<-m$Patient;cells<-plan$targets
stopifnot(length(patients)==12L,identical(m$group,c(rep('Low',5),rep('High',7))))
labels<-c('CD8⁺ T cells','CD4⁺ T cells','Cycling T cells','NK cells','NK / γδ T cells\n(unresolved)','B cells','Plasma cells','Conventional DCs','pDCs','Monocytes','Macrophages','Neutrophils','Mast cells')
blue<-'#2C6CA2';red<-'#BA4249';low<-'#2E8798';high<-'#C3466C';ink<-'#243842';muted<-'#637783';hc<-c(Adeno='#677DA5',Squamous='#C3A269')
pal<-colorRamp(c(blue,'#FAFAF8',red));colz<-function(z){x<-pmax(0,pmin(1,(z+2.5)/5));rgb(pal(x),maxColorValue=255)}
txt<-function(label,x,y,size=10,col=ink,bold=FALSE,just='centre'){grid.text(label,x,y,just=just,gp=gpar(fontfamily='Arial',fontsize=size,col=col,fontface=if(bold)'bold' else 'plain',lineheight=.98))}
rect<-function(x,y,w,h,fill,border=NA,lwd=.5){grid.rect(x,y,w,h,gp=gpar(fill=fill,col=border,lwd=lwd))}
line<-function(x,y,col,lwd=1){grid.lines(x,y,gp=gpar(col=col,lwd=lwd))}
draw<-function(){
 grid.newpage();txt('Immune-cell composition',.035,.955,19,bold=TRUE,just='left')
 txt('Post-treatment NSCLC  |  12 patients  |  Frozen whole-tumour degradation groups',.035,.918,10.5,muted,just='left')
 x0<-.20;hw<-.435;ytop<-.78;hh<-.555;tw<-hw/12;th<-hh/13
 for(i in seq_along(cells)){
  cell<-cells[i];y<-ytop-(i-.5)*th;txt(labels[i],x0-.012,y,10.5,bold=cell %in% c('CD8 T cell','Conventional DC'),just='right')
  d<-v[cell_type==cell][match(patients,Patient)]
  for(j in seq_along(patients))rect(x0+(j-.5)*tw,y,tw,th,colz(d$row_z[j]),'white',.7)
  st<-s[cell_type==cell];txt(sprintf('%.2f',st$mean_Low_percent),.682,y,10);txt(sprintf('%.2f',st$mean_High_percent),.734,y,10)
  xx<-.846+st$difference_pp/36*.122;cc<-if(st$difference_pp>0)red else blue
  line(c(.846,xx),c(y,y),cc,1.3);grid.points(unit(xx,'npc'),unit(y,'npc'),pch=16,size=unit(2.1,'mm'),gp=gpar(col=cc))
  txt(sprintf('%.3f',st$BH13_q),.956,y,10,col=if(st$BH13_q<.05)ink else muted,bold=st$BH13_q<.05)
 }
 rect(x0+5*tw,ytop-hh/2,.004,hh,'white')
 # zero guide is drawn as small segments between rows to keep all points visible
 for(i in 0:12){yy<-ytop-i*th;line(c(.846,.846),c(yy-th*.25,yy-th*.75),'#B7C2C9',.65)}
 for(j in seq_along(patients)){
  x<-x0+(j-.5)*tw
  txt(patients[j],x,.203,9)
  rect(x,.816,tw,.016,if(m$group[j]=='Low')low else high)
  rect(x,.797,tw,.016,hc[m$histology[j]])
 }
 rect(x0+5*tw,.8065,.004,.036,'white')
 txt('Degradation',x0-.012,.816,8.8,just='right');txt('Histology',x0-.012,.797,8.8,just='right')
 txt('Low (n = 5)',x0+2.5*tw,.853,11,low,TRUE);txt('High (n = 7)',x0+8.5*tw,.853,11,high,TRUE)
 txt('Mean fraction (%)',.708,.853,10,bold=TRUE);txt('Low',.682,.812,9.5,low);txt('High',.734,.812,9.5,high)
 txt('High − Low',.846,.853,10,bold=TRUE);txt('percentage points',.846,.812,8.7,muted)
 txt('BH q',.956,.853,10,bold=TRUE);txt('13 tests',.956,.812,8.5,muted)
 line(c(.785,.907),c(.225,.225),'#ABB5BC',.6)
 for(z in c(-15,0,15)){x<-.846+z/36*.122;line(c(x,x),c(.225,.22),'#ABB5BC',.6);txt(if(z>0)paste0('+',z)else as.character(z),x,.203,8.5,muted)}
 txt('Relative fraction within each cell type (row z score)',x0,.165,9.5,just='left')
 cbx<-seq(x0,x0+.19,length.out=257)
 for(k in 1:256)rect((cbx[k]+cbx[k+1])/2,.136,.19/256+.0001,.015,colz(-2.5+(k-.5)/256*5))
 for(z in c(-2.5,0,2.5))txt(if(z>0)paste0('+',z)else as.character(z),x0+.19*(z+2.5)/5,.112,8.5,muted)
 for(j in seq_along(hc)){x<-.49+(j-1)*.10;rect(x,.136,.012,.015,hc[j]);txt(names(hc)[j],x+.012,.136,9,just='left')}
 txt('Fractions use all captured singlets; group means weight each patient equally. Colours saturate beyond ±2.5 z.',.035,.077,9,muted,just='left')
 txt('BH q: histology-stratified permutation tests on log-ratios (300 allocations), adjusted across 13 cell types.',.035,.052,9,muted,just='left')
 txt('Group differences are descriptive percentage points; bold q < 0.05. Relative capture composition, not absolute tissue abundance.',.035,.027,8.6,muted,just='left')
}
cairo_pdf(file.path(base,'figures/H.pdf'),width=12.5,height=7.9,family='Arial',bg='white');draw();dev.off()
ragg::agg_png(file.path(base,'figures/H.png'),width=12.5,height=7.9,units='in',res=600,background='white');draw();dev.off()
ragg::agg_png(file.path(base,'qa/H_preview.png'),width=12.5,height=7.9,units='in',res=150,background='white');draw();dev.off()
