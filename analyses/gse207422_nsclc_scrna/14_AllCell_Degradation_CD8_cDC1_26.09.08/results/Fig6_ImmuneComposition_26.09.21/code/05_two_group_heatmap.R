#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(grid);library(digest)})
arg<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]);base<-normalizePath(file.path(dirname(arg),'..'))
v<-fread(file.path(base,'tables/patient_values.csv'));s<-fread(file.path(base,'tables/statistics.csv'));plan<-fromJSON(file.path(base,'ANALYSIS_PLAN.json'));cells<-plan$targets
inputs<-file.path(base,c('tables/patient_values.csv','tables/statistics.csv','ANALYSIS_PLAN.json','figures/H.png','figures/H.pdf'));hash<-vapply(inputs,function(p)digest(file=p,algo='sha256',serialize=FALSE),character(1))
stopifnot(nrow(v)==156,!anyDuplicated(v[,.(Patient,cell_type)]),uniqueN(v$Patient)==12,all(s$n_Low==5),all(s$n_High==7))
g<-v[,.(n=.N,mean_percent=mean(percent),mean_patient_z=mean(row_z)),by=.(cell_type,group)];g[,row_order:=match(cell_type,cells)];g[,col_order:=match(group,c('Low','High'))];setorder(g,row_order,col_order)
for(ct in cells){d<-v[cell_type==ct];z<-(d$percent-mean(d$percent))/sd(d$percent);stopifnot(max(abs(z-d$row_z))<1e-12);for(gr in c('Low','High'))stopifnot(abs(g[cell_type==ct&group==gr,mean_percent]-s[cell_type==ct,get(paste0('mean_',gr,'_percent'))])<1e-10)}
g<-merge(g,s[,.(cell_type,BH13_q)],by='cell_type',sort=FALSE);setorder(g,row_order,col_order)
fwrite(g,file.path(base,'tables/two_group_heatmap.csv'),bom=TRUE)
labels<-c('CD8+ T cells','CD4+ T cells','Cycling T cells','NK cells','NK / γδ T cells\n(unresolved)','B cells','Plasma cells','Conventional DCs','pDCs','Monocytes','Macrophages','Neutrophils','Mast cells')
blue<-'#2C6CA2';red<-'#BA4249';ink<-'#243842';muted<-'#637783';low<-'#2E8798';high<-'#C3466C';lim<-max(1,ceiling(max(abs(g$mean_patient_z))*4)/4)
pal<-colorRamp(c(blue,'#FAFAF8',red));colz<-function(z)rgb(pal(pmax(0,pmin(1,(z+lim)/(2*lim)))),maxColorValue=255)
txt<-function(label,x,y,size=12,col=ink,bold=FALSE,just='centre'){grid.text(label,x,y,just=just,gp=gpar(fontfamily='Arial',fontsize=size,col=col,fontface=if(bold)'bold'else'plain',lineheight=1.04))}
rect<-function(x,y,w,h,fill,border=NA,lwd=.5)grid.rect(x,y,w,h,gp=gpar(fill=fill,col=border,lwd=lwd))
draw<-function(){
 grid.newpage();txt('Immune-cell composition',.045,.955,21,bold=TRUE,just='left');txt('Whole-tumour Degradation Low vs High',.045,.915,12.5,muted,just='left')
 txt('Post-treatment NSCLC | Patient-level group means',.045,.885,10.5,muted,just='left')
 x0<-.395;tw<-.177;ytop<-.781;th<-.0425
 txt('Low',x0+.5*tw,.839,15,low,TRUE);txt('n = 5',x0+.5*tw,.812,11,muted)
 txt('High',x0+1.5*tw,.839,15,high,TRUE);txt('n = 7',x0+1.5*tw,.812,11,muted)
 txt('BH q',.882,.839,13,bold=TRUE);txt('13 cell types',.882,.812,10,muted)
 for(i in seq_along(cells)){
  ct<-cells[i];y<-ytop-(i-.5)*th;st<-s[cell_type==ct]
  txt(labels[i],x0-.018,y,12.5,bold=ct%in%c('CD8 T cell','Conventional DC'),just='right')
  for(j in 1:2){row<-g[cell_type==ct&col_order==j];zz<-row$mean_patient_z;rect(x0+(j-.5)*tw,y,tw,th,colz(zz),'white',1.2);txt(sprintf('%.2f%%',row$mean_percent),x0+(j-.5)*tw,y,13,col=if(abs(zz)/lim>.70)'white'else ink,bold=TRUE)}
  txt(sprintf('%.3f',st$BH13_q),.882,y,12.5,col=if(st$BH13_q<.05)ink else muted,bold=st$BH13_q<.05)
 }
 txt('Colour: mean patient-level row z score',.395,.189,10.5,just='left')
 for(k in 1:256)rect(.395+(k-.5)/256*.30,.164,.30/256+.0001,.017,colz(-lim+(k-.5)/256*2*lim))
 for(z in c(-lim,0,lim))txt(if(z>0)paste0('+',z)else as.character(z),.395+.30*(z+lim)/(2*lim),.138,10,muted)
 txt('Numbers: mean fraction of all captured singlets; each patient weighted equally.',.045,.099,9.5,muted,just='left')
 txt('Colours average the original 12-patient z scores; no two-column restandardization.',.045,.074,9.5,muted,just='left')
 txt('Original histology-stratified permutation tests; BH across 13 types. Bold q < 0.05.',.045,.049,9.5,muted,just='left')
 txt('Relative captured-cell composition; not absolute tissue abundance.',.045,.024,9.5,muted,just='left')
}
cairo_pdf(file.path(base,'figures/H2.pdf'),width=8.5,height=8.6,family='Arial',bg='white');draw();dev.off()
ragg::agg_png(file.path(base,'figures/H2.png'),width=8.5,height=8.6,units='in',res=600,background='white');draw();dev.off()
ragg::agg_png(file.path(base,'qa/H2_preview.png'),width=8.5,height=8.6,units='in',res=140,background='white');draw();dev.off()
stopifnot(all(hash==vapply(inputs,function(p)digest(file=p,algo='sha256',serialize=FALSE),character(1))))
write_json(list(type='Display-only two-group summary',inputs=data.frame(path=inputs,sha256=hash),color_scale=c(-lim,lim),group_n=c(Low=5,High=7),statistical_reanalysis=FALSE,checks='26 arithmetic group mean percentages match original; original row-z formula matches; all source files unchanged',visual_QA='pending'),file.path(base,'qa/H2_QA.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Two-column plot complete; source values and original figures unchanged\n')
