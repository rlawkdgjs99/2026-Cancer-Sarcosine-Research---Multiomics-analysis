#!/usr/bin/env Rscript
# Display all species in either criterion in >=1 cohort; no statistical refitting.
# Orange means criterion membership, not organism presence or KO carriage.
suppressPackageStartupMessages({library(data.table);library(grid);library(ragg)})
script <- normalizePath(sub('^--file=', '', grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
OUT <- dirname(script)
CO <- c('PRJNA751792','PRJNA1023797','PRJEB22863')
source <- fread(file.path(OUT,'csv','Venn_3cohorts_membership.csv'))
stopifnot(!anyDuplicated(source[,.(Category,Species_full)]),
          all(source$n_cohorts == rowSums(as.matrix(source[,..CO]))))
PANELS <- list(a=c('Production-associated','NR-enriched'),
               b=c('Degradation-associated','R-enriched'))
ORANGE <- '#EB9438'
plotrows <- rbindlist(lapply(names(PANELS), function(panel){
 cats <- PANELS[[panel]]
 ids <- unique(source[Category %in% cats & n_cohorts>0, Species_full])
 d <- source[Species_full %in% ids & Category %in% cats]
 ranking <- merge(d[Category==cats[1],.(Species_full,Species,n_first=n_cohorts)],
                  d[Category==cats[2],.(Species_full,n_second=n_cohorts)],by='Species_full')
 setorderv(ranking,c('n_first','n_second','Species'),c(-1,-1,1))
 ranking[,Row_order:=seq_len(.N)]
 d <- merge(d,ranking[,.(Species_full,Row_order)],by='Species_full')
 z <- melt(d,id.vars=c('Category','Species_full','Species','Row_order'),
           measure.vars=CO,variable.name='Cohort',value.name='Meets_criterion',variable.factor=FALSE)
 z[,`:=`(Panel=panel,Criterion_order=match(Category,cats),Cohort_order=match(Cohort,CO),
         Species_label=gsub('_',' ',Species,fixed=TRUE),
         Square_fill=ifelse(Meets_criterion,ORANGE,'#FFFFFF'))]
 setorder(z,Row_order,Criterion_order,Cohort_order)
 z
}))
setcolorder(plotrows,c('Panel','Row_order','Species_full','Species','Species_label','Category','Criterion_order','Cohort','Cohort_order','Meets_criterion','Square_fill'))
stopifnot(nrow(plotrows)==102,uniqueN(plotrows[Panel=='a',Species_full])==3,
          uniqueN(plotrows[Panel=='b',Species_full])==14,
          all(plotrows[, .N,by=.(Panel,Species_full)]$N==6))
fwrite(plotrows,file.path(OUT,'csv','Species_4criteria_cohort_membership_plot.csv'),bom=TRUE)
WIDTH <- 9.35
ROW_HEIGHT <- .30
HEADER_HEIGHT <- 1.42
BLOCK_LEFT <- c(4.39,6.77)
BLOCK_WIDTH <- 1.92
CENTERS <- c(.34,.96,1.58)
SQUARE <- .18
FAMILY <- 'Arial'
text_at <- function(label,x,top,size=11,face='plain',just='left',rot=0,color='#111111') {
 grid.text(label,x=unit(x,'in'),y=unit(HEIGHT-top,'in'),just=just,rot=rot,
           gp=gpar(fontfamily=FAMILY,fontsize=size,fontface=face,col=color))
}
rect_at <- function(x,top,w,h,fill='white',line='#333333',lwd=.8) {
 grid.rect(x=unit(x,'in'),y=unit(HEIGHT-top,'in'),width=unit(w,'in'),height=unit(h,'in'),
           just=c('left','top'),gp=gpar(fill=fill,col=line,lwd=lwd))
}
panel_height <- function(panel) HEADER_HEIGHT + uniqueN(plotrows[Panel==panel,Species_full])*ROW_HEIGHT + .06
draw_panel <- function(panel,top){
 d <- plotrows[Panel==panel];cats <- PANELS[[panel]]
 n <- uniqueN(d$Species_full);rowtop <- top+HEADER_HEIGHT
 text_at(panel,.28,top+.06,18,'bold')
 text_at(paste(n,'species'),.28,top+.34,10.5,color='#555555')
 labels <- unique(d[,.(Row_order,Species_label)])
 for (i in 1:2){
  x0 <- BLOCK_LEFT[i]
  title <- if(cats[i]=='NR-enriched') 'ICI-NR-enriched' else if(cats[i]=='R-enriched') 'ICI-R-enriched' else cats[i]
  nspecies <- uniqueN(d[Category==cats[i] & Meets_criterion==TRUE,Species_full])
  text_at(title,x0+BLOCK_WIDTH/2,top+.05,13,'bold',just='centre')
  text_at(paste0('(',nspecies,' species)'),x0+BLOCK_WIDTH/2,top+.29,10,color='#555555',just='centre')
  for(j in 1:3) text_at(CO[j],x0+CENTERS[j],top+.89,9.5,rot=55,just='centre')
  rect_at(x0,rowtop-.02,BLOCK_WIDTH,n*ROW_HEIGHT+.04,fill=NA,line='#111111',lwd=1)
 }
 for(j in seq_len(n)){
  label <- labels[Row_order==j,Species_label]
  face <- if(startsWith(label,'GGB')) 'plain' else 'italic'
  text_at(label,4.16,rowtop+(j-.5)*ROW_HEIGHT,11.4,face,just='right')
 }
 for(i in seq_len(nrow(d))){
  x <- BLOCK_LEFT[d$Criterion_order[i]]+CENTERS[d$Cohort_order[i]]
  y <- rowtop+(d$Row_order[i]-.5)*ROW_HEIGHT
  rect_at(x-SQUARE/2,y-SQUARE/2,SQUARE,SQUARE,fill=d$Square_fill[i],line='#333333',lwd=.85)
 }
 invisible(top+panel_height(panel))
}
draw_footer <- function(top){
 rect_at(.28,top+.015,.14,.14,ORANGE)
 text_at('Meets criterion',.49,top+.085,10)
 rect_at(2.08,top+.015,.14,.14,'white')
 text_at('Does not meet criterion',2.29,top+.085,10)
 text_at('Species shown: union across all three cohorts and the two criteria in each panel.',.28,top+.35,8.4,color='#444444')
 text_at('Association: Spearman rho > 0.3 and BH q < 0.05 with the 2-KO production / 5-KO degradation score.',.28,top+.53,8.4,color='#444444')
 text_at('Enrichment: full-species DA BH q < 0.05 and log2FC(NR/R) > 1 (NR) or < -1 (R).',.28,top+.71,8.4,color='#444444')
 text_at('Species prevalence >=10% within cohort; white squares do not imply absence of the organism.',.28,top+.89,8.4,color='#444444')
}
render <- function(panels,name){
 HEIGHT <<- .82 + sum(vapply(panels,panel_height,numeric(1))) + .48*(length(panels)-1) + 1.14
 draw <- function(){
  grid.newpage()
  text_at('NSCLC species membership by cohort',.28,.26,16,'bold')
  text_at('All qualifying species included, regardless of overlap between cohorts',.28,.56,10.5,color='#555555')
  top <- .88
  for(panel in panels){top <- draw_panel(panel,top);if(panel!=tail(panels,1))top<-top+.48}
  draw_footer(top+.10)
 }
 png <- file.path(OUT,'figures',paste0(name,'.png'))
 pdf <- file.path(OUT,'figures',paste0(name,'.pdf'))
 ragg::agg_png(png,width=WIDTH,height=HEIGHT,units='in',res=300,background='white');draw();dev.off()
 grDevices::cairo_pdf(pdf,width=WIDTH,height=HEIGHT,family=FAMILY,onefile=TRUE);draw();dev.off()
 cat(name,WIDTH,HEIGHT,'inches\n')
}
render(c('a','b'),'Main_species_4criteria_cohort_membership')
render('a','Main_species_Production_NR_cohort_membership')
render('b','Main_species_Degradation_R_cohort_membership')
cat('Plotted square counts by criterion/cohort:\n')
print(plotrows[,.(filled=sum(Meets_criterion),cells=.N),by=.(Panel,Category,Cohort)])
