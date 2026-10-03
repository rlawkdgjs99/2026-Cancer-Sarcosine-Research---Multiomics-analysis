#!/usr/bin/env Rscript
# Re-display the existing Figure 6 e/g/i GSEA results. No inferential rerun.
# Design: post-treatment patients, frozen whole-tumour degradation groups,
# histology-adjusted patient pseudobulk contrasts. Original BH families retained.
# Workflow: validate exact IDs/rosters/full results -> export 12 source rows -> plot.
# Deterministic figure; seed is recorded for reproducibility, not used for inference.
set.seed(260917)
cat('Loading plotting dependencies.\n')
suppressPackageStartupMessages({library(data.table); library(grid)})
arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
script <- normalizePath(sub('^--file=', '', arg))
out <- dirname(dirname(script))
base <- dirname(out)
root <- out
while (!file.exists(file.path(root, 'PROJECT_HANDOFF.md'))) {
  stopifnot(dirname(root) != root)
  root <- dirname(root)
}
rel <- function(p) substring(normalizePath(p), nchar(root)+2L)
options(device=function(...) grDevices::pdf(file=file.path(tempdir(), 'sc_heatmap_fallback.pdf')))
dirs <- c(Whole_tumour=out,
          CD8=file.path(base, '14_AllCell_Degradation_CD8_cDC1_26.09.08'),
          cDC=file.path(base, '17_cDC_AllCell_Degradation_26.09.16'))
ids <- c('REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION', 'REACTOME_TCR_SIGNALING',
         'REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY', 'HALLMARK_INTERFERON_GAMMA_RESPONSE')
labels <- c('Antigen cross-presentation\n(MHC-I)', 'TCR signaling',
            'TNFR2-related\nnoncanonical NF-κB', 'IFN-γ response')
source_files <- c(file.path(dirs[1], 'results/tables/plotdata_Fig8e_Exact_Fig7_Pathways.csv'),
  file.path(dirs[2], 'results/Fig8_CellType_Pathways_26.09.14/results/tables/plotdata_Exact_Fig7_Pathways.csv'),
  file.path(dirs[3], 'results/tables/plotdata_Exact_display.csv'))
full_files <- file.path(dirs, 'results/tables/02_GSEA_all_scopes.csv')
roster_files <- c(file.path(dirs[1], 'results/tables/01_frozen_patient_groups.csv'),
  file.path(dirs[2], 'results/tables/01_FROZEN_allcell_patient_scores_groups.csv'),
  file.path(dirs[3], 'results/tables/01_patient_eligibility.csv'))
q_fields <- c('q_global','q_global','q_global_cDC')
family_n <- c(4958L, 5959L, 3797L)
family_labels <- c('All whole-tumour gene sets', 'CD8 and cDC1 gene-set tests',
                   'All conventional-DC gene sets')
sample_n <- c(12L,12L,8L); high_n <- c(7L,7L,5L); low_n <- c(5L,5L,3L)
checks <- list()
check <- function(name, ok) {
  checks[[name]] <<- isTRUE(ok)
  if (!isTRUE(ok)) stop('Validation failed: ', name)
}
input_files <- c(source_files, full_files, roster_files)
hashes <- rbindlist(lapply(input_files,function(p) {
  cat('Reading input:',basename(dirname(dirname(p))),basename(p),'\n'); flush.console()
  data.table(path=rel(p),sha256=digest::digest(file=p,algo='sha256'))
}))
cat('Source hashes recorded; reading frozen rosters.\n')
rosters <- lapply(roster_files, fread)
w <- rosters[[1]]; a <- rosters[[2]]; cdc <- rosters[[3]]
check('Unique whole-tumour patients', nrow(w)==15L && !anyDuplicated(w$Patient))
check('Whole-tumour and CD8 frozen labels identical',
      identical(w[order(Patient),.(Patient,group,degradation_mean_z,cutoff)],
                a[order(Patient),.(Patient,group,degradation_mean_z,cutoff)]))
check('Score is the original mean of SARDH/PIPOX gene z values',
      max(abs(w$degradation_mean_z-(w$SARDH_z+w$PIPOX_z)/2)) < 1e-12)
check('Original median cutpoint and group labels',
      max(abs(w$cutoff-median(w$degradation_mean_z))) < 1e-12 &&
      all(w$group==ifelse(w$degradation_mean_z>w$cutoff,'High','Low')))
post <- w[treatment=='Post']
cp <- cdc[treatment=='Post' & toupper(as.character(eligible))=='TRUE']
check('Whole-tumour and CD8 post-treatment roster',
      nrow(post)==12L && sum(post$group=='High')==7L && sum(post$group=='Low')==5L && all(post$CD8_cells>=50))
check('Primary cell counts', sum(post$all_cells)==78192L && sum(post$CD8_cells)==16229L)
check('Conventional DC eligibility and cell counts',
      nrow(cp)==8L && sum(cp$group=='High')==5L && sum(cp$group=='Low')==3L &&
      all(cp$cDC_cells>=50) && sum(cp$cDC_cells)==772L)
check('Conventional DC reuses whole-tumour group labels',
      all(cp$group==w$group[match(cp$Patient,w$Patient)]))
parts <- list()
for (j in seq_along(dirs)) {
  cat('Validating compartment',j,'\n')
  d <- fread(source_files[j]); f <- fread(full_files[j])[scope=='post_group']
  check(paste(j,'original BH family size'), nrow(f)==family_n[j])
  # Independent diagnostic only: reproduce original BH over the COMPLETE family.
  # Do not replace saved q values, and never adjust only the 12 displayed cells.
  check(paste(j,'original full-family BH values'),
        max(abs(p.adjust(f$pval,method='BH')-f[[q_fields[j]]])) < 1e-12)
  if (j==2L) {d<-d[cell_type=='CD8']; f<-f[cell_type=='CD8']}
  check(paste(j,'four exact gene sets'), nrow(d)==4L && !anyDuplicated(d$pathway) && setequal(d$pathway,ids))
  d <- d[match(ids,pathway)]; f <- f[match(ids,pathway)]
  check(paste(j,'NES and q equal full primary result'),
        all(d$NES==f$NES) && all(d[[q_fields[j]]]==f[[q_fields[j]]]))
  check(paste(j,'n and group counts'), all(d$n==sample_n[j]) && all(d$High==high_n[j]) && all(d$Low==low_n[j]))
  check(paste(j,'finite NES and valid q'), all(is.finite(d$NES)) && all(is.finite(d[[q_fields[j]]])) &&
        all(d[[q_fields[j]]]>=0 & d[[q_fields[j]]]<=1))
  if ('scope' %in% names(d)) check(paste(j,'primary scope'),all(d$scope=='post_group'))
  parts[[j]] <- data.table(pathway=ids, program=gsub('\n',' ',labels),
    compartment=names(dirs)[j], original_panel=c('6e','6g','6i')[j],
    scope='post_group', NES=d$NES, BH_q=d[[q_fields[j]]], n=d$n, High=d$High, Low=d$Low,
    cells=c(78192L,16229L,772L)[j], bh_family_n=family_n[j], bh_family=family_labels[j],
    source_table=rel(source_files[j]), source_q_field=q_fields[j],
    eligible_patients=paste(sort(if(j==3L) cp$Patient else post$Patient),collapse=';'))
}
data <- rbindlist(parts)
check('12 unique pathway by compartment cells', nrow(unique(data[,.(pathway,compartment)]))==12L)
check('Signed NES not clipped', all(abs(data$NES)<3.5))
qfmt <- function(q) if(q<.001) formatC(q,format='e',digits=1) else sprintf('%.3f',q)
data[,nes_label:=paste0(sprintf('%+.2f',NES),ifelse(BH_q<.05,'*',''))]
data[,q_label:=paste0('q = ',vapply(BH_q,qfmt,''))]
# Shared, zero-centred diverging scale. Positive = degradation-High enrichment.
palette <- colorRamp(c('#326995','#FBFAF8','#B6453E'),space='Lab')
colour <- function(x) rgb(palette((x+3.5)/7)/255)
data[,fill:=colour(NES)]
text_colour <- function(hex) {
  s <- col2rgb(hex)/255
  lin <- ifelse(s<=.04045,s/12.92,((s+.055)/1.055)^2.4)
  l <- as.numeric(c(.2126,.7152,.0722)%*%lin)
  dark_l <- ((16/255+.055)/1.055)^2.4
  ifelse(1.05/(l+.05) > (l+.05)/(dark_l+.05),'#FFFFFF','#101010')
}
data[,text_colour:=text_colour(fill)]
csv <- file.path(out,'results/tables/Fig6_scRNA_GSEA_heatmap.csv')
fwrite(data,csv)
cat('Source table exported; rendering PNG and PDF.\n')
ink <- '#202C35'; muted <- '#64717B'; width <- 10.4; height <- 6.7
left <- .405; right <- .970; top <- .735; bottom <- .285
cw <- (right-left)/3; rh <- (top-bottom)/4
xs <- left+cw*(.5+0:2); ys <- top-rh*(.5+0:3)
txt <- function(label,x,y,size=15,font=1,colour=ink,just='centre',...) {
  grid.text(label,x=x,y=y,just=just,gp=gpar(fontfamily='Arial',fontsize=size,
    fontface=font,col=colour,lineheight=1.08),...)
}
draw <- function() {
  grid.newpage(); grid.rect(gp=gpar(fill='white',col=NA))
  txt('Shared immune programs in NSCLC',.035,.952,23,2,just='left')
  txt('Post-treatment tumour scRNA-seq (GSE207422)',.035,.899,13.5,colour=muted,just='left')
  headers <- c('Whole-tumour\npseudobulk','CD8⁺ T cells','Conventional\nDCs')
  for(j in 1:3) {
    # Plotmath uses an ordinary '+' glyph at superscript position; Unicode U+207A
    # is not reliably embedded by the macOS Quartz PDF device with Arial.
    header <- if(j==2L) expression(bold(CD8^'+'~'T cells')) else headers[j]
    txt(header,xs[j],.822,16.5,2)
    txt(sprintf('High %d / Low %d',high_n[j],low_n[j]),xs[j],.764,11.3,colour=muted)
  }
  for(i in 1:4) {
    txt(labels[i],left-.022,ys[i],16.4,just='right')
    for(j in 1:3) {
      z <- data[pathway==ids[i] & compartment==names(dirs)[j]]
      grid.rect(x=xs[j],y=ys[i],width=cw-.005,height=rh-.007,
                gp=gpar(fill=z$fill,col=NA))
      txt(z$nes_label,xs[j],ys[i]+.016,19,2,colour=z$text_colour)
      txt(z$q_label,xs[j],ys[i]-.022,11.8,colour=z$text_colour)
    }
  }
  txt('NES (Degradation-High vs Low)',(left+right)/2,.227,13.5)
  nx <- 512L; barleft <- left+.055; barright <- right-.055
  cols <- colour(seq(-3.5,3.5,length.out=nx))
  grid.rect(x=barleft+(seq_len(nx)-.5)/nx*(barright-barleft),y=.185,
    width=(barright-barleft)/nx+.00003,height=.019,
    gp=gpar(col=NA,fill=cols))
  for(v in c(-3.5,0,3.5)) {
    xp <- barleft+(v+3.5)/7*(barright-barleft)
    grid.lines(x=c(xp,xp),y=c(.172,.178),gp=gpar(col=muted,lwd=.6))
    txt(if(v>0) paste0('+',v) else as.character(v),xp,.153,11.5)
  }
  txt('Low-enriched',barleft,.12,11.5,colour='#326995',just='left')
  txt('High-enriched',barright,.12,11.5,colour='#A63B35',just='right')
  txt('* BH q < 0.05; original correction families retained.',.035,.064,11.3,colour=muted,just='left')
  txt('Patient-level pseudobulk; histology-adjusted GSEA. Column counts are patients.',
      .035,.029,10.5,colour=muted,just='left')
}
png <- file.path(out,'results/figures/SC.png')
pdf <- file.path(out,'results/figures/SC.pdf')
ragg::agg_png(png,width=width,height=height,units='in',res=600,background='white')
draw(); invisible(dev.off())
if (isTRUE(capabilities('aqua'))) {
  # Native macOS vector PDF: avoids an optional XQuartz-dependent Cairo library.
  grDevices::quartz(type='pdf',file=pdf,width=width,height=height,family='Arial',bg='white')
} else {
  grDevices::cairo_pdf(pdf,width=width,height=height,family='Arial',bg='white')
}
draw(); invisible(dev.off())
capture.output(sessionInfo(),file=file.path(out,'logs/sessionInfo_scRNA_heatmap.txt'))
hashes[,sha256_after:=vapply(file.path(root,path),
  function(p) digest::digest(file=p,algo='sha256'),'')]
if (any(hashes$sha256_after!=hashes$sha256)) print(hashes[sha256_after!=sha256])
check('All source inputs unchanged',all(hashes$sha256_after==hashes$sha256))
receipt <- list(task='Existing Figure 6 e/g/i GSEA re-display', checks=checks,inputs=hashes,
  packages=list(R=R.version.string,ragg=as.character(packageVersion('ragg')),
                data_table=as.character(packageVersion('data.table'))),
  plot=list(width_inches=width,height_inches=height,dpi=600,limits=c(-3.5,3.5),
            x_centres=xs,y_centres=ys),source_table=rel(csv),
  outputs=c(rel(png),rel(pdf)),seed=260917L)
jsonlite::write_json(receipt,file.path(out,'results/qa/SC_build.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
cat('PASS:',length(checks),'source/design checks.\n',rel(png),'\n',rel(pdf),'\n')
