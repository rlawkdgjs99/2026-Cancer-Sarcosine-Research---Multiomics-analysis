# Recreate CRC Figure 3a/b from frozen criterion membership, 2026-09-11.
# No abundance detection, statistical refitting or new species selection.
# Row order is extracted from slide4 image32/image33 in the author deck:
# Sarcosine FIgures 구상 - 26.09.11.pptx, SHA256
# 13af5bd49ee6b0c36a4df050ba8addb5eefd1f6ae08caea972ca3ede5a9e9313.
# Run: Rscript plot_Fig3ab_cohort_membership_26.09.11.R
# Optional positional arguments: integrated analysis directory, output directory.
args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- sub("^--file=", "", script_arg[1])
analysis_dir <- if(length(args)>=1) normalizePath(args[1]) else dirname(dirname(normalizePath(script)))
crc_dir <- dirname(analysis_dir)
output_dir <- if(length(args)>=2) args[2] else file.path(analysis_dir,"results_integrated","CRC_WGS_Fig3ab_membership_26.09.11")
source <- Sys.glob(file.path(crc_dir,"*","CRC_WGS_4_Venn_species_membership_26.09.01","CRC_WGS_4_Venn_species_membership.csv"))
stopifnot(length(source)==1L)
d <- read.csv(source, check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM")
stopifnot(nrow(d)==516L, !anyDuplicated(paste(d$Venn_category,d$Species)))
cohorts <- c("PRJEB6070","PRJEB10878","PRJEB27928","PRJNA429097")
# The original Windows column order is not asserted: the new output explicitly
# labels this consistent order, matching other current CRC manuscript panels.
colors <- c("Production-associated"="#2166AC", "Cancer-enriched"="#D73027",
            "Degradation-associated"="#F28E2B", "Healthy-enriched"="#70BDE7")
titles <- c("Production-associated"="Production-associated", "Cancer-enriched"="CRC-enriched",
            "Degradation-associated"="Degradation-associated", "Healthy-enriched"="Healthy-enriched")
row_order <- list(
  a=c("Hungatella_hathewayi", "Clostridium_symbiosum", "Ruthenibacterium_lactatiformans", "Eisenbergiella_tayi", "Enterocloster_aldenensis", "Enterocloster_lavalensis", "Erysipelatoclostridium_ramosum", "GGB33512_SGB15201", "Anaeromassilibacillus_sp_An250", "Anaerotruncus_rubiinfantis", "Candidatus_Pararuminococcus_gallinarum", "Clostridiales_bacterium", "Clostridium_SGB4750", "Dysosmobacter_sp_NSJ_60", "Enterocloster_bolteae", "Enterocloster_citroniae", "GGB2653_SGB3574", "GGB2982_SGB3964", "GGB33586_SGB53517", "GGB34900_SGB14891", "GGB9818_SGB15459", "Intestinimonas_butyriciproducens"),
  b=c("Clostridiales_bacterium_KLE1615", "Coprococcus_eutactus", "Lachnospira_eligens", "Anaerobutyricum_hallii", "Clostridium_sp_AF27_2AA", "Faecalibacillus_intestinalis", "Fusicatenibacter_saccharivorans", "GGB51441_SGB71759", "Lachnospiraceae_bacterium_AM48_27BH", "Lachnospiraceae_bacterium_Marseille_Q4251", "Roseburia_faecis", "Adlercreutzia_equolifaciens", "Anaerostipes_hadrus", "Blautia_stercoris", "Blautia_wexlerae", "Clostridiaceae_bacterium_AF18_31LB", "Clostridiaceae_bacterium_Marseille_Q4143", "Clostridium_sp_AF34_13", "GGB79734_SGB15291", "Oscillospiraceae_bacterium_Marseille_Q3528", "Roseburia_inulinivorans", "Agathobaculum_butyriciproducens", "Blautia_glucerasea", "Blautia_massiliensis", "Clostridiaceae_unclassified_SGB4771", "Clostridium_sp_AF34_10BH", "Clostridium_sp_AM33_3", "Desulfovibrio_SGB5077", "Dysosmobacter_sp_BX15", "Eubacterium_ramulus", "Eubacterium_rectale", "Eubacterium_sp_AF34_35BH", "Faecalicatena_fissicatena", "GGB3746_SGB5089", "GGB9614_SGB15049", "Lachnospira_sp_NSJ_43", "Lachnospiraceae_bacterium", "Roseburia_intestinalis", "Ruminococcus_bromii")
)
group_pairs <- list(a=c("Production-associated","Cancer-enriched"),
                    b=c("Degradation-associated","Healthy-enriched"))
# Frozen input table contains the full union for each criterion. A species
# absent from a criterion's union is non-member in all four cohorts.
make_flags <- function(species, group) {
  r <- d[d$Species==species & d$Venn_category==group,,drop=FALSE]
  stopifnot(nrow(r)<=1L)
  if(nrow(r)==0L) return(rep(FALSE,4))
  flags <- as.logical(unlist(r[paste0(cohorts,"_in_set")],use.names=FALSE))
  stopifnot(!anyNA(flags),sum(flags)==r$n_cohorts_in_set)
  flags
}
# Independent checks against the frozen full per-criterion membership tables.
source_rel <- c("Production-associated"="cross_cohort/cross_cohort_prod_assoc_membership_matrix.csv",
 "Cancer-enriched"="cross_cohort/cross_cohort_CRC_enriched_membership_matrix.csv",
 "Degradation-associated"="cross_cohort/cross_cohort_deg_assoc_membership_matrix.csv",
 "Healthy-enriched"="SupFig3de_HEALTHY_ENRICHED_26.08.27/healthy_enriched_membership_matrix.csv")
for(group in names(source_rel)) {
  raw <- read.csv(file.path(analysis_dir,"results_integrated",source_rel[[group]]),check.names=FALSE)
  m <- d[d$Venn_category==group,,drop=FALSE]
  stopifnot(setequal(raw$Species,m$Species),!anyDuplicated(raw$Species))
  raw <- raw[match(m$Species,raw$Species),,drop=FALSE]
  suffix <- if(group=="Healthy-enriched") "_strict" else "_sig"
  for(co in cohorts) stopifnot(identical(as.logical(m[[paste0(co,"_in_set")]]),as.logical(raw[[paste0(co,suffix)]])))
}
width_pt <- 548
row_step <- 15.4
cell <- 11.8
centers <- list(c(301,324,347,370),c(434,457,480,503))
render_panel <- function(panel, ext, filename) {
  ids <- row_order[[panel]]; n <- length(ids); height_pt <- 122+n*row_step+30
  if(ext=="png") png(filename,width=width_pt/72,height=height_pt/72,units="in",res=600,type="cairo",bg="white")
  else cairo_pdf(filename,width=width_pt/72,height=height_pt/72,family="Arial",onefile=TRUE,fallback_resolution=600)
  on.exit(dev.off())
  par(mar=c(0,0,0,0),oma=c(0,0,0,0),family="Arial",xaxs="i",yaxs="i",cex=1)
  plot.new();plot.window(xlim=c(0,width_pt),ylim=c(0,height_pt),asp=1)
  text(13,height_pt-18,panel,adj=c(0,.5),font=2,cex=15/12)
  # Header title is category-colored; every cohort is explicitly labelled.
  for(g in 1:2) {
    group <- group_pairs[[panel]][g];xs<-centers[[g]]
    text(mean(xs),height_pt-18,titles[[group]],font=2,cex=10.2/12,col=colors[[group]])
    for(j in 1:4) text(xs[j],height_pt-103,cohorts[j],srt=90,adj=c(0,.5),cex=8/12)
  }
  for(i in seq_along(ids)) {
    y <- height_pt-122-(i-1)*row_step
    label <- gsub("_"," ",ids[i],fixed=TRUE)
    font <- if(grepl("^GGB",label)) 1 else 3
    # Fixed compact typography and exact unabridged source identifiers.
    stopifnot(strwidth(label,units="user",cex=9.1/12,font=font)<263)
    text(275,y,label,adj=c(1,.5),cex=9.1/12,font=font)
    for(g in 1:2) {
      group<-group_pairs[[panel]][g];flags<-make_flags(ids[i],group)
      for(j in 1:4) {
        x<-centers[[g]][j]
        rect(x-cell/2,y-cell/2,x+cell/2,y+cell/2,
             col=if(flags[j]) colors[[group]] else "white",border="#303030",lwd=.7)
      }
    }
  }
  text(275,14,"Colored: criterion met     White: criterion not met",adj=c(0.5,.5),cex=8/12,col="#444444")
}
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
for(panel in c("a","b")) {
  stem<-paste0("Fig3",panel,"_",if(panel=="a") "Production_CRC" else "Degradation_Healthy","_cohort_membership")
  for(ext in c("png","pdf")) render_panel(panel,ext,file.path(output_dir,paste0(stem,".",ext)))
  message("Fig3",panel,": ",length(row_order[[panel]])," species; ",length(row_order[[panel]])*8," cells")
}
message("Verified frozen input membership; no statistical analysis rerun. R ",getRversion())
