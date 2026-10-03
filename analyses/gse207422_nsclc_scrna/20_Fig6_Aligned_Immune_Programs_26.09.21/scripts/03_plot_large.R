args <- commandArgs(FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)]))
out <- dirname(dirname(script))
suppressPackageStartupMessages(library(grid))
d <- read.csv(file.path(out,"tables/plot_data.csv"),check.names=FALSE)
stopifnot(nrow(d)==8,all(d$NES>0),sum(d$BH_q>=.05)==1)
ids <- c("REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION","REACTOME_TCR_SIGNALING","REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY","HALLMARK_INTERFERON_GAMMA_RESPONSE")
labels <- list("MHC-I antigen\ncross-presentation","TCR signaling","TNFR2-related\nnoncanonical NF-κB",expression(paste("IFN-",gamma," response")))
ys <- c(.76,.51,.375,.24)
ink <- "#233741"; muted <- "#647780"; red <- "#BC554C"; line <- "#DCE3E6"
txt <- function(label,x,y,size=11,col=ink,fontface="plain",just="left")grid.text(label,x,y,just=just,gp=gpar(fontfamily="Arial",fontsize=size,col=col,fontface=fontface))
seg <- function(x0,x1,y,col=line,lwd=.6)grid.lines(c(x0,x1),c(y,y),gp=gpar(col=col,lwd=lwd))
qexpr <- function(q){if(q>=.001) bquote(italic(q)==.(formatC(q,format="f",digits=3))) else {e<-floor(log10(q));m<-formatC(q/10^e,format="f",digits=1);bquote(italic(q)==.(m)%*%10^.(e))}}
draw_combined <- function(){
 grid.newpage();grid.rect(gp=gpar(fill="white",col=NA))
 txt("Immune program",.018,.928,17,fontface="bold")
 txt("Whole-tumour pseudobulk",.272,.936,21,fontface="bold")
 txt("12 patients  |  High 7 / Low 5",.272,.88,14,col=muted)
 txt("Conventional DCs",.692,.936,21,fontface="bold")
 txt("8 patients  |  High 5 / Low 3",.692,.88,14,col=muted)
 txt(expression(bold("CD8")^"+"~bold("T cells")),.692,.655,21,fontface="bold")
 txt("12 patients  |  High 7 / Low 5",.692,.602,14,col=muted)
 # Subtle row separation replaces arrows and braces.
 seg(.018,.985,.685)
 for(i in seq_along(ids)){
  txt(labels[[i]],.018,ys[i],17)
  for(side in c("Whole_tumour",if(i==1)"cDC"else"CD8")){
   r <- d[d$pathway==ids[i]&d$compartment==side,];stopifnot(nrow(r)==1)
   left <- if(side=="Whole_tumour") .272 else .692
   width <- .211; sx <- left+width+.014
   w <- width*r$NES/3.5
   grid.rect(left+w/2,ys[i],width=w,height=.043,gp=gpar(fill=if(r$BH_q<.05)red else "white",col=red,lwd=.9))
   txt(sprintf("%.2f",r$NES),sx,ys[i]+.017,16,fontface="bold")
   txt(qexpr(r$BH_q),sx,ys[i]-.025,13.5,col=muted)
  }
 }
 for(left in c(.272,.692)){
  width <- .211
  seg(left,left+width,.153,col=muted,lwd=.7)
  for(v in seq(0,3.5,.5)){
   x<-left+width*v/3.5
   grid.lines(c(x,x),c(.153,.145),gp=gpar(col=muted,lwd=.6))
   if(v%%1==0)txt(as.character(v),x,.125,13.5,col=muted,just="centre")
  }
  txt("NES (Degradation High vs Low)",left+width/2,.086,13,col=ink,just="centre")
 }
 grid.rect(.023,.035,width=.010,height=.015,gp=gpar(fill=red,col=red))
 txt(expression(paste("BH ",italic(q)," < 0.05")),.034,.035,11,col=muted)
 grid.rect(.153,.035,width=.010,height=.015,gp=gpar(fill="white",col=red))
 txt(expression(paste("BH ", italic(q) >= 0.05)),.164,.035,11,col=muted)
 txt("Whole-tumour-defined groups; original full-family BH q values.",.315,.035,11,col=muted)
}

# Standalone panels use the same physical font and NES scales.
draw_single <- function(comp, height){
 grid.newpage();grid.rect(gp=gpar(fill="white",col=NA))
 pushViewport(viewport(xscale=c(0,8.8),yscale=c(0,height)))
 t <- function(label,x,y,size=17,col=ink,fontface="plain",just="left")
  grid.text(label,unit(x,"native"),unit(y,"native"),just=just,gp=gpar(fontfamily="Arial",fontsize=size,col=col,fontface=fontface))
 ln <- function(x0,x1,y,col=line,lwd=.7) grid.lines(unit(c(x0,x1),"native"),unit(c(y,y),"native"),gp=gpar(col=col,lwd=lwd))
 title <- switch(comp,Whole_tumour="Whole-tumour pseudobulk",cDC="Conventional DCs",CD8=expression(bold("CD8")^"+"~bold("T cells")))
 sub <- if(comp=="cDC") "8 patients  |  High 5 / Low 3" else "12 patients  |  High 7 / Low 5"
 t(title,.15,height-.28,22,fontface="bold");t(sub,.15,height-.66,14,col=muted)
 sel <- if(comp=="Whole_tumour") 1:4 else if(comp=="cDC") 1 else 2:4
 for(k in seq_along(sel)){
  i<-sel[k];y<-height-1.2-(k-1)*.65
  r<-d[d$pathway==ids[i]&d$compartment==comp,];stopifnot(nrow(r)==1)
  t(labels[[i]],.15,y,17)
  w<-3.25*r$NES/3.5
  grid.rect(unit(3.45+w/2,"native"),unit(y,"native"),width=unit(w,"native"),height=unit(.23,"native"),gp=gpar(fill=if(r$BH_q<.05)red else "white",col=red,lwd=1.2))
  t(sprintf("%.2f",r$NES),6.92,y+.12,17,fontface="bold")
  t(qexpr(r$BH_q),6.92,y-.16,14,col=muted)
 }
 ln(3.45,6.7,.96,col=muted)
 for(v in seq(0,3.5,.5)){
  x<-3.45+3.25*v/3.5
  grid.lines(unit(c(x,x),"native"),unit(c(.96,.90),"native"),gp=gpar(col=muted,lwd=.7))
  if(v%%1==0)t(as.character(v),x,.73,13,col=muted,just="centre")
 }
 t("NES (Degradation High vs Low)",5.075,.44,14,just="centre")
 t(expression(paste("Filled: BH ",italic(q)," < 0.05   |   Outline: BH ",italic(q)>=0.05)),.15,.13,11,col=muted)
 popViewport()
}
export <- function(name,width,height,drawfun){
 ragg::agg_png(file.path(out,paste0("figures/",name,".png")),width=width,height=height,units="in",res=600,background="white");drawfun();invisible(dev.off())
 grDevices::quartz(file=file.path(out,paste0("figures/",name,".pdf")),type="pdf",width=width,height=height,family="Arial",bg="white");drawfun();invisible(dev.off())
}
export("A2",14.5,6.4,draw_combined)
export("W",8.8,4.6,function()draw_single("Whole_tumour",4.6))
export("D",8.8,2.65,function()draw_single("cDC",2.65))
export("T",8.8,3.95,function()draw_single("CD8",3.95))
