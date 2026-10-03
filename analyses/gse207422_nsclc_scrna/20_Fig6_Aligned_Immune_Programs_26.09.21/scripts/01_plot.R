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
draw <- function(){
 grid.newpage();grid.rect(gp=gpar(fill="white",col=NA))
 txt("Immune program",.018,.928,13,fontface="bold")
 txt("Whole-tumour pseudobulk",.272,.936,16,fontface="bold")
 txt("12 patients  |  High 7 / Low 5",.272,.88,10.5,col=muted)
 txt("Conventional DCs",.692,.936,16,fontface="bold")
 txt("8 patients  |  High 5 / Low 3",.692,.88,10.5,col=muted)
 txt(expression(bold("CD8")^"+"~bold("T cells")),.692,.655,16,fontface="bold")
 txt("12 patients  |  High 7 / Low 5",.692,.602,10.5,col=muted)
 # Subtle row separation replaces arrows and braces.
 seg(.018,.985,.685)
 for(i in seq_along(ids)){
  txt(labels[[i]],.018,ys[i],12)
  for(side in c("Whole_tumour",if(i==1)"cDC"else"CD8")){
   r <- d[d$pathway==ids[i]&d$compartment==side,];stopifnot(nrow(r)==1)
   left <- if(side=="Whole_tumour") .272 else .692
   width <- .211; sx <- left+width+.014
   w <- width*r$NES/3.5
   grid.rect(left+w/2,ys[i],width=w,height=.043,gp=gpar(fill=if(r$BH_q<.05)red else "white",col=red,lwd=.9))
   txt(sprintf("%.2f",r$NES),sx,ys[i]+.017,11.5,fontface="bold")
   txt(qexpr(r$BH_q),sx,ys[i]-.025,9.5,col=muted)
  }
 }
 for(left in c(.272,.692)){
  width <- .211
  seg(left,left+width,.153,col=muted,lwd=.7)
  for(v in seq(0,3.5,.5)){
   x<-left+width*v/3.5
   grid.lines(c(x,x),c(.153,.145),gp=gpar(col=muted,lwd=.6))
   if(v%%1==0)txt(as.character(v),x,.125,9.5,col=muted,just="centre")
  }
  txt("NES (Degradation High vs Low)",left+width/2,.086,10,col=ink,just="centre")
 }
 grid.rect(.023,.035,width=.010,height=.015,gp=gpar(fill=red,col=red))
 txt(expression(paste("BH ",italic(q)," < 0.05")),.034,.035,9,col=muted)
 grid.rect(.153,.035,width=.010,height=.015,gp=gpar(fill="white",col=red))
 txt(expression(paste("BH ", italic(q) >= 0.05)),.164,.035,9,col=muted)
 txt("Same whole-tumour-defined patient groups; original full-family BH q values retained.",.315,.035,9,col=muted)
}
ragg::agg_png(file.path(out,"figures/A.png"),width=14.5,height=5.0,units="in",res=600,background="white");draw();invisible(dev.off())
grDevices::quartz(file=file.path(out,"figures/A.pdf"),type="pdf",width=14.5,height=5.0,family="Arial",bg="white");draw();invisible(dev.off())
capture.output(sessionInfo(),file=file.path(out,"qa/sessionInfo.txt"))
