suppressPackageStartupMessages({library(png);library(digest);library(data.table)})
a<-commandArgs(TRUE)[1];fs<-list.files(file.path(a,'results/figures/PPT_insert'),pattern='[.]png$',full.names=TRUE)
r<-rbindlist(lapply(fs,function(f){x<-readPNG(f);data.table(file=basename(f),height=dim(x)[1],width=dim(x)[2],channels=dim(x)[3],RGB_pixel_sha256=digest(as.raw(round(aperm(x,c(3,2,1))*255)),algo='sha256',serialize=FALSE))}))
fwrite(r,file.path(a,'results/qa/R_libpng_pixel_checks.csv'));print(r)
