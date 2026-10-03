#!/usr/bin/env Rscript
arg<-sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))]);out<-normalizePath(file.path(dirname(arg),'..'))
suppressPackageStartupMessages({library(png);library(digest);library(jsonlite)})
a<-fromJSON(file.path(out,'results/qa/transport_manifest.json'));checks<-list()
for(i in seq_len(nrow(a))){x<-readPNG(a$transport[i]);stopifnot(length(dim(x))==3,dim(x)[3]==3);px<-as.raw(round(as.vector(aperm(x,c(3,2,1)))*255));h<-digest(px,algo='sha256',serialize=FALSE);ok<-dim(x)[2]==a$width[i]&&dim(x)[1]==a$height[i]&&h==a$pixel_sha256[i];stopifnot(ok);checks[[i]]<-list(file=basename(a$transport[i]),pass=ok,pixel_sha256=h)}
write_json(checks,file.path(out,'results/qa/R_decoder_checks.json'),pretty=TRUE,auto_unbox=TRUE);cat('Independent libpng decoder passed',length(checks),'files\n')
