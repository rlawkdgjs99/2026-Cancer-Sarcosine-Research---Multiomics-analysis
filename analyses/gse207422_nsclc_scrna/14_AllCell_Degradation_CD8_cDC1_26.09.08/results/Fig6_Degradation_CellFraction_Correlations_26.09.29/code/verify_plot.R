args <- commandArgs(trailingOnly=TRUE)
o <- args[1]
d <- read.csv(file.path(o,"Patient_data.csv"),check.names=FALSE)
s <- read.csv(file.path(o,"Correlation_statistics.csv"),check.names=FALSE)
keys <- c("CD8","cDC","Epithelial")
# Independent base-R coefficient, BH and exhaustive restricted-permutation validation.
perms <- function(v) {if(length(v)==1) return(matrix(v,nrow=1)); do.call(rbind,lapply(seq_along(v),function(j) cbind(v[j],perms(v[-j]))))}
checks <- list()
for(scope in c("Post12","All15")) {
 dd <- if(scope=="Post12") subset(d,treatment=="Post") else d
 z <- apply(dd[c("degradation_mean_z",paste0(keys,"_percent"))],2,rank)
 cov <- if(scope=="Post12") model.matrix(~histology,dd) else model.matrix(~histology+treatment,dd)
 res <- qr.resid(qr(cov),z)
 for(j in 1:3) {
  a <- s[s$scope==scope & s$key==keys[j],]
  stopifnot(abs(cor(dd$degradation_mean_z,dd[[paste0(keys[j],"_percent")]],method="spearman")-a$rho[a$analysis=="Spearman"])<1e-12)
  rr <- cor(res[,1],res[,j+1]); stopifnot(abs(rr-a$rho[a$analysis=="Partial_Spearman"])<1e-12)
  strata <- if(scope=="Post12") dd$histology else interaction(dd$histology,dd$treatment,drop=TRUE)
  sums <- 0
  for(ix in split(seq_len(nrow(dd)),strata)) {
    pp <- perms(ix); contrib <- as.numeric(matrix(res[pp,1],nrow=nrow(pp)) %*% res[ix,j+1])
    sums <- as.vector(outer(sums,contrib,"+"))
  }
  rperm <- sums/sqrt(sum(res[,1]^2)*sum(res[,j+1]^2))
  exactp <- mean(abs(rperm)>=abs(rr)-1e-12)
  stopifnot(abs(exactp-a$p[a$analysis=="Partial_Spearman"])<1e-12)
  checks[[length(checks)+1]] <- data.frame(scope=scope,key=keys[j],rho_R=cor(z[,1],z[,j+1]),partial_rho_R=rr,partial_exact_p_R=exactp,passed=TRUE)
 }
 for(method in unique(s$analysis)) {
   a<-s[s$scope==scope & s$analysis==method,];stopifnot(max(abs(p.adjust(a$p,"BH")-a$BH_q))<1e-12)
 }
}
write.csv(do.call(rbind,checks),file.path(o,"Independent_R_verification.csv"),row.names=FALSE)
library(ggplot2);library(patchwork)
dir.create(file.path(o,"figures"),showWarnings=FALSE)
titles <- c("CD8+ T cells","Conventional DCs","Epithelial cells")
colors <- c("#078B98","#7542B6","#D45264")
for(scope in c("Post12","All15")) {
 dd<-if(scope=="Post12") subset(d,treatment=="Post") else d
 plots<-list()
 for(j in 1:3) {
  a<-s[s$scope==scope & s$analysis=="Spearman" & s$key==keys[j],]
  dd$percent<-dd[[paste0(keys[j],"_percent")]]
  dd$Histology<-factor(dd$histology,levels=c("Adeno","Squamous"),labels=c("Adenocarcinoma","Squamous"))
  p<-ggplot(dd,aes(degradation_mean_z,percent))+
    geom_point(aes(shape=Histology),size=3.8,color=colors[j],alpha=.95,stroke=.7)+
    scale_shape_manual(values=c(16,17))+
    guides(shape=guide_legend(override.aes=list(colour="#465E69",alpha=1)))+
    scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(.025,.09)))+
    labs(title=titles[j],subtitle=sprintf("%s | n = %d\nSpearman rho = %+.3f | BH q = %.3f",if(scope=="Post12") "Post-treatment" else "All patients",nrow(dd),a$rho,a$BH_q),x="Sarcosine degradation score",y="Fraction of captured cells (%)")+
    theme_classic(base_size=16,base_family="Arial")+
    theme(plot.title=element_text(face="bold",size=22,color="#203744",margin=margin(b=9)),plot.subtitle=element_text(size=15,lineheight=1.2,color="#556B77",margin=margin(b=16)),axis.title=element_text(size=16),axis.title.x=element_text(margin=margin(t=12)),axis.title.y=element_text(margin=margin(r=10)),axis.text=element_text(size=14,color="#203744"),axis.line=element_line(color="#657C87",linewidth=.5),legend.position="bottom",legend.title=element_blank(),legend.text=element_text(size=12),plot.margin=margin(18,18,12,18))
  plots[[j]]<-p
  if(scope=="Post12") {
   ggsave(file.path(o,"figures",paste0(keys[j],".png")),p,width=6.2,height=5.8,dpi=320,device=ragg::agg_png,bg="white")
   ggsave(file.path(o,"figures",paste0(keys[j],".pdf")),p,width=6.2,height=5.8,device=cairo_pdf,bg="white")
  }
 }
 comb<-wrap_plots(plots,nrow=1,guides="collect") & theme(legend.position="bottom")
 ggsave(file.path(o,"figures",paste0(scope,".png")),comb,width=17.5,height=5.5,dpi=320,device=ragg::agg_png,bg="white")
 ggsave(file.path(o,"figures",paste0(scope,".pdf")),comb,width=17.5,height=5.5,device=cairo_pdf,bg="white")
}
writeLines(capture.output(sessionInfo()),file.path(o,"R_sessionInfo.txt"))
cat("Independent R checks passed; PNG/PDF figures saved.\n")
