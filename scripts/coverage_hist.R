id_out <- args[1]
depth_hist <- args[2]

sequencing_depth <- read.delim(depth_hist,header=F)

colnames(sequencing_depth) <- c("rname","position","depth")
sequencing_depth$logdepth <- log10(sequencing_depth$depth)

png(paste0(id_out,"_log_depth.png"))
plot(sequencing_depth$position,sequencing_depth$logdepth,type="l",lty=1,xlab="Position",
     ylab="log10(depth)",main=paste0(id_out, "sequencing depth"), ylim=c(0,6))
abline(h=log10(20),col="red")
dev.off()