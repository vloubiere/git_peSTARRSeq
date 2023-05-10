dat <- data.table(file= list.files("db/umi_counts/", full.names = T))
dat[, c("library", "cdition", "rep"):= tstrsplit(gsub(".txt", "", basename(file)), "_", 1:3)]
dat[, c("N_pairs", "total_counts", "umi_counts"):={
  fread(file)[, .(.N, sum(total_counts), sum(umi_counts))]
}, file]
dat[, Cc:= sapply(cdition, switch, "input"= "grey", "screen"= "tomato")]

pdf("pdf/paper_v1")
vl_par(mfrow= c(2,2), oma= c(0,0,2,0))
dat[, {
  barplot(N_pairs)
  text(grconvertX(.5, "ndc", "user"),
       grconvertY(1, "ndc", "user"),
       library, 
       xpd= NA,
       pos= 1,
       cex= 2)
  bar <- barplot(total_counts,
                 col= adjustcolor(Cc, 0.5))
  vl_tilt_xaxis(bar, 
                labels = paste(cdition, rep))
  barplot(umi_counts, 
          col= NA,
          density= 20,
          add= T,
          yaxt= "n")
}, .(library)]
dev.off()