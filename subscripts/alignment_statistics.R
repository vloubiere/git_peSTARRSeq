setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

lib <- as.data.table(readxl::read_xlsx("../../exp_data/vl_libraries.xlsx"))
dat <- data.table(file= list.files("db/umi_counts/", full.names = T))
dat[, c("library", "cdition", "rep"):= tstrsplit(gsub(".txt", "", basename(file)), "_", 1:3)]
dat[lib, total_pairs:= i.Complexity, on= "library==lib_ID"]
dat[, c("N_pairs", "total_counts", "collapsed_counts", "umi_counts"):= {
  fread(file)[, .(.N, 
                  sum(total_counts), 
                  sum(umi_counts), 
                  .(c(rep(0, total_pairs-.N), umi_counts)))]
}, .(file, total_pairs)]
dat[, min_reads:= formatC(min(collapsed_counts), format = "e", digits =  1), .(library, cdition)]
dat[, Cc:= sapply(cdition, switch, "input"= "grey", "screen"= "tomato")]

pdf("pdf/draft/alignment_statistics.pdf",
    width = 14,
    height = 4)
par(mfrow= c(1,4), 
    oma= c(0,0,2,10), 
    mgp= c(3.5,0.5,0),
    mar= c(5,5,2,2),
    cex= 1,
    las=1,
    tcl= -0.2)
dat[, {
  # Barplot number of pairs
  bar <- barplot(N_pairs, 
                 ylab= "Number of pairs",
                 ylim= c(0, total_pairs))
  abline(h= total_pairs, lty= 2)
  vl_tilt_xaxis(bar, 
                labels = paste(cdition, rep))
  # violin plot reads per pair
  vl_boxplot(umi_counts,
             names= paste(cdition, rep),
             ylab= "Number of reads", 
             tilt.names= T)
  # Density
  plot(NA, 
       xlim= c(0,15), 
       ylim= c(0,1),
       ylab= "Density",
       xlab= "Supporting reads+1 (log2)")
  lapply(seq(umi_counts), function(i) {
    lines(density(log2(umi_counts[[i]]+1), 
                  from= 0, 
                  to= 15), 
          col= Cc[i],
          lwd= 2)
  })
  unique(.SD[, .(Cc, cdition)])[, {
    legend("topleft",
           lty= c(1, 1),
           legend= cdition,
           col= Cc,
           bty= "n")
  }]
  # Library name title
  text(grconvertX(.5, "ndc", "user"),
       grconvertY(1, "ndc", "user"),
       library[1], 
       xpd= NA,
       pos= 1,
       cex= 2)
  # Barplot number of reads (before and after collapsing)
  bar <- barplot(total_counts,
                 col= adjustcolor(Cc, 0.5),
                 ylab= "Number of reads")
  vl_tilt_xaxis(bar, 
                labels = paste(cdition, rep))
  barplot(collapsed_counts, 
          col= NA,
          density= 20,
          add= T,
          yaxt= "n")
  unique(.SD[, .(cdition, min_reads, Cc)])[, {
    legend(par("usr")[2],
           par("usr")[4],
           c("Total reads", "Collapsed reads", paste0(cdition, " (min= ", min_reads, ")")),
           fill= c("white", NA, adjustcolor(Cc, 0.5)),
           density= c(NA, 20, NA, NA),
           bty= "n",
           xpd= NA)
  }]
  .SD
}, .(library, total_pairs)]
dev.off()