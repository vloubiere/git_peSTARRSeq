setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

peaks <- data.table(file= list.files("db/peaks/", full.names = T))
peaks[, cdition := gsub("_peaks.txt", "", basename(file))]
peaks <- peaks[, fread(file, colClasses = c("character", rep("numeric", 5))), (peaks)]
# Check that K4me3 peaks correspond to annotated promoters
tss <- as.data.table(import("/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf"))
peaks[, check:= ifelse(cdition!="H3K4me3", T, tss[.SD, .N>0, .EACHI, on= c("seqnames", "start<end", "start>start")]$V1)]
peaks <- peaks[(check), !"check"]
# Compute center
peaks[, center:= rowMeans(.SD), .SDcols= c("start", "end")]


# PLOT object
pl <- list(prom_enh200_distance= peaks[cdition=="STARR_DSCP_200"][peaks[cdition=="H3K4me3"], min(abs(center-i.center)), .EACHI, on= "seqnames"]$V1,
           prom_enh600_distance= peaks[cdition=="STARR_DSCP_600"][peaks[cdition=="H3K4me3"], min(abs(center-i.center)), .EACHI, on= "seqnames"]$V1,
           enh200_enh200_distance= peaks[cdition=="STARR_DSCP_200"][peaks[cdition=="STARR_DSCP_200"], min(abs(center-i.center)[abs(center-i.center)>0]), .EACHI, on= "seqnames"]$V1,
           enh600_enh600_distance= peaks[cdition=="STARR_DSCP_600"][peaks[cdition=="STARR_DSCP_600"], min(abs(center-i.center)[abs(center-i.center)>0]), .EACHI, on= "seqnames"]$V1)
pl <- rbindlist(lapply(pl, as.data.table), idcol = T)
pl$.id <- factor(pl$.id, 
                 levels = c("prom_enh200_distance", 
                            "prom_enh600_distance", 
                            "enh200_enh200_distance", 
                            "enh600_enh600_distance"))


# PLOT
pdf(2, height = 10)
par(mfrow= c(4, 1))
pl[, {
  .c <- V1[is.finite(V1)]
  hist(.c, 
       breaks = seq(0, max(.c)+100, 100), 
       include.lowest = T, 
       xlim= c(0, 1e4),
       border= "black", 
       col= "black",
       main= .id, 
       las= 1,
       xlab= "genomic distance")
  marks <- c(300, 2000)
  abline(v= marks, col= "red")
  text(x= c(300, 2000), 
       y= grconvertY(1, "npc", "user"),
       col= "red", 
       labels = marks,
       xpd= T, 
       pos= 3)
  print("")
}, .id]
dev.off()

