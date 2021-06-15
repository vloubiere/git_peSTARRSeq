dat <- readRDS("/groups/stark/lorbeer/vincent/STAP_oligo.RDS")

FL_clones <- data.table(coor= c("chr2R_18090097_18090230_-",
                                "chr2R_2520526_2520659_-",
                                "chr3L_297736_297869_+",
                                "chr2R_8392163_8392296_-",
                                "chr2R_12579143_12579276_+",
                                "chr3R_16783401_16783534_+",
                                "chr2L_3478367_3478500_+",
                                "chr3L_217307_217440_+",
                                "chr2R_15057666_15057799_+"))
FL_clones[, c("seqnames", "start", "end", "strand"):= tstrsplit(coor, "_"), coor]

plot(log2(zfh1_enh_avg+1)~log2(no_enh_Rep1+1), 
     dat, 
     col= "grey", 
     pch = 19, 
     xlim= c(-1,15),
     ylim= c(-1,15))
abline(0, 1, col= "red")
points(dat[oligo_id %in% FL_clones$coor, .(log2(no_enh_Rep1+1), log2(zfh1_enh_avg+1))], col= "red")
points(dat[46983, .(log2(no_enh_Rep1+1), log2(zfh1_enh_avg+1))], col= "red")

# identify(log2(dat$no_enh_Rep1+1), log2(dat$zfh1_enh_avg+1))
candidates <- dat[c(46983, 2477)]
seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3, GRanges(candidates[, .(seqnames, start, end)]))

actCP1_F <- substr(seqs[1], 1, 24)
actCP1_R <- as.character(reverseComplement(DNAString(substr(seqs[1], nchar(seqs[1])-21, nchar(seqs[1])))))
actCP2_F <- substr(seqs[2], 1, 22)
actCP2_R <- as.character(reverseComplement(DNAString(substr(seqs[2], nchar(seqs[2])-21, nchar(seqs[2])))))

vl_oligo_Tm(actCP1_F)
vl_oligo_Tm(actCP1_R)
vl_oligo_Tm(actCP2_F)
vl_oligo_Tm(actCP2_R)
