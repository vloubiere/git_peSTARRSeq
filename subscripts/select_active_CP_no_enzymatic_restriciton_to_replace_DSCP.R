setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(BSgenome.Dmelanogaster.UCSC.dm3)

#-------------------#
# Import 72k STAP-Seq library
#-------------------#
lib <- readRDS("/groups/stark/lorbeer/vincent/STAP_oligo.RDS")
lib[, sequence:= as.character(getSeq(BSgenome.Dmelanogaster.UCSC.dm3, GRanges(seqnames, IRanges(start, end))))]

#-------------------#
# Import constructs
#-------------------#
# Import my constructs with varying CPs
cons <- fread("/groups/stark/vloubiere/exp_data/vl_constructs_sequences.txt")
cons <- cons[Class=="CP" & grepl("^actCP_", name), !c("Class", "detail", "sequence")]
cons[, oligo_id:= paste0(c(seqnames, start, end, strand), collapse= "_"), (cons)]
cons[, owner:= "vl"]
# Add Franzi's ones
cons <- rbind(cons, 
              data.table(oligo_id= c("chr2R_18090097_18090230_-",
                                     "chr2R_2520526_2520659_-",
                                     "chr3L_297736_297869_+",
                                     "chr2R_8392163_8392296_-",
                                     "chr2R_12579143_12579276_+",
                                     "chr3R_16783401_16783534_+",
                                     "chr2L_3478367_3478500_+",
                                     "chr3L_217307_217440_+",
                                     "chr2R_15057666_15057799_+"),
                         owner= "fl"), 
              fill= T)
suppressWarnings(cons[owner=="fl", c("seqnames", "start", "end", "strand"):= tstrsplit(oligo_id, "_"), oligo_id])
cons[owner=="fl", name:= paste0("actCP_fl", seq(.N))]

#-------------------#
# Color data based on overlaps
#-------------------#
dat <- copy(lib)
dat[, x:= log2(no_enh_Rep1+1)]
dat[, y:= log2(zfh1_enh_avg+1)]
dat[, z:= log2(ssp3_enh_Rep1+1)]
# Present in Franzis or my  selection
dat[cons, c("Cc", "name"):= .(ifelse(i.owner=="fl", "red", "green"),
                              i.name), on= "oligo_id"]
# Overlaps Rps12 promoter
RPS12_tss <- data.table(seqnames= "chr3L", 
                        start= 13016456, 
                        end= 13016456, strand= "+")
dat[RPS12_tss, c("Cc", "name"):= .("gold", "RpS12"), on= c("seqnames", "start<end", "end>start", "strand")]
#Check for the presence of restriction sites incompatible with STARR-Seq
dat[!is.na(Cc) & grepl("ACCGGT|GTCGAC", sequence), Cc:= "blue"] # AgeI|SalI

#-------------------#
# Make region design object
#-------------------#
design <- data.table(xleft= c(1.6, 2.1, 2.1, 5.1)-0.5,
                     ybottom= c(8.5, 1.9, 5.5, 11.25)-0.5,
                     xright= c(2.2, 2.7, 2.7, 5.85)+0.4,
                     ytop= c(9.2, 2.7, 6.3, 12)+0.5)
design[order(-ytop, xright), name:= paste0("Design_act", seq(4))]
dat[design, design:= ifelse(is.na(Cc), i.name, as.numeric(NA)), on= c("x>=xleft", "x<=xright", "y>=ybottom", "y<=ytop")]
dat[!is.na(design), F_primer:= substr(sequence, 1, 22)]
dat[!is.na(design), c("F_Tm", "F_GC"):= vl_oligo_Tm(F_primer), F_primer]
dat[!is.na(design), F_Clamp:= grepl("C$|G$", F_primer), F_primer]
dat[!is.na(design), R_primer:= substr(vl_revComp(sequence), 1, 22), sequence]
dat[!is.na(design), c("R_Tm", "R_GC"):= vl_oligo_Tm(R_primer), R_primer]
dat[!is.na(design), R_Clamp:= grepl("C$|G$", R_primer), R_primer]
dat[!is.na(design) & 
        F_Tm>56 & R_Tm>56 & abs(F_Tm-R_Tm)<1 &
        between(F_GC, 35, 65) & between(R_GC, 35, 65) &
        F_Clamp & R_Clamp, check:= ifelse(seq(.N)<3, T, NA), design]
setorderv(dat, c("y", "x"), order = c(-1, 1))
dat[(check), c("Cc", "name"):= .("purple", 
                                 paste0("actCP_", seq(.N)+3))]

#-------------------#
# PLOT
#-------------------#
# Order so that colored dots appear first
pdf("pdf/design/CPs_responsiveness.pdf")
par(pty= "s", 
    las= 1)
# Plot developmental induction
smoothScatter(dat$x,
              dat$y,
              colramp = colorRampPalette(c("white", "grey30", "grey20", "grey10")),
              xlab= "No enhancer",
              ylab= "ZFH1 enhancer")
points(dat[!is.na(Cc), x],
       dat[!is.na(Cc), y], 
       col= dat[!is.na(Cc), Cc],
       pch= 19)
text(dat[!is.na(Cc), x],
     dat[!is.na(Cc), y], 
     labels = dat[!is.na(Cc), name], 
     pos= 4,
     cex= 0.2,
     col= dat[!is.na(Cc), Cc])
legend("bottomright", 
       bty= "n",
       pch= 19, 
       col= c("red", "green", "gold", "blue"),
       legend = c("Franzi", "vl", "RpS12", "AgeI/SalI"))
# Design areas
design[, {
    rect(xleft, 
         ybottom,
         xright, 
         ytop, 
         lty= 2)
    text(xright,
         ytop,
         name,
         pos= 4)
}, (design)]
# Plot housekeeping induction
smoothScatter(dat$x,
              dat$z,
              colramp = colorRampPalette(c("white", "grey30", "grey20", "grey10")),
              xlab= "No enhancer",
              ylab= "SSP3 enhancer")
points(dat[!is.na(Cc), x],
       dat[!is.na(Cc), z], 
       col= dat[!is.na(Cc), Cc],
       pch= 19)
text(dat[!is.na(Cc), x],
     dat[!is.na(Cc), z], 
     labels = dat[!is.na(Cc), name], 
     pos= 4,
     cex= 0.2,
     col= dat[!is.na(Cc), Cc])
legend("bottomright", 
       bty= "n",
       pch= 19, 
       col= c("red", "green", "gold", "blue"),
       legend = c("Franzi", "vl", "RpS12", "AgeI/SalI"))
dev.off()

#------------------------------#
# SAVE selected primers
#------------------------------#
fwrite(dat[Cc=="purple", .(oligo_id, name, seqnames, start, end, strand, F_primer, R_primer, CP_sequence= sequence)], 
       "db/library_design/alternative_CPs/selected_CPs_primers_#2.txt", 
       sep= "\t")
