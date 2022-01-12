setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)

input <- fread("Rdata/metadata_processed.txt")
input <- input[vllib %in% c("vllib017", "vllib019", "vllib021") 
           & file.exists(pairs_counts) 
           & cdition=="input"]
input[, .(DESeq2_pseudo_rep, comment), vllib]
input[is.na(DESeq2_pseudo_rep), DESeq2_pseudo_rep:= 2]
input <- input[, fread(pairs_counts), .(vllib, pairs_counts, cdition, DESeq2_pseudo_rep)]

screen <- fread("Rdata/metadata_processed.txt")
screen <- screen[vllib=="vllib021"
               & file.exists(pairs_counts) 
               & cdition=="screen"]
screen[, .(DESeq2_pseudo_rep, comment), vllib]
screen[!is.na(DESeq2_pseudo_rep), DESeq2_pseudo_rep:= 1]
screen[is.na(DESeq2_pseudo_rep), DESeq2_pseudo_rep:= 2]
screen <- screen[, fread(pairs_counts), .(vllib, pairs_counts, cdition, DESeq2_pseudo_rep)]

dat <- rbind(screen, input)


pdf("pdf/alignment/PCC_replicates_w_wo_intron.pdf")
dat[, {
  title <- paste(vllib, cdition)
  mat <- dcast(.SD, 
               L+R~DESeq2_pseudo_rep, 
               value.var= "umi_counts",
               fun.aggregate= sum)
  x <- log2(mat[, `1`]+1)
  y <- log2(mat[, `2`]+1)
  smoothScatter(x, 
                y,
                main= title,
                las=1,
                xlab= "intron removed",
                ylab= "intron not removed")
  legend("topleft", 
         paste0("PCC= ", round(cor.test(x, y)$estimate, 2)),
         bty= "n")
  print(title)
}, .(vllib, cdition)]
dev.off()


file.show("pdf/alignment/PCC_replicates_w_wo_intron.pdf")