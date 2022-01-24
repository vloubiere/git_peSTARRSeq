setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)

# Import
dat <- fread("Rdata/metadata_processed.txt")
dat <- dat[vllib=="vllib021"
           & file.exists(pairs_counts) 
           & cdition=="screen"]
dat[, cdition:= ifelse(grepl("UNSPLICED", comment), "UNSPLICED", "SPLICED")]
dat <- dat[, fread(pairs_counts), .(vllib, pairs_counts, cdition, DESeq2_pseudo_rep)]
# Aggregate different sequencing runs
dat <- dat[, .(umi_counts= sum(umi_counts)), .(vllib, cdition, L, R)]
# norm
dat[, norm_counts:= log2((umi_counts+1)/sum(umi_counts)), .(vllib, cdition)]
dat <- dcast(dat, 
             vllib+L+R~cdition,
             value.var = "norm_counts")

pdf("pdf/alignment/PCC_replicates_w_wo_intron_SCREEN.pdf", width = 8.1)
par(mar= c(7,4,2,2))
layout(matrix(1:2, ncol= 2), widths = c(1,0.25))
dat[, {
  smoothScatter(UNSPLICED, 
                SPLICED,
                main= vllib,
                las=1,
                xlab= "'nascent' UNSPLICED (log2 lib. norm. counts)",
                ylab= "'mature'SPLICED (log2 lib. norm. counts)")
  abline(0, 1)
  legend("topleft", 
         paste0("PCC= ", round(cor.test(UNSPLICED, SPLICED)$estimate, 2)),
         bty= "n")
  print(title)
}, vllib]
boxplot(dat[, .(UNSPLICED, SPLICED)], 
        notch= T, 
        las= 2,
        ylab= "lib. norm. UMI counts",
        lty=1,
        staplewex= NA,
        pch= 16,
        outline= F)
dev.off()



