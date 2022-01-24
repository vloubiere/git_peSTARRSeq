setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)

# Import data
dat <- fread("Rdata/metadata_processed.txt")
dat <- dat[vllib %in% c("vllib017", "vllib019", "vllib021") 
           & file.exists(pairs_counts) 
           & cdition=="input"]
dat[, cdition:= ifelse(grepl("no_spacer_digestion", comment), "intron", "no_intron")]
dat <- dat[, fread(pairs_counts), .(vllib, pairs_counts, cdition, DESeq2_pseudo_rep)]
# Merge all sequencing runs
dat <- dat[, .(umi_counts= sum(umi_counts)), .(vllib, cdition, L, R)]
# Normalize
dat[, norm_counts:= log2((umi_counts+1)/sum(umi_counts)), .(vllib, cdition)]
# cast
dat <- dcast(dat, 
             vllib+L+R~cdition,
             value.var = "norm_counts")

pdf("pdf/alignment/PCC_replicates_w_wo_intron_INPUT.pdf")
dat[, {
  smoothScatter(intron,
                no_intron,
                main= vllib,
                las=1,
                ylab= "intron not removed (log2 lib. norm. counts)",
                xlab= "intron removed (log2 lib. norm. counts)")
  abline(0, 1)
  legend("topleft", 
         paste0("PCC= ", round(cor.test(intron, no_intron)$estimate, 2)),
         bty= "n")
  print(title)
}, vllib]
dev.off()

