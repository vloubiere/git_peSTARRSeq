setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)

dat <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
dat <- as.data.table(dat)
cols <- colnames(dat)
dat[, (cols):= lapply(.SD, function(x) ifelse(x=="NA", NA, x)), .SDcols= cols]

dat <- dat[, .(file= list.files("db/merged_counts/", paste0(DESeq2_group, "_", cdition, "_", CP, "_rep", DESeq2_pseudo_rep, ".*merged.txt$"), full.names = T)), 
           .(DESeq2_group, cdition, DESeq2_pseudo_rep)]
dat <- dat[, fread(file), (dat)]
dat <- dat[type!="switched"]

PCC <- dcast(dat, DESeq2_group+cdition+L+R~DESeq2_pseudo_rep,
             value.var= "umi_counts")
names(PCC)[(length(PCC)-1):length(PCC)] <- paste0("rep", names(PCC)[(length(PCC)-1):length(PCC)])

pdf("pdf/alignment/PCC_replicates.pdf")
PCC[L!=R & !is.na(rep1) & !is.na(rep2), {
  x <- log2(rep1)
  y <- log2(rep2)
  smoothScatter(x, 
                y,
                main= paste(cdition, DESeq2_group),
                las=1)
  legend("topleft", 
         paste0("PCC= ", round(cor.test(x, y)$estimate, 2)),
         bty= "n")
  print("")
}, .(DESeq2_group, cdition)]
dev.off()
