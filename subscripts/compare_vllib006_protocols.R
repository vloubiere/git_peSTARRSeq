setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)
require(readxl)

# Import data
if(!exists("dat"))
{
  dat <- list(old= "db/umi_counts/vllib006_rev-pe-STARR-Seq_DSCP_T8_revSPA1_2000_775_T_input_vllib011_0.1_200ng_rep1__HK2HFBGXH_1_20210205B_20210206.txt",
              new_r1= "db/merged_counts/vllib006_rev-pe-STARR-Seq_DSCP_T8_revSPA1_2000_input_rep1_merged.txt",
              new_r2= "db/merged_counts//vllib006_rev-pe-STARR-Seq_DSCP_T8_revSPA1_2000_input_rep2_merged.txt")
  dat <- lapply(dat, fread)
  dat[[1]] <- dat[[1]][, .(umi_counts= .N), .(L, R)]
  dat <- rbindlist(dat, fill= T, idcol = T)
  # Only keep pairs that are not switched
  dat <- dat[unique(dat[type!="switched", .(L, R)]), on= c("L", "R")]
  # Add merge new_r1 and new_r2
  dat <- rbind(dat, 
               dat[grepl("new", .id), .(umi_counts= sum(umi_counts), .id= "new_merged"), .(L, R)], 
               fill= T)
}
w_spikein <- dcast(dat, L+R~.id, value.var = "umi_counts")
wo_spikein <- dcast(dat[unique(dat[type=="pair", .(L, R)]), on= c("L", "R")], L+R~.id, value.var = "umi_counts")

pdf("pdf/alignment/PCC_40x_20x_2ndStrand_revSTARRSeq_inputs.pdf")
chart.Correlation(log2(as.matrix(w_spikein[, new_merged:old])), 
                  main= "With Spike-in")
chart.Correlation(log2(as.matrix(wo_spikein[, new_merged:old])), 
                  main= "Without Spike-in")
dev.off()
