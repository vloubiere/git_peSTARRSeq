setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Selection (Do not run again!)
#-----------------------------------------------#
# dat <- readRDS("db/FC_tables_DESeq2/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_DESeq2_final_oe.rds")
# dat <- merge(unique(dat[, .(ID= L, indL)]),
#              unique(dat[, .(ID= R, indR)]))
# # Add TWIST data
# lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
# twist <- fread("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/BA_300bp_TWIST_STARRSeq.txt")
# lib[twist, dev_log2FC_TWIST:= i.dev_log2FoldChange, on= "BA_ID==ID"]
# dat[lib, dev_log2FC_TWIST:= i.dev_log2FC_TWIST, on= "ID==ID_vl"]
# dat <- dat[!is.na(dev_log2FC_TWIST)]
# 
# set.seed(1)
# weak_candidates <- dat[between(indL, 1, 1.5) & abs(indL-indR)<0.5 & dev_log2FC_TWIST>2][sample(.N, 3)]
# set.seed(1)
# saturating_candidates <- dat[between(indL, 6, Inf) & abs(indL-indR)<0.5][sample(.N, 2)]
# sel <- rbindlist(list(weak_candidates= weak_candidates,
#                       saturating_candidates= saturating_candidates), 
#                  idcol = "group")
# fwrite(sel,
#        "db/library_design/enhancers_upstream_CP/selected_enhancers.txt")

#----------------------#
# Design primers
#----------------------#
dat <- fread("db/library_design/enhancers_upstream_CP/selected_enhancers.txt")
dat[as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds")), seq:= enh_sequence, on= "ID==ID_vl"]

# fw primers
i <- 20
fw <- stringr::str_sub(dat$seq, 1 , i)
tm <- rep(0, length(fw))
while(any(!grepl("C$|G$", fw)) | any(tm<58))
{
  idx <- which(!grepl("C$|G$", fw) | any(tm<58))
  fw[idx] <- stringr::str_sub(dat$seq[idx], 1 , i)
  tm[idx] <- sapply(fw[idx], function(x) vl_oligo_Tm(x)[1])
  i <- i+1
}

# rev primers
i <- 20
rev <- sapply(stringr::str_sub(dat$seq, 249-i+1 , 249), vl_revComp)
tm <- rep(0, length(rev))
while(any(!grepl("C$|G$", rev)) | any(tm<58))
{
  idx <- which(!grepl("C$|G$", rev) | any(tm<58))
  rev[idx] <- sapply(stringr::str_sub(dat$seq[idx], 249-i+1 , 249), vl_revComp)
  tm[idx] <- sapply(rev[idx], function(x) vl_oligo_Tm(x)[1])
  i <- i+1
}
