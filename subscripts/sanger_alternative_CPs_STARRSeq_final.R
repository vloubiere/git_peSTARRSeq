setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
dir <- "db/sanger_sequencing/alternative_CPs_STARRSeq_final/"
dir.create(dir, showWarnings = F)
require(vlfunctions)

#---------------------------------------------------#
# IMPORT FILES AND COPY THEM LOCALLY
#---------------------------------------------------#
# source("../../exp_data/update_files.R")
dat <- as.data.table(read_xlsx(path = "/groups/stark/vloubiere/exp_data/vl_sanger_sequencing.xlsx"))
dat <- dat[Date==211122]
dat[, molBio:= list.files("/groups/MolBioService/CollaborationData/SEQUENCING RESULTS-new/2021/LOUB/2021-11-22_PCW-IMP-232_0160/",
                          paste0("_", seqID, "_"),
                          full.names = T), seqID]
dat[, file:= list.files(dir,
                        paste0("_", seqID, "_"),
                        full.names = T), seqID]
dat[is.na(file) & file.exists(molBio), file:= paste0(dir, basename(molBio))]
dat[!file.exists(file), file.copy(molBio, file), .(molBio, file)]

#---------------------------------------------------#
# IMPORT sequences
#---------------------------------------------------#
pl <- as.data.table(readxl::read_xlsx("../../exp_data/vl_plasmids.xlsx"), 
                    key= "ID")
constructs <- fread("../../exp_data/vl_constructs_sequences.txt", 
                    key= "name")
dat[, refseq:= pl[.BY, Sequence, on= "ID==current"], .(current= gsub("(.*)_.*", "\\1", name))]
# dat[, refseq:= pl[ID=="pGL3_IF_DSCPII_mhcI_GFP_AgeI-SalI", Sequence]]
dat[, refseq:= paste0(substr(refseq, 4000, nchar(refseq)), substr(refseq, 1, 3999))]

#---------------------------------------------------#
# IMPORT sequences
#-------------------------------------------------------#
pdf("pdf/design/sanger/sanger_alternative_CPs_FINAL_20211122.pdf")
dat[, {
  vl_sanger_align(refseq, 
                  abfiles = file, 
                  revcomp = as.logical(revComp), 
                  feat_sequences = constructs["actCP_3", sequence])
  mtext(name, line = 3)
}, .(name, refseq)]
dev.off()
