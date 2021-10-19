setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
dir <- "db/sanger_sequencing/actCP3_DSCP_20210921/"
dir.create(dir, showWarnings = F)
require(vlfunctions)

#---------------------------------------------------#
# IMPORT FILES AND COPY THEM LOCALLY
#---------------------------------------------------#
# source("../../exp_data/update_files.R")
dat <- as.data.table(read_xlsx(path = "/groups/stark/vloubiere/exp_data/vl_sanger_sequencing.xlsx"))
dat <- dat[grepl("pVL159_actCP3", name) & Date==210921 & project=="pe_STARRSeq"]
dat[, molBio:= list.files("/groups/MolBioService/CollaborationData/SEQUENCING RESULTS-new/2021/LOUB/2021-09-21_PCW-IMP-232_0371/",
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
pl <- fread("../../exp_data/vl_plasmids.txt", 
            key= "ID")
constructs <- fread("../../exp_data/vl_constructs_sequences.txt", 
                    key= "name")
dat[, refseq:= pl["pVL159", Sequence]]
dat[, refseq:= paste0(substr(refseq, 4000, nchar(refseq)), substr(refseq, 1, 3999))]

#---------------------------------------------------#
# IMPORT sequences
#-------------------------------------------------------#
pdf("pdf/sanger_sequencing/sanger_actCP3_DSCP_20210921.pdf")
dat[, {
  vl_sanger_align(refseq, 
                  abfiles = file, 
                  revcomp = as.logical(revComp), 
                  feat_sequences = constructs["actCP_3", sequence])
  mtext(name, line = 3)
}, .(name, refseq)]
dev.off()
