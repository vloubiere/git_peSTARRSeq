setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
dir.create("db/sanger_sequencing/actCP_DSCP_20210714/", showWarnings = F)
require(data.table)
require(sangerseqR)
require(vlfunctions)
require(stringr)
require(readxl)

#---------------------------------------------------#
# IMPORT FILES 
#---------------------------------------------------#
dat <- as.data.table(read_xlsx(path = "/groups/stark/vloubiere/exp_data/vl_sanger_sequencing.xlsx"))
dat <- dat[grepl("actCP2", name) & Date==210714 & project=="pe_STARRSeq"]
# files <- list.files(c("/groups/MolBioService/CollaborationData/SEQUENCING RESULTS-new/2021/LOUB/2021-07-14_PCW-IMP-253_0231/",
#                       "/groups/MolBioService/CollaborationData/SEQUENCING RESULTS-new/2021/LOUB/2021-07-14_PCW-IMP-253_0232/"),
#                     full.names = T)
# files <- files[grepl(paste0(dat$seqID, collapse= "|"), files)]
# file.copy(files,
#           paste0("db/sanger_sequencing/actCP_DSCP_20210714/", basename(files)))

#---------------------------------------------------#
# IMPORT sequences
#---------------------------------------------------#
# source("/groups/stark/vloubiere/exp_data/update_files.R")
constructs <- fread("../../exp_data/vl_constructs_sequences.txt", 
                    key= "name")
pl <- fread("../../exp_data/vl_plasmids.txt", 
            key= "ID")

#---------------------------------------------------#
# Refseq
#---------------------------------------------------#
dat[, file:= list.files("db/sanger_sequencing/actCP_DSCP_20210714/", 
                        seqID, 
                        full.names = T), seqID]
dat[, refseq:= pl[ID=="pVL145", Sequence]]
dat[, refseq:= substr(refseq, 500, nchar(refseq)-250)]

#-------------------------------------------------------#
# Final check and plot
#-------------------------------------------------------#
pdf("pdf/design/sanger/sanger_actCP_DSCP_20210714.pdf")
# par(mar= c(5,10,5,5))
dat[, vl_sanger_align(refseq, 
                      abfiles = file, 
                      revcomp = as.logical(revComp), 
                      feat_sequences = constructs["actCP_2", sequence]), refseq]
dev.off()
