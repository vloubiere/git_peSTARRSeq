setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(readxl)

# Assemble metadata files
dat <- as.data.table(read_excel("old_versions/20210221_82cdc9a_backup/Rdata/luciferase_validations/clean_stocks.xlsx"))
colnames(dat)[colnames(dat)=="enh_L"] <- "L"
colnames(dat)[colnames(dat)=="enh_R"] <- "R"
dat <- dat[, .(seqID= grep("NA", unlist(c(seq_ID_1, seq_ID_2)), invert = T, value = T), Comments), Sample_ID:R]
dat[, Sequencing_date:= tstrsplit(seqID, "_", keep = 2)]
dat[, c("sanger_L", "sanger_R"):= {
  folder <- grep(Sequencing_date, list.dirs("db/sanger_sequencing/luc_validations/"), value= T)
  .(list.files(folder, paste0("LOUB_", Sample_ID, "_", ".*CASeq001"), full.names = T),
    list.files(folder, paste0("LOUB_", Sample_ID, "_", ".*LOH011"), full.names = T))
}, .(Sample_ID, Sequencing_date)]
dat[, c("Sample_ID", "Colony_ID"):= .(unlist(tstrsplit(Sample_ID, "c", keep= 1)), Sample_ID)]
dat$seqID <- NULL

# Generate refseq
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
constructs <- fread("/groups/stark/vloubiere/exp_data/vl_constructs_sequences.txt", key= "name")
backbone <- vl_digest(constructs["DSCP_STARRSeq", sequence], c("BshTI", "SalI"))
dat[, refseq:= paste0(backbone[1], 
                      lib[.BY, oligo_full_sequence, on= "ID_vl==L"], 
                      constructs["CGCov_F", sequence],
                      constructs["SCR1", sequence],
                      constructs["CGCov_R", sequence],
                      lib[.BY, oligo_full_sequence, on= "ID_vl==R"], 
                      backbone[3]), .(L, R)] 
dat[, refseq:= toupper(refseq)]

# Import other sequences
feat <- fread("/groups/stark/vloubiere/exp_data/vl_constructs_sequences.txt", key = "name")
feat <- feat[c("Flink_+0", "R1link+0", "R2link+0", "R3link+0", "SCR1")]

pdf("pdf/sanger_sequencing/peSTARRSeq_SCR1_luciferase_validations_constructs.pdf", height = 2, width = 10)
par(mar= c(1, 15, 5, 10))
for(i in seq(nrow(dat)))
{
  files <- c(dat[i, sanger_L], dat[i, sanger_R])
  if(!any(is.na(files)))
  {
    vl_sanger_align(refseq = dat[i, refseq], 
                    revcomp = c(F, T),
                    abfiles = c(dat[i, sanger_L], dat[i, sanger_R]), 
                    feat_sequences = feat$sequence, 
                    feat_names = feat$name,
                    feat_cols = c("red", "orange", "yellow", "green", "cornflowerblue"))
    mtext(paste(dat[i, L], "x", dat[i, R]), line = 3)
    mtext(paste("sample:", dat[i, Sample_ID], "| colony:", dat[i, Colony_ID]), line = 1)
  }
}
dev.off()
