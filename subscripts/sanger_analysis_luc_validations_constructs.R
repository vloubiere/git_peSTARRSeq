dat <- readRDS("Rdata/validations_luciferase_constructs_metadata.rds")

feat <- fread("/groups/stark/vloubiere/exp_data/constructs_sequences.txt", key = "name")
feat <- feat[c("Flink_+0", "R1link+0", "R2link+0", "R3link+0", "SCR1")]

pdf("pdf/sanger_luciferase_constructs.pdf", height = 2, width = 10)
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
