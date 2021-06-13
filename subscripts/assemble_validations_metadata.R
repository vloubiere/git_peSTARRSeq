# Assemble metadata files
res <- as.data.table(read_excel("../pe_STARRSeq/Rdata/luciferase_validations/clean_stocks.xlsx"))
colnames(res)[colnames(res)=="enh_L"] <- "L"
colnames(res)[colnames(res)=="enh_R"] <- "R"
res <- res[, .(seqID= grep("NA", unlist(c(seq_ID_1, seq_ID_2)), invert = T, value = T), Comments), Sample_ID:R]
res[, Sequencing_date:= tstrsplit(seqID, "_", keep = 2)]
res[, c("sanger_L", "sanger_R"):= {
   folder <- grep(Sequencing_date, list.dirs("db/sanger"), value= T)
   .(list.files(folder, paste0("LOUB_", Sample_ID, "_", ".*CASeq001"), full.names = T),
     list.files(folder, paste0("LOUB_", Sample_ID, "_", ".*LOH011"), full.names = T))
}, .(Sample_ID, Sequencing_date)]
res[, c("Sample_ID", "Colony_ID"):= .(unlist(tstrsplit(Sample_ID, "c", keep= 1)), Sample_ID)]
res$seqID <- NULL

# Generate refseq
lib <- as.data.table(readRDS("Rdata/vl_library_112019.rds"))
constructs <- fread("/groups/stark/vloubiere/exp_data/constructs_sequences.txt", key= "name")
res[, refseq:= paste0(constructs["DSCP_upstream_OLDcloning_site", sequence], 
                      lib[.BY, oligo_full_sequence, on= "ID_vl==L"], 
                      constructs["CGCov_F", sequence],
                      constructs["SCR1", sequence],
                      constructs["CGCov_R", sequence],
                      lib[.BY, oligo_full_sequence, on= "ID_vl==R"], 
                      constructs["DSCP_downstream_OLDcloning_site", sequence]), .(L, R)] 
res[, refseq:= toupper(refseq)]

saveRDS(res, "Rdata/validations_luciferase_constructs_metadata.rds")
