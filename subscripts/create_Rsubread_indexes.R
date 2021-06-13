# import
idx <- fread("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.txt")
idx <- idx[grepl("002|012|013|014", Sample_ID), .(Sample_ID= gsub("^([^_]*)_.*", "\\1", Sample_ID))]
idx <- unique(idx[, .(lib= gsub("^([^_]*)_.*", "\\1", Sample_ID))])

# Directories
idx[, dir:= paste0("/groups/stark/vloubiere/genomes/Rsubread_genome/dm3/custom_peSTARRSeq_", lib, "/")]
idx[, if(!dir.exists(dir)){dir.create(dir)}, dir]


# Import lib
lib <- as.data.table(readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/library/vl_library_112019.rds"))
lib <- lib[, .(ID_vl, oligo_full_sequence)]

# add template switch sequences
add <- fread("/groups/stark/vloubiere/exp_data/constructs_sequences.txt", key= "name")

# assemble libs
libs <- rbind(lib, 
              data.table(ID_vl= c("ts_SCR2_01002", "ts_HAM1_01003", "ts_SUP1_01004"), 
                         oligo_full_sequence= c(paste0(add[c("Flink_lib", "SCR2", "R1link_lib"), sequence], collapse= ""),
                                                paste0(add[c("Flink_lib", "HAM1", "R1link_lib"), sequence], collapse= ""),
                                                paste0(add[c("Flink_lib", "SUP1", "R1link_lib"), sequence], collapse= ""))))
# Assemble libs
idx <- idx[, {
  if(lib=="vllib002")
  {
    libs
  }else if(lib=="vllib012")
  {
    libs[grepl("_C_|ts_SCR2_01002", ID_vl)]
  }else  if(lib=="vllib013")
  {
    libs[grepl("_A_|_C_|ts_SCR2_01002", ID_vl)]
  }else  if(lib=="vllib014")
  {
    libs[grepl("_B_|_C_|ts_SCR2_01002", ID_vl)]
  }
}, idx]

# save fasta
idx[, fa:= paste0(dir, "peSTARRSeq_", lib, "_sequences.fa")]
idx[, {
  if(!file.exists(fa))
  {
    sequences <- DNAStringSet(oligo_full_sequence)
    names(sequences) <- ID_vl
    writeXStringSet(sequences, fa)
  }
  print("DONE")
}, fa]

# make index
idx[, basename:= paste0(dirname(fa), "/index/", gsub("sequences.fa$", "idx", basename(fa)))]
idx[, if(!dir.exists(dirname(basename))){dir.create(dirname(basename))}, basename]
idx[, buildindex(basename= basename, 
                 reference = fa), .(basename, fa)]
