setwd("/groups/stark/vloubiere/exp_data/")
system("sh update_files.sh")
setwd("/groups/stark/vloubiere/projects/pe_STARRSeq_2/")
fq <- fread("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.txt")
fq <- fq[grepl("002|012|013|014", Sample_ID)]
fq[, fq_prefix:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq_2/db/fastq/", `Flowcell ID`, "_", Sample_ID)]
fq[, i5:= as.character(reverseComplement(DNAStringSet(unlist(tstrsplit(i5, "-", keep=2)))))]
fq <- fq[, .(fq_file= paste0(fq_prefix, c("_1.fq.gz", "_2.fq.gz"))), fq]
fq <- fq[!file.exists(fq_file), .(cmd= extract_reads(BAM_path, fq_prefix, i5)), .(BAM_path, fq_prefix, i5)]
mclapply(fq$cmd, system, mc.preschedule = F, mc.cores = getDTthreads()-1)
