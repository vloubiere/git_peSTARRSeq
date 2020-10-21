setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
require(data.table)
require(Biostrings)

# Load metadata
dat <- fread("db/metadata/vl_sequencing_metadata.txt")
dat[, c("Sample_ID", "condition"):= .(paste0(`Flowcell ID`, "_", Sample_ID), Sample_ID)]
dat <- dat[grepl("libvl002", Sample_ID) & (SEQUENCING_WORKED)]

#-------------------------------------------------#
# Process sequencing data
#-------------------------------------------------#
# Extract reads from seq bam -> fastq
dat[, fq_prefix:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/fastq/", Sample_ID)]
dat[, i5:= as.character(reverseComplement(DNAStringSet(unlist(tstrsplit(i5, "-", keep=2)))))]
dat <- dat[, .(fq_file= paste0(fq_prefix, c("_1.fq.gz", "_2.fq.gz"))), c(colnames(dat))]
dat[!file.exists(fq_file), extract_cmd:= extract_reads(BAM_path, fq_prefix, i5), .(BAM_path, fq_prefix, i5)]
# Trim reads1 (3 bp)
dat[grepl("_1.fq.gz$", fq_file), trimmed_fq1_file := gsub("_1.fq.gz", "_1.33bp_3prime.fq.gz", fq_file)]
dat[ifelse(is.na(trimmed_fq1_file), F, !file.exists(trimmed_fq1_file)), trim_fq1_cmd:= trim3_fastq(fq_file, keep= 33, "/groups/stark/vloubiere/projects/pe_STARRSeq/db/fastq/", cores = 8), fq_file]
# Alignment
dat[, sam_file := paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/sam/", Sample_ID, ".sam")]
dat[, bwt_idx:= "/groups/stark/vloubiere/genomes/Custom_peSTARRSeq_1/index_bowtie1/custom_peSTARR1"]
dat[!file.exists(sam_file), align_cmd:= bowtie1_align(fq1= paste0(trimmed_fq1_file[1], ",", fq_file[2]), index_prefix= bwt_idx, sam_output= sam_file, cores= 12), Sample_ID]
dat[!file.exists(sam_file), align_cmd:= paste0(align_cmd, " --sam-nohead")]
# Process sam files using custom functions
dat[, rdsall_file:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/read_counts/", condition, ".all.rds"), condition]
dat[, rdsuniq_file:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/read_counts/", condition, ".uniq.UMI.rds"), condition]
dat[!file.exists(rdsall_file) | !file.exists(rdsuniq_file), 
    UMI_collapsing_cmd:= Rsub("/groups/stark/vloubiere/pipelines/peSTARRSeq_sam_counts_UMI_collapsing.R", c(unique(sam_file), rdsall_file, rdsuniq_file)), 
    .(condition, rdsall_file, rdsuniq_file)]

# RUN
run <- melt(dat, measure= patterns(cmd= "cmd$"))
run[!is.na(value), bsub(paste0(unique(value), collapse= ";"), cores = 12, o= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/", name = paste0(Sample_ID, "_vl")), Sample_ID]

