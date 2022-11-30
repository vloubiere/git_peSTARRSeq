setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(readxl)
require(data.table)

#--------------------------------------------------------------#
# Update exp data
# Fetch dropbox folder containing my metadata and update local files
#--------------------------------------------------------------#
if(F)
  source("/groups/stark/vloubiere/exp_data/update_files.R")
meta <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)
cols <- colnames(meta)
meta[, (cols):= lapply(.SD, function(x) ifelse(x=="NA", NA, x)), .SDcols= cols]
meta[, output_prefix:= paste0("/", my_ID, "__", gsub(".bam$", "", basename(BAM_path)))]
meta <- meta[(DESeq2)]
# Check unique prefixes!!
if(any(meta[, .N, output_prefix]$N>1))
  stop(paste0("Some output prefixes are not unique and cannot be used. Check metadata table (replicates?)!\n"))

#-------------------------------------------------------------#
# Output names
#-------------------------------------------------------------#
meta[, fq_prefix:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/fastq/", output_prefix)]
meta[, fq1:= paste0(fq_prefix, "_1.fq.gz")]
meta[, fq2:= paste0(fq_prefix, "_2.fq.gz")]
meta[, bam:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/bam/", output_prefix, ".bam")]
meta[, bam_summary:= gsub(".bam$", ".bam.summary", bam)]
meta[, umi_counts:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/umi_counts/", output_prefix, ".txt.gz")]
meta[, pairs_counts:= {
  dir <- paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/merged_counts/", group)
  if(DESeq2)
    dir <- paste0(dir, "_", cdition, "_rep", DESeq2_pseudo_rep)
  paste0(dir, "_merged_pair_counts.txt")
}, .(group, cdition, DESeq2, DESeq2_pseudo_rep)]
meta[(DESeq2), dds_file:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/dds/", group, ".dds")]
meta[(DESeq2), FC_file_DESeq2:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/FC_tables_DESeq2/", group, "_DESeq2_final_oe.rds")]
meta[(DESeq2), FC_file_ratio:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/FC_tables_ratio/", group, "_counts_norm_final_oe.rds")]
fwrite(meta, "Rdata/metadata_processed.txt", na= NA)

#-------------------------------------------------------------#
# PARALLELIZATION
#-------------------------------------------------------------#
cols <- c("fq1", "fq2", "bam", "bam_summary", "umi_counts", "pairs_counts", "dds_file", "FC_file_DESeq2", "FC_file_ratio")
meta[, check_exists:= all(file.exists(na.omit(unlist(.SD)))), .(group, DESeq2), .SDcols= cols]
meta <- meta[!(check_exists)]
meta[, {
  # Save .R script
  tmp <- tempfile(tmpdir = "/groups/stark/vloubiere/projects/pe_STARRSeq/logs/", fileext = ".txt")
  fwrite(cbind(.SD, DESeq2, group), tmp)
  # Bsub
  Rcmd <- paste("module load r/4.1.2-foss-2021b; Rscript", 
                normalizePath("git_peSTARRSeq/subscripts/pipeline_3.0.R"),
                tmp)
  if(DESeq2)
  {
    bsub_cmd <- paste("/groups/stark/software-all/shell/bsub",
                      "-C 8", # N cpus
                      "-m 32", # memory
                      paste("-n", group), #name
                      "-T '2-00:00:00'", #name
                      "-o /groups/stark/vloubiere/projects/pe_STARRSeq/logs/", #stdo
                      "-e /groups/stark/vloubiere/projects/pe_STARRSeq/logs/") #stde 
    
  }else{
    bsub_cmd <- paste("/groups/stark/software-all/shell/bsub",
                      "-C 4", # N cpus
                      "-m 16", # memory
                      paste("-n", group), #name
                      "-T '08:00:00'", #name
                      "-o /groups/stark/vloubiere/projects/pe_STARRSeq/logs/", #stdo
                      "-e /groups/stark/vloubiere/projects/pe_STARRSeq/logs/") #stde
  }
  # Wrap and Submit
  bsub_cmd <- paste0(bsub_cmd, " \"", Rcmd, "\"")
  Sys.unsetenv("SBATCH_RESERVATION")
  Sys.unsetenv("SBATCH_WCKEY")
  job_ID <- system(bsub_cmd, intern = T)
  # Return bsub ID
  unlist(job_ID[2])
  print("Submitted!")
}, .(group, DESeq2)]