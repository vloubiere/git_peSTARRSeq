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
# IMPORTANT!!
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
meta[(DESeq2), FC_file_DESeq:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/FC_tables/DESeq2/", 
                                      group, "_counts_norm_final_oe.txt"), .(group, cdition, DESeq2)]
meta[(DESeq2), dds_file:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/FC_tables/DESeq2/", 
                                 group, ".dds"), .(group, cdition, DESeq2)]
meta[(DESeq2), FC_file_ratio:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/FC_tables/ratio/", 
                                      group, "_counts_norm_final_oe.txt"), .(group, cdition, DESeq2)]
fwrite(meta, "Rdata/metadata_processed.txt", na= NA)

#-------------------------------------------------------------#
# File removal?
#-------------------------------------------------------------#
pattern <- "vllib023|vllib024|vllib025|vllib026|vllib027|vllib028"
files <- list.files("db/fastq/", pattern, full.names = T)
files <- list.files("db/bam/", pattern, full.names = T)
files <- list.files("db/umi_counts/", pattern, full.names = T)
files <- list.files("db/merged_counts/", pattern, full.names = T)
files <- list.files("db/FC_tables/", pattern, full.names = T)

#-------------------------------------------------------------#
# PARALLELIZATION
#-------------------------------------------------------------#
cols <- c("fq1", "fq2", "bam", "bam_summary", "umi_counts", "pairs_counts", "FC_file_DESeq", "dds_file", "FC_file_ratio")
# meta <- meta[vllib=="vllib015" & DESeq2] # Example
meta[, check_exists:= all(file.exists(na.omit(unlist(.SD)))), .(group, DESeq2), .SDcols= cols]
meta <- meta[!(check_exists)]
meta[, {
  # Save as a .R script
  tmp <- tempfile(tmpdir = "/groups/stark/vloubiere/projects/pe_STARRSeq/logs/", fileext = ".txt")
  fwrite(cbind(.SD, DESeq2, group), tmp)
  # Bsub
  Rcmd <- paste("module load build-env/2020; module load r/3.6.2-foss-2018b; /software/2020/software/r/3.6.2-foss-2018b/bin/Rscript", 
                normalizePath("git_peSTARRSeq/functions/pipeline_3.0.R"),
                tmp)
  if(DESeq2)
  {
    bsub_cmd <- paste("/groups/stark/software-all/shell/bsub",
                      "-C 12", # N cpus
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

