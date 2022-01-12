setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(readxl)
require(data.table)

# PATHS
dir.create("/groups/stark/vloubiere/projects/pe_STARRSeq/db/fastq/", showWarnings = F)
dir.create("/groups/stark/vloubiere/projects/pe_STARRSeq/db/sam/", showWarnings = F)
dir.create("/groups/stark/vloubiere/projects/pe_STARRSeq/db/umi_counts/", showWarnings = F)
dir.create("/groups/stark/vloubiere/projects/pe_STARRSeq/db/merged_counts/", showWarnings = F)
dir.create("/groups/stark/vloubiere/projects/pe_STARRSeq/db/dds/", showWarnings = F)
dir.create("/groups/stark/vloubiere/projects/pe_STARRSeq/db/FC_tables/", showWarnings = F)
dir.create("/groups/stark/vloubiere/projects/pe_STARRSeq/db/final_tables_exp_model/", showWarnings = F)
dir.create("/groups/stark/vloubiere/projects/pe_STARRSeq/db/final_tables_exp_model/counts_norm/", showWarnings = F)
dir.create("/groups/stark/vloubiere/projects/pe_STARRSeq/db/final_tables_exp_model/replicates_counts_norm/", showWarnings = F)

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
meta[, sam:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/sam/", output_prefix, ".sam")]
meta[, sam_summary:= gsub(".sam$", ".sam.summary", sam)]
meta[, umi_counts:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/umi_counts/", output_prefix, ".txt")]
meta[, umi_summary:= gsub(".txt$", "_summary.txt", umi_counts)]
meta[, c("summary_counts", "pairs_counts", "spike_counts", "switched_counts"):= {
  dir <- paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/merged_counts/", group)
  if(DESeq2)
    dir <- paste0(dir, "_", cdition, "_rep", DESeq2_pseudo_rep)
  .(paste0(dir, "_summary.txt"),
    paste0(dir, "_merged_pair_counts.txt"),
    paste0(dir, "_merged_spikein_counts.txt"),
    paste0(dir, "_merged_switched_counts.txt"))
}, .(group, cdition, DESeq2, DESeq2_pseudo_rep)]
meta[(DESeq2), FC_file:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/final_tables_exp_model/counts_norm/", 
                                group, "_counts_norm_final_oe.txt"), .(group, cdition, DESeq2)]
fwrite(meta, "Rdata/metadata_processed.txt", na= NA)

#-------------------------------------------------------------#
# File removal?
#-------------------------------------------------------------#
pattern <- "vllib015|vllib016|vllib017|vllib018|vllib019|vllib020|vllib021|vllib022"
files <- list.files("db/fastq/", pattern, full.names = T)
files <- c(files, list.files("db/sam/", pattern, full.names = T))
files <- c(files, list.files("db/umi_counts/", pattern, full.names = T))
files <- c(files, list.files("db/merged_counts/", pattern, full.names = T))
files <- c(files, list.files("db/final_tables_exp_model/counts_norm/", pattern, full.names = T))

#-------------------------------------------------------------#
# PARALLELIZATION
#-------------------------------------------------------------#
cols <- c("fq1", "fq2", "sam", "sam_summary", "umi_counts", "umi_summary", 
          "summary_counts", "pairs_counts", "spike_counts", "switched_counts", "FC_file")
meta$check_exists <- meta[, apply(.SD, 1, function(x) all(file.exists(x[!is.na(x)]))), .SDcols= cols]
meta[, {
  if(any(!check_exists))
  {
    # Save as a .R script
    test <<- tmp <- tempfile(tmpdir = "/groups/stark/vloubiere/projects/pe_STARRSeq/logs/", fileext = ".txt")
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
                        "-m 8", # memory
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
  }
  print("Submitted!")
}, .(group, DESeq2)]

