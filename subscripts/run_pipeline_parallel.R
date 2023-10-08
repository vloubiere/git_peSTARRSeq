setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(readxl)
require(data.table)

# Update exp data (dropbox folder) and import metadata ----
if(F)
  source("/groups/stark/vloubiere/exp_data/update_files.R")
meta <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)
cols <- names(meta)
meta[, (cols):= lapply(.SD, function(x) ifelse(x=="NA", NA, x)), .SDcols= cols]

# Select used libraries ----
sel <- c("vllib002", # Large dev library
         "vllib006", # Reverse pSTARR-Seq
         "vllib015", # Restricted dev library
         "vllib016", # Restricted hk library
         "vllib025", # highHk CP
         "vllib026", # lowHk CP
         "vllib027", # lowDev CP
         "vllib028", # highDev CP
         "vllib029", # mutant lib
         "vllib030") # DHS lib
meta <- as.data.table(meta)[(DESeq2) & vllib %in% sel]

# Input and output filenames ---
meta[, fq1:= list.files("/scratch/stark/vloubiere/fastq/",
                        paste0(gsub(".bam$", "", basename(BAM_path)), ".*", i5, "_1.fq.gz"),
                        full.names = T), .(BAM_path, i5)]
meta[, fq2:= list.files("/scratch/stark/vloubiere/fastq/",
                        paste0(gsub(".bam$", "", basename(BAM_path)), ".*", i5, "_2.fq.gz"),
                        full.names = T), .(BAM_path, i5)]
meta[, index:= {
  fcase(vllib %in% c("vllib002", "vllib006"), "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist8_lib/twist8",
        vllib %in% c("vllib029", "vllib030"), "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist15_lib/twist15",
        default = "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist12_lib/twist12")
}, vllib]
meta[, bam:= paste0("/scratch/stark/vloubiere/bam/", vllib, "_", cdition, "_", DESeq2_pseudo_rep, ".bam")]
meta[, umi_counts:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/umi_counts/", vllib, "_", cdition, "_", DESeq2_pseudo_rep, ".txt")]

# Align and compute counts ----
meta[, cmd:= {
  # [required] 1/ The type of assay. Can be one of 'pe-STARR-Seq' or 'rev-pe-STARR-Seq' \n
  # [required] 2/ The index to use for the alignment \n
  # [required] 3/ A list of coma-separated _1.fq files containing read 1 \n
  # [required] 4/ A list of coma-separated _2.fq files containing read 2 \n
  # [required] 5/ Output bam file (.bam)\n
  # [required] 6/ Output .txt file where UMI-collapsed counts will be stored \n")
  if(!file.exists(umi_counts))
  {
    paste("module load build-env/2020; module load r/3.6.2-foss-2018b; Rscript git_peSTARRSeq/functions/pSTARRSeq_pipeline.R",
          type,
          index,
          paste0(unique(fq1), collapse = ","),
          paste0(unique(fq2), collapse = ","),
          bam,
          umi_counts)
  }else
    as.character(NA)
}, .(type, index, bam, umi_counts)]

# Submit ----
cores <- 8
mem <- 32
run <- meta[!is.na(cmd)]
if(nrow(run))
{
  run[, {
    vl_bsub(cmd, 
            cores= cores, 
            m = mem, 
            name = "vlloub", 
            t = '1-00:00:00',
            o= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/",
            e= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/")
  }, cmd]
}