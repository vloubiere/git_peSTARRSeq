setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(readxl)
require(data.table)

read_meta <- function(x)
{
  meta <- read_xlsx(x)
  meta <- as.data.table(meta)
  meta <- as.data.table(lapply(meta, function(x) ifelse(x=="NA", NA, x)))
  meta[, output_prefix:= paste0(my_ID, "__", gsub(".bam$", "", basename(BAM_path)))]
  return(meta)
}

# source("/groups/stark/vloubiere/exp_data/update_files.R")
new <- read_meta("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
old <- read_meta("/groups/stark/vloubiere/exp_data/old/vl_sequencing_metadata_20211104.xlsx")
all <- merge(new,
             old,
             suffixes= c("_new", "_old"),
             by= c("Flowcell ID", 
                   "Lane number", 
                   "Submission_date", 
                   "Sample_CSF_ID", 
                   "Request ID",
                   "Sample_ID",
                   "Replicate", 
                   "i5",
                   "idx",
                   "Owner",
                   "Method",
                   "Pool",
                   "Reads_needed_mln",
                   "Genome",
                   "Spike-ins_Genome",
                   "Source",
                   "Fragment Size",
                   "Read_Type",
                   "Platform",
                   "Sequencing_date",
                   "BAM_path"))

# PATHS
dir_fq <- normalizePath("db/fastq/", mustWork = F)
dir_sam <- normalizePath("db/sam/", mustWork = F)
dir_umiCounts <- normalizePath("db/umi_counts/", mustWork = F)
dir_mergedCounts <- normalizePath("db/merged_counts/", mustWork = F)
dir_dds <- normalizePath("db/dds/", mustWork = F)
dir_FC <- normalizePath("db/FC_tables/", mustWork = F)
dir_final <- normalizePath("db/final_tables_exp_model/", mustWork = F)

# Rename function
rename_pipe <- function(dir, old_pattern, new_pattern)
{
  .c <- list.files(dir, 
                   pattern = old_pattern, 
                   full.names = T)
  if(length(.c))
    file.rename(.c, gsub(old_pattern, new_pattern, .c))
}
  

#-------------------#
# RENAME
#-------------------#
all[, {
  rename_pipe(dir_fq, output_prefix_old, output_prefix_new) # fastq
  rename_pipe(dir_sam, output_prefix_old, output_prefix_new) # sam
  rename_pipe(dir_umiCounts, output_prefix_old, output_prefix_new) # umi counts
  # merged counts
  if(!is.na(DESeq2_group_old))
  {
    rename_pipe(dir_mergedCounts, DESeq2_group_old, DESeq2_group_new) # Merged counts
    rename_pipe(paste0(dir_final, "/counts_norm"), DESeq2_group_old, DESeq2_group_new) # Counts normalized exp/obs FC
    rename_pipe(dir_dds, DESeq2_group_old, DESeq2_group_new) # DESeq2 dds
    rename_pipe(dir_FC, DESeq2_group_old, DESeq2_group_new) # DESeq2 FC
    rename_pipe(paste0(dir_final, "/DESeq2"), DESeq2_group_old, DESeq2_group_new) # DESeq2 normalized exp/obs FC
  }
}, .(output_prefix_old, output_prefix_new, DESeq2_group_old, DESeq2_group_new)]





