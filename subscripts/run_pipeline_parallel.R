setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(readxl)
require(data.table)

# Import metadata with fq files
meta <- readRDS("Rdata/metadata.rds")

# Output filenames ---
meta[, index:= {
  fcase(library=="T8", "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist8_lib/twist8", # WT oligo pool
        library=="T12", "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist12_lib/twist12", # Focused oligo pool
        vllib=="vllib029", "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist15_mutant_lib/twist15", # Mutated oligo pool
        vllib=="vllib030", "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist15_DHS_lib/twist15") # DHS oligo pool
}]
meta[, Trim3:= 0]
meta[, Trim5:= 0]
meta[, nMismMatch:= ifelse(vllib=="vllib029", 5, 3)]
meta[, bam:= paste0("/scratch/stark/vloubiere/bam/", screen, "_", cdition, "_", rep, ".bam")]
meta[, umi_counts:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/umi_counts/", screen, "_", cdition, "_", rep, ".txt")]

# Save processed metadata ----
saveRDS(meta,
        "Rdata/metadata_processed.rds")

# Align and compute counts ----
meta[, cmd:= {
  # [required] 1/ The type of assay. Can be one of 'pe-STARR-Seq' or 'rev-pe-STARR-Seq' \n
  # [required] 2/ The index to use for the alignment \n
  # [required] 3/ A list of coma-separated _1.fq files containing read 1 \n
  # [required] 4/ A list of coma-separated _2.fq files containing read 2 \n
  # [required] 5/ Output bam file (.bam)\n
  # [required] 6/ Output .txt file where UMI-collapsed counts will be stored \n
  # [required] 7/ Number of allowed mismatches \n
  # [required] 8/ Trim5 \n
  # [required] 9/ Trim3 \n")
  if(!file.exists(umi_counts))
  {
    paste("/software/f2022/software/r/4.3.0-foss-2022b/bin/Rscript /groups/stark/vloubiere/projects/pe_STARRSeq/git_peSTARRSeq/functions/pSTARRSeq_pipeline.R",
          type,
          index,
          paste0(unique(fq1), collapse = ","),
          paste0(unique(fq2), collapse = ","),
          bam,
          umi_counts,
          nMismMatch,
          Trim5,
          Trim3)
  }else
    as.character(NA)
}, .(type, index, bam, umi_counts, nMismMatch, Trim5, Trim3)]

# Submit ----
cores <- 12
mem <- 64
run <- meta[!is.na(cmd)]
if(nrow(run))
{
  run[, {
    vl_bsub(cmd,
            cores= cores,
            m = mem,
            name = "vlloub", 
            t = '3-00:00:00',
            o= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/",
            e= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/logs/")
  }, cmd]
}